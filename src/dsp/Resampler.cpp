/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*- vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band Library
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2018 Particular Programs Ltd.

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.

    Alternatively, if you have a valid commercial licence for the
    Rubber Band Library obtained by agreement with the copyright
    holders, you may redistribute and/or modify it under the terms
    described in that licence.

    If you wish to distribute code using the Rubber Band Library
    under terms other than those of the GNU General Public License,
    you must obtain a valid commercial licence before doing so.
*/

#include "Resampler.h"
#include "base/Profiler.h"

#include <cstdlib>
#include <cmath>

#include <iostream>
#include <algorithm>

#include "system/Allocators.h"

#ifdef HAVE_IPP
#include <ippversion.h>
#if (IPP_VERSION_MAJOR < 7)
#include <ipps.h>
#include <ippsr.h>
#include <ippac.h>
#else
#include <ipps.h>
#endif
#endif

#ifdef HAVE_LIBSAMPLERATE
#include <samplerate.h>
#endif

#ifdef HAVE_LIBRESAMPLE
#include <libresample.h>
#endif

#ifdef USE_SPEEX
#include "speex/speex_resampler.h"
#endif

#ifndef HAVE_IPP
#ifndef HAVE_LIBSAMPLERATE
#ifndef HAVE_LIBRESAMPLE
#ifndef USE_SPEEX
#error No resampler implementation selected!
#endif
#endif
#endif
#endif

namespace RubberBand {

class ResamplerImpl
{
public:
    virtual ~ResamplerImpl() { }
    
    virtual int resample(const float *const R__ *const R__ in, 
                         float *const R__ *const R__ out,
                         int incount,
                         float ratio,
                         bool final) = 0;
    
    virtual int resampleInterleaved(const float *const R__ in, 
                                    float *const R__ out,
                                    int incount,
                                    float ratio,
                                    bool final) = 0;

    virtual int getChannelCount() const = 0;

    virtual void reset() = 0;
};

namespace Resamplers {

#ifdef HAVE_IPP

class D_IPP : public ResamplerImpl
{
public:
    D_IPP(Resampler::Quality quality, int channels, int maxBufferSize,
          int debugLevel);
    ~D_IPP();

    int resample(const float *const R__ *const R__ in,
                 float *const R__ *const R__ out,
                 int incount,
                 float ratio,
                 bool final);

    int resampleInterleaved(const float *const R__ in,
                            float *const R__ out,
                            int incount,
                            float ratio,
                            bool final = false);

    int getChannelCount() const { return m_channels; }

    void reset();

protected:
    // to m_outbuf
    int doResample(int outcount, double ratio, bool final);

    IppsResamplingPolyphase_32f **m_state;
    float **m_inbuf;
    size_t m_inbufsz;
    float **m_outbuf;
    size_t m_outbufsz;
    int m_bufsize;
    int m_channels;
    int m_window;
    float m_factor;
    int m_history;
    int *m_lastread;
    double *m_time;
    int m_debugLevel;
    
    void setBufSize(int);
};

D_IPP::D_IPP(Resampler::Quality quality, int channels, int maxBufferSize,
             int debugLevel) :
    m_state(0),
    m_channels(channels),
    m_debugLevel(debugLevel)
{
    if (m_debugLevel > 0) {
        std::cerr << "Resampler::Resampler: using IPP implementation"
                  << std::endl;
    }

    int nStep = 16;
    IppHintAlgorithm hint = ippAlgHintFast;

    switch (quality) {

    case Resampler::Best:
        m_window = 64;
        nStep = 80;
        hint = ippAlgHintAccurate;
        break;

    case Resampler::FastestTolerable:
        nStep = 16;
        m_window = 16;
        hint = ippAlgHintFast;
        break;

    case Resampler::Fastest:
        m_window = 24;
        nStep = 64;
        hint = ippAlgHintFast;
        break;
    }

    m_factor = 8; // initial upper bound on m_ratio, may be amended later
    m_history = int(m_window * 0.5 * std::max(1.0, 1.0 / m_factor)) + 1;

    m_state = new IppsResamplingPolyphase_32f *[m_channels];

    m_lastread = new int[m_channels];
    m_time = new double[m_channels];

    m_inbufsz = 0;
    m_outbufsz = 0;
    m_inbuf = 0;
    m_outbuf = 0;
    m_bufsize = 0;

    setBufSize(maxBufferSize + m_history);

    if (m_debugLevel > 1) {
        std::cerr << "bufsize = " << m_bufsize << ", window = " << m_window << ", nStep = " << nStep << ", history = " << m_history << std::endl;
    }

#if (IPP_VERSION_MAJOR >= 7)
    int specSize = 0;
    ippsResamplePolyphaseGetSize_32f(float(m_window),
                                     nStep,
                                     &specSize,
                                     hint);
    if (specSize == 0) {
#ifndef NO_EXCEPTIONS
        throw Resampler::ImplementationError;
#else        
        abort();
#endif
    }
#endif

    for (int c = 0; c < m_channels; ++c) {
#if (IPP_VERSION_MAJOR < 7)
        ippsResamplePolyphaseInitAlloc_32f(&m_state[c],
                                           float(m_window),
                                           nStep,
                                           0.95f,
                                           9.0f,
                                           hint);
#else
        m_state[c] = (IppsResamplingPolyphase_32f *)ippsMalloc_8u(specSize);
        ippsResamplePolyphaseInit_32f(float(m_window),
                                      nStep,
                                      0.95f,
                                      9.0f,
                                      m_state[c],
                                      hint);
#endif
        
        m_lastread[c] = m_history;
        m_time[c] = m_history;
    }

    if (m_debugLevel > 1) {
        std::cerr << "Resampler init done" << std::endl;
    }
}

D_IPP::~D_IPP()
{
#if (IPP_VERSION_MAJOR < 7)
    for (int c = 0; c < m_channels; ++c) {
        ippsResamplePolyphaseFree_32f(m_state[c]);
    }
#else
    for (int c = 0; c < m_channels; ++c) {
        ippsFree(m_state[c]);
    }
#endif

    deallocate_channels(m_inbuf, m_channels);
    deallocate_channels(m_outbuf, m_channels);

    delete[] m_lastread;
    delete[] m_time;
    delete[] m_state;
}

void
D_IPP::setBufSize(int sz)
{
    if (m_debugLevel > 1) {
        std::cerr << "resize bufsize " << m_bufsize << " -> ";
    }

    m_bufsize = sz;

    if (m_debugLevel > 1) {
        std::cerr << m_bufsize << std::endl;
    }

    int n1 = m_bufsize + m_history + 2;
    int n2 = lrintf(ceil((m_bufsize - m_history) * m_factor + 2));

    if (m_debugLevel > 1) {
        std::cerr << "(outbufsize = " << n2 << ")" << std::endl;
    }

    m_inbuf = reallocate_and_zero_extend_channels
        (m_inbuf, m_channels, m_inbufsz, m_channels, n1);

    m_outbuf = reallocate_and_zero_extend_channels
        (m_outbuf, m_channels, m_outbufsz, m_channels, n2);
            
    m_inbufsz = n1;
    m_outbufsz = n2;
}

int
D_IPP::resample(const float *const R__ *const R__ in,
                float *const R__ *const R__ out,
                int incount,
                float ratio,
                bool final)
{
    if (ratio > m_factor) {
        m_factor = ratio;
        m_history = int(m_window * 0.5 * std::max(1.0, 1.0 / m_factor)) + 1;
    }

    for (int c = 0; c < m_channels; ++c) {
        if (m_lastread[c] + incount + m_history > m_bufsize) {
            setBufSize(m_lastread[c] + incount + m_history);
        }
    }

    for (int c = 0; c < m_channels; ++c) {
        for (int i = 0; i < incount; ++i) {
            m_inbuf[c][m_lastread[c] + i] = in[c][i];
        }
        m_lastread[c] += incount;
    }

    int got = doResample(int(round(incount * ratio)), ratio, final);
    
    for (int c = 0; c < m_channels; ++c) {
        v_copy(out[c], m_outbuf[c], got);
    }

    return got;
}

int
D_IPP::resampleInterleaved(const float *const R__ in,
                           float *const R__ out,
                           int incount,
                           float ratio,
                           bool final)
{
    if (ratio > m_factor) {
        m_factor = ratio;
        m_history = int(m_window * 0.5 * std::max(1.0, 1.0 / m_factor)) + 1;
    }

    for (int c = 0; c < m_channels; ++c) {
        if (m_lastread[c] + incount + m_history > m_bufsize) {
            setBufSize(m_lastread[c] + incount + m_history);
        }
    }

    for (int c = 0; c < m_channels; ++c) {
        for (int i = 0; i < incount; ++i) {
            m_inbuf[c][m_lastread[c] + i] = in[i * m_channels + c];
        }
        m_lastread[c] += incount;
    }

    int got = doResample(int(round(incount * ratio)), ratio, final);

    v_interleave(out, m_outbuf, m_channels, got);

    return got;
}

int
D_IPP::doResample(int outspace, double ratio, bool final)
{
    int outcount = 0;
    
    for (int c = 0; c < m_channels; ++c) {

        int n = m_lastread[c] - m_history - int(m_time[c]);

        if (c == 0 && m_debugLevel > 2) {
            std::cerr << "at start, lastread = " << m_lastread[c] << ", history = "
                 << m_history << ", time = " << m_time[c] << ", therefore n = "
                 << n << std::endl;
        }

        if (n <= 0) {
            if (c == 0 && m_debugLevel > 1) {
                std::cerr << "not enough input samples to do anything" << std::endl;
            }
            continue;
        }
        
        if (c == 0 && m_debugLevel > 2) {
            std::cerr << "before resample call, time = " << m_time[c] << std::endl;
        }

        // We're committed to not overrunning outspace, so we need to
        // offer the resampler only enough samples to ensure it won't

        int limit = int(floor(outspace / ratio));
        if (n > limit) {
            if (c == 0 && m_debugLevel > 1) {
                std::cerr << "trimming input samples from " << n << " to " << limit
                     << " to avoid overrunning " << outspace << " at output"
                     << std::endl;
            }
            n = limit;
        }
        
#if (IPP_VERSION_MAJOR < 7)
        ippsResamplePolyphase_32f(m_state[c],
                                  m_inbuf[c],
                                  n,
                                  m_outbuf[c],
                                  ratio,
                                  1.0f,
                                  &m_time[c],
                                  &outcount);
#else
        ippsResamplePolyphase_32f(m_inbuf[c],
                                  n,
                                  m_outbuf[c],
                                  ratio,
                                  1.0f,
                                  &m_time[c],
                                  &outcount,
                                  m_state[c]);
#endif

        int t = int(round(m_time[c]));
        
        if (c == 0 && m_debugLevel > 2) {
            std::cerr << "converted " << n << " samples to " << outcount
                 << ", time advanced to " << t << std::endl;
            std::cerr << "will move " << m_lastread[c] + m_history - t
                 << " unconverted samples back from index " << t - m_history
                 << " to 0" << std::endl;
        }

        v_move(m_inbuf[c],
               m_inbuf[c] + t - m_history,
               m_lastread[c] + m_history - t);

        m_lastread[c] -= t - m_history;
        m_time[c] -= t - m_history;

        if (c == 0 && m_debugLevel > 2) {
            std::cerr << "lastread reduced to " << m_lastread[c]
                 << ", time reduced to " << m_time[c]
                 << std::endl;
        }
        
        if (final && n < limit) {

            // Looks like this actually produces too many samples
            // (additionalcount is a few samples too large).

            // Also, we aren't likely to have enough space in the
            // output buffer as the caller won't have allowed for
            // all the samples we're retrieving here.

            // What to do?

            int additionalcount = 0;

            if (c == 0 && m_debugLevel > 2) {
                std::cerr << "final call, padding input with " << m_history
                     << " zeros (symmetrical with m_history)" << std::endl;
            }
            
            for (int i = 0; i < m_history; ++i) {
                m_inbuf[c][m_lastread[c] + i] = 0.f;
            }

            if (c == 0 && m_debugLevel > 2) {
                std::cerr << "before resample call, time = " << m_time[c] << std::endl;
            }

            int nAdditional = m_lastread[c] - int(m_time[c]);

            if (n + nAdditional > limit) {
                if (c == 0 && m_debugLevel > 1) {
                    std::cerr << "trimming final input samples from " << nAdditional
                         << " to " << (limit - n)
                         << " to avoid overrunning " << outspace << " at output"
                         << std::endl;
                }
                nAdditional = limit - n;
            }
            
#if (IPP_VERSION_MAJOR < 7)
            ippsResamplePolyphase_32f(m_state[c],
                                      m_inbuf[c],
                                      nAdditional,
                                      m_outbuf[c],
                                      ratio,
                                      1.0f,
                                      &m_time[c],
                                      &additionalcount);
#else
            ippsResamplePolyphase_32f(m_inbuf[c],
                                      nAdditional,
                                      m_outbuf[c],
                                      ratio,
                                      1.0f,
                                      &m_time[c],
                                      &additionalcount,
                                      m_state[c]);
#endif
        
            if (c == 0 && m_debugLevel > 2) {
                std::cerr << "converted " << n << " samples to " << additionalcount
                     << ", time advanced to " << m_time[c] << std::endl;
                std::cerr << "outcount = " << outcount << ", additionalcount = " << additionalcount << ", sum " << outcount + additionalcount << std::endl;
            }

            if (c == 0) {
                outcount += additionalcount;
            }
        }
    }

    if (m_debugLevel > 2) {
        std::cerr << "returning " << outcount << " samples" << std::endl;
    }
    
    return outcount;
}

void
D_IPP::reset()
{
    //!!!
}

#endif /* HAVE_IPP */

#ifdef HAVE_LIBSAMPLERATE

class D_SRC : public ResamplerImpl
{
public:
    D_SRC(Resampler::Quality quality, int channels, int maxBufferSize,
          int m_debugLevel);
    ~D_SRC();

    int resample(const float *const R__ *const R__ in,
                 float *const R__ *const R__ out,
                 int incount,
                 float ratio,
                 bool final);

    int resampleInterleaved(const float *const R__ in,
                            float *const R__ out,
                            int incount,
                            float ratio,
                            bool final = false);

    int getChannelCount() const { return m_channels; }

    void reset();

protected:
    SRC_STATE *m_src;
    float *m_iin;
    float *m_iout;
    float m_lastRatio;
    int m_channels;
    int m_iinsize;
    int m_ioutsize;
    int m_debugLevel;
};

D_SRC::D_SRC(Resampler::Quality quality, int channels, int maxBufferSize,
             int debugLevel) :
    m_src(0),
    m_iin(0),
    m_iout(0),
    m_lastRatio(1.f),
    m_channels(channels),
    m_iinsize(0),
    m_ioutsize(0),
    m_debugLevel(debugLevel)
{
    if (m_debugLevel > 0) {
        std::cerr << "Resampler::Resampler: using libsamplerate implementation"
                  << std::endl;
    }

    int err = 0;
    m_src = src_new(quality == Resampler::Best ? SRC_SINC_BEST_QUALITY :
                    quality == Resampler::Fastest ? SRC_LINEAR :
                    SRC_SINC_FASTEST,
                    channels, &err);

    if (err) {
        std::cerr << "Resampler::Resampler: failed to create libsamplerate resampler: " 
                  << src_strerror(err) << std::endl;
#ifndef NO_EXCEPTIONS
        throw Resampler::ImplementationError;
#endif
    }

    if (maxBufferSize > 0 && m_channels > 1) {
        m_iinsize = maxBufferSize * m_channels;
        m_ioutsize = maxBufferSize * m_channels * 2;
        m_iin = allocate<float>(m_iinsize);
        m_iout = allocate<float>(m_ioutsize);
    }

    reset();
}

D_SRC::~D_SRC()
{
    src_delete(m_src);
    deallocate(m_iin);
    deallocate(m_iout);
}

int
D_SRC::resample(const float *const R__ *const R__ in,
                float *const R__ *const R__ out,
                int incount,
                float ratio,
                bool final)
{
    SRC_DATA data;

    int outcount = lrintf(ceilf(incount * ratio));

    if (m_channels == 1) {
        data.data_in = const_cast<float *>(*in); //!!!???
        data.data_out = *out;
    } else {
        if (incount * m_channels > m_iinsize) {
            m_iin = reallocate<float>(m_iin, m_iinsize, incount * m_channels);
            m_iinsize = incount * m_channels;
        }
        if (outcount * m_channels > m_ioutsize) {
            m_iout = reallocate<float>(m_iout, m_ioutsize, outcount * m_channels);
            m_ioutsize = outcount * m_channels;
        }
        v_interleave(m_iin, in, m_channels, incount);
        data.data_in = m_iin;
        data.data_out = m_iout;
    }

    data.input_frames = incount;
    data.output_frames = outcount;
    data.src_ratio = ratio;
    data.end_of_input = (final ? 1 : 0);

    int err = src_process(m_src, &data);

    if (err) {
        std::cerr << "Resampler::process: libsamplerate error: "
                  << src_strerror(err) << std::endl;
#ifndef NO_EXCEPTIONS
        throw Resampler::ImplementationError;
#endif
    }

    if (m_channels > 1) {
        v_deinterleave(out, m_iout, m_channels, data.output_frames_gen);
    }

    m_lastRatio = ratio;

    return data.output_frames_gen;
}

int
D_SRC::resampleInterleaved(const float *const R__ in,
                           float *const R__ out,
                           int incount,
                           float ratio,
                           bool final)
{
    SRC_DATA data;

    int outcount = lrintf(ceilf(incount * ratio));

    data.data_in = const_cast<float *>(in);
    data.data_out = out;

    data.input_frames = incount;
    data.output_frames = outcount;
    data.src_ratio = ratio;
    data.end_of_input = (final ? 1 : 0);

    int err = src_process(m_src, &data);

    if (err) {
        std::cerr << "Resampler::process: libsamplerate error: "
                  << src_strerror(err) << std::endl;
#ifndef NO_EXCEPTIONS
        throw Resampler::ImplementationError;
#endif
    }

    m_lastRatio = ratio;

    return data.output_frames_gen;
}

void
D_SRC::reset()
{
    src_reset(m_src);
}

#endif /* HAVE_LIBSAMPLERATE */

#ifdef HAVE_LIBRESAMPLE

class D_Resample : public ResamplerImpl
{
public:
    D_Resample(Resampler::Quality quality, int channels, int maxBufferSize,
          int m_debugLevel);
    ~D_Resample();

    int resample(const float *const R__ *const R__ in,
                 float *const R__ *const R__ out,
                 int incount,
                 float ratio,
                 bool final);

    int resampleInterleaved(const float *const R__ in,
                            float *const R__ out,
                            int incount,
                            float ratio,
                            bool final);

    int getChannelCount() const { return m_channels; }

    void reset();

protected:
    void *m_src;
    float *m_iin;
    float *m_iout;
    float m_lastRatio;
    int m_channels;
    int m_iinsize;
    int m_ioutsize;
    int m_debugLevel;
};

D_Resample::D_Resample(Resampler::Quality quality, int channels, int maxBufferSize,
             int debugLevel) :
    m_src(0),
    m_iin(0),
    m_iout(0),
    m_lastRatio(1.f),
    m_channels(channels),
    m_iinsize(0),
    m_ioutsize(0),
    m_debugLevel(debugLevel)
{
    if (m_debugLevel > 0) {
        std::cerr << "Resampler::Resampler: using libresample implementation"
                  << std::endl;
    }

    float min_factor = 0.125f;
    float max_factor = 8.0f;

    m_src = resample_open(quality == Resampler::Best ? 1 : 0, min_factor, max_factor);

    if (!m_src) {
        std::cerr << "Resampler::Resampler: failed to create libresample resampler: " 
                  << std::endl;
        throw Resampler::ImplementationError; //!!! of course, need to catch this!
    }

    if (maxBufferSize > 0 && m_channels > 1) {
        m_iinsize = maxBufferSize * m_channels;
        m_ioutsize = maxBufferSize * m_channels * 2;
        m_iin = allocate<float>(m_iinsize);
        m_iout = allocate<float>(m_ioutsize);
    }

    reset();
}

D_Resample::~D_Resample()
{
    resample_close(m_src);
    if (m_iinsize > 0) {
        deallocate(m_iin);
    }
    if (m_ioutsize > 0) {
        deallocate(m_iout);
    }
}

int
D_Resample::resample(const float *const R__ *const R__ in,
                     float *const R__ *const R__ out,
                     int incount,
                     float ratio,
                     bool final)
{
    float *data_in;
    float *data_out;
    int input_frames, output_frames, end_of_input, source_used;
    float src_ratio;

    int outcount = lrintf(ceilf(incount * ratio));

    if (m_channels == 1) {
        data_in = const_cast<float *>(*in); //!!!???
        data_out = *out;
    } else {
        if (incount * m_channels > m_iinsize) {
            m_iin = reallocate<float>(m_iin, m_iinsize, incount * m_channels);
            m_iinsize = incount * m_channels;
        }
        if (outcount * m_channels > m_ioutsize) {
            m_iout = reallocate<float>(m_iout, m_ioutsize, outcount * m_channels);
            m_ioutsize = outcount * m_channels;
        }
        v_interleave(m_iin, in, m_channels, incount);
        data_in = m_iin;
        data_out = m_iout;
    }

    input_frames = incount;
    output_frames = outcount;
    src_ratio = ratio;
    end_of_input = (final ? 1 : 0);

    int output_frames_gen = resample_process(m_src,
                                             src_ratio,
                                             data_in,
                                             input_frames,
                                             end_of_input,
                                             &source_used,
                                             data_out,
                                             output_frames);

    if (output_frames_gen < 0) {
        std::cerr << "Resampler::process: libresample error: "
                  << std::endl;
        throw Resampler::ImplementationError; //!!! of course, need to catch this!
    }

    if (m_channels > 1) {
        v_deinterleave(out, m_iout, m_channels, output_frames_gen);
    }

    m_lastRatio = ratio;

    return output_frames_gen;
}

int
D_Resample::resampleInterleaved(const float *const R__ in,
                                float *const R__ out,
                                int incount,
                                float ratio,
                                bool final)
{
    int input_frames, output_frames, end_of_input, source_used;
    float src_ratio;

    int outcount = lrintf(ceilf(incount * ratio));

    input_frames = incount;
    output_frames = outcount;
    src_ratio = ratio;
    end_of_input = (final ? 1 : 0);

    int output_frames_gen = resample_process(m_src,
                                             src_ratio,
                                             const_cast<float *>(in),
                                             input_frames,
                                             end_of_input,
                                             &source_used,
                                             out,
                                             output_frames);

    if (output_frames_gen < 0) {
        std::cerr << "Resampler::process: libresample error: "
                  << std::endl;
        throw Resampler::ImplementationError; //!!! of course, need to catch this!
    }

    m_lastRatio = ratio;

    return output_frames_gen;
}

void
D_Resample::reset()
{
}

#endif /* HAVE_LIBRESAMPLE */

#ifdef USE_SPEEX
    
class D_Speex : public ResamplerImpl
{
public:
    D_Speex(Resampler::Quality quality, int channels, int maxBufferSize,
            int debugLevel);
    ~D_Speex();

    int resample(const float *const R__ *const R__ in,
                 float *const R__ *const R__ out,
                 int incount,
                 float ratio,
                 bool final);

    int resampleInterleaved(const float *const R__ in,
                            float *const R__ out,
                            int incount,
                            float ratio,
                            bool final = false);

    int getChannelCount() const { return m_channels; }

    void reset();

protected:
    SpeexResamplerState *m_resampler;
    float *m_iin;
    float *m_iout;
    int m_channels;
    int m_iinsize;
    int m_ioutsize;
    float m_lastratio;
    bool m_initial;
    int m_debugLevel;

    void setRatio(float);
};

D_Speex::D_Speex(Resampler::Quality quality, int channels, int maxBufferSize,
                 int debugLevel) :
    m_resampler(0),
    m_iin(0),
    m_iout(0),
    m_channels(channels),
    m_iinsize(0),
    m_ioutsize(0),
    m_lastratio(1),
    m_initial(true),
    m_debugLevel(debugLevel)
{
    int q = (quality == Resampler::Best ? 10 :
             quality == Resampler::Fastest ? 0 : 4);

    if (m_debugLevel > 0) {
        std::cerr << "Resampler::Resampler: using Speex implementation with q = "
                  << q 
                  << std::endl;
    }

    int err = 0;
    m_resampler = speex_resampler_init_frac(m_channels,
                                            1, 1,
                                            48000, 48000, // irrelevant
                                            q,
                                            &err);
    

    if (err) {
        std::cerr << "Resampler::Resampler: failed to create Speex resampler" 
                  << std::endl;
#ifndef NO_EXCEPTIONS
        throw Resampler::ImplementationError;
#endif
    }

    if (maxBufferSize > 0 && m_channels > 1) {
        m_iinsize = maxBufferSize * m_channels;
        m_ioutsize = maxBufferSize * m_channels * 2;
        m_iin = allocate<float>(m_iinsize);
        m_iout = allocate<float>(m_ioutsize);
    }
}

D_Speex::~D_Speex()
{
    speex_resampler_destroy(m_resampler);
    deallocate<float>(m_iin);
    deallocate<float>(m_iout);
}

void
D_Speex::setRatio(float ratio)
{
    // Speex wants a ratio of two unsigned integers, not a single
    // float.  Let's do that.

    unsigned int big = 272408136U; 
    unsigned int denom = 1, num = 1;

    if (ratio < 1.f) {
        denom = big;
        double dnum = double(big) * double(ratio);
        num = (unsigned int)dnum;
    } else if (ratio > 1.f) {
        num = big;
        double ddenom = double(big) / double(ratio);
        denom = (unsigned int)ddenom;
    }
    
    if (m_debugLevel > 1) {
        std::cerr << "D_Speex: Desired ratio " << ratio << ", requesting ratio "
                  << num << "/" << denom << " = " << float(double(num)/double(denom))
                  << std::endl;
    }
    
    int err = speex_resampler_set_rate_frac
        (m_resampler, denom, num, 48000, 48000);
    //!!! check err
    
    speex_resampler_get_ratio(m_resampler, &denom, &num);
    
    if (m_debugLevel > 1) {
        std::cerr << "D_Speex: Desired ratio " << ratio << ", got ratio "
                  << num << "/" << denom << " = " << float(double(num)/double(denom))
                  << std::endl;
    }
    
    m_lastratio = ratio;

    if (m_initial) {
        speex_resampler_skip_zeros(m_resampler);
        m_initial = false;
    }
}

int
D_Speex::resample(const float *const R__ *const R__ in,
                  float *const R__ *const R__ out,
                  int incount,
                  float ratio,
                  bool final)
{
    if (ratio != m_lastratio) {
        setRatio(ratio);
    }

    unsigned int uincount = incount;
    unsigned int outcount = lrintf(ceilf(incount * ratio)); //!!! inexact now

    float *data_in, *data_out;

    if (m_channels == 1) {
        data_in = const_cast<float *>(*in);
        data_out = *out;
    } else {
        if (incount * m_channels > m_iinsize) {
            m_iin = reallocate<float>(m_iin, m_iinsize, incount * m_channels);
            m_iinsize = incount * m_channels;
        }
        if (outcount * m_channels > m_ioutsize) {
            m_iout = reallocate<float>(m_iout, m_ioutsize, outcount * m_channels);
            m_ioutsize = outcount * m_channels;
        }
        v_interleave(m_iin, in, m_channels, incount);
        data_in = m_iin;
        data_out = m_iout;
    }

    int err = speex_resampler_process_interleaved_float(m_resampler,
                                                        data_in,
                                                        &uincount,
                                                        data_out,
                                                        &outcount);

//    if (incount != int(uincount)) {
//        std::cerr << "Resampler: NOTE: Consumed " << uincount
//                  << " of " << incount << " frames" << std::endl;
//    }

//    if (outcount != lrintf(ceilf(incount * ratio))) {
//        std::cerr << "Resampler: NOTE: Obtained " << outcount
//                  << " of " << lrintf(ceilf(incount * ratio)) << " frames"
//                  << std::endl;
//    }
        
    //!!! check err, respond appropriately


    if (m_channels > 1) {
        v_deinterleave(out, m_iout, m_channels, outcount);
    }

    return outcount;
}

int
D_Speex::resampleInterleaved(const float *const R__ in,
                             float *const R__ out,
                             int incount,
                             float ratio,
                             bool final)
{
    if (ratio != m_lastratio) {
        setRatio(ratio);
    }

    unsigned int uincount = incount;
    unsigned int outcount = lrintf(ceilf(incount * ratio)); //!!! inexact now

    float *data_in = const_cast<float *>(in);
    float *data_out = out;

    int err = speex_resampler_process_interleaved_float(m_resampler,
                                                        data_in,
                                                        &uincount,
                                                        data_out,
                                                        &outcount);

//    std::cerr << "D_SPEEX: incount " << incount << " ratio " << ratio << " req " << lrintf(ceilf(incount * ratio)) << " final " << final << " output_frames_gen " << outcount << std::endl;

    return outcount;
}

void
D_Speex::reset()
{
    m_lastratio = -1.0; // force reset of ratio
    m_initial = true;
    speex_resampler_reset_mem(m_resampler);
}

#endif

} /* end namespace Resamplers */

Resampler::Resampler(Resampler::Quality quality, int channels,
                     int maxBufferSize, int debugLevel)
{
    m_method = -1;
    
    switch (quality) {

    case Resampler::Best:
#ifdef HAVE_IPP
        m_method = 0;
#endif
#ifdef USE_SPEEX
        m_method = 2;
#endif
#ifdef HAVE_LIBRESAMPLE
        m_method = 3;
#endif
#ifdef HAVE_LIBSAMPLERATE
        m_method = 1;
#endif
        break;

    case Resampler::FastestTolerable:
#ifdef HAVE_IPP
        m_method = 0;
#endif
#ifdef HAVE_LIBRESAMPLE
        m_method = 3;
#endif
#ifdef HAVE_LIBSAMPLERATE
        m_method = 1;
#endif
#ifdef USE_SPEEX
        m_method = 2;
#endif
        break;

    case Resampler::Fastest:
#ifdef HAVE_IPP
        m_method = 0;
#endif
#ifdef HAVE_LIBRESAMPLE
        m_method = 3;
#endif
#ifdef USE_SPEEX
        m_method = 2;
#endif
#ifdef HAVE_LIBSAMPLERATE
        m_method = 1;
#endif
        break;
    }

    if (m_method == -1) {
        std::cerr << "Resampler::Resampler(" << quality << ", " << channels
                  << ", " << maxBufferSize << "): No implementation available!"
                  << std::endl;
        abort();
    }

    switch (m_method) {
    case 0:
#ifdef HAVE_IPP
        d = new Resamplers::D_IPP(quality, channels, maxBufferSize, debugLevel);
#else
        std::cerr << "Resampler::Resampler(" << quality << ", " << channels
                  << ", " << maxBufferSize << "): No implementation available!"
                  << std::endl;
        abort();
#endif
        break;

    case 1:
#ifdef HAVE_LIBSAMPLERATE
        d = new Resamplers::D_SRC(quality, channels, maxBufferSize, debugLevel);
#else
        std::cerr << "Resampler::Resampler(" << quality << ", " << channels
                  << ", " << maxBufferSize << "): No implementation available!"
                  << std::endl;
        abort();
#endif
        break;

    case 2:
#ifdef USE_SPEEX
        d = new Resamplers::D_Speex(quality, channels, maxBufferSize, debugLevel);
#else
        std::cerr << "Resampler::Resampler(" << quality << ", " << channels
                  << ", " << maxBufferSize << "): No implementation available!"
                  << std::endl;
        abort();
#endif
        break;

    case 3:
#ifdef HAVE_LIBRESAMPLE
        d = new Resamplers::D_Resample(quality, channels, maxBufferSize, debugLevel);
#else
        std::cerr << "Resampler::Resampler(" << quality << ", " << channels
                  << ", " << maxBufferSize << "): No implementation available!"
                  << std::endl;
        abort();
#endif
        break;
    }

    if (!d) {
        std::cerr << "Resampler::Resampler(" << quality << ", " << channels
                  << ", " << maxBufferSize
                  << "): Internal error: No implementation selected"
                  << std::endl;
        abort();
    }
}

Resampler::~Resampler()
{
    delete d;
}

int 
Resampler::resample(const float *const R__ *const R__ in,
                    float *const R__ *const R__ out,
                    int incount, float ratio, bool final)
{
    Profiler profiler("Resampler::resample");
    return d->resample(in, out, incount, ratio, final);
}

int 
Resampler::resampleInterleaved(const float *const R__ in,
                               float *const R__ out,
                               int incount, float ratio, bool final)
{
    Profiler profiler("Resampler::resample");
    return d->resampleInterleaved(in, out, incount, ratio, final);
}

int
Resampler::getChannelCount() const
{
    return d->getChannelCount();
}

void
Resampler::reset()
{
    d->reset();
}

}
