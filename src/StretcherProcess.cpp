/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2008 Chris Cannam.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "StretcherImpl.h"
#include "PercussiveAudioCurve.h"
#include "HighFrequencyAudioCurve.h"
#include "ConstantAudioCurve.h"
#include "StretchCalculator.h"
#include "StretcherChannelData.h"
#include "Resampler.h"
#include "Profiler.h"

#include <cassert>
#include <cmath>
#include <set>
#include <map>
#include <deque>


using std::cerr;
using std::endl;

namespace RubberBand {

RubberBandStretcher::Impl::ProcessThread::ProcessThread(Impl *s, size_t c) :
    m_s(s),
    m_channel(c),
    m_dataAvailable(std::string("data ") + char('A' + c)),
    m_abandoning(false)
{ }

void
RubberBandStretcher::Impl::ProcessThread::run()
{
    if (m_s->m_debugLevel > 1) {
        cerr << "thread " << m_channel << " getting going" << endl;
    }

    ChannelData &cd = *m_s->m_channelData[m_channel];

    while (cd.inputSize == -1 ||
           cd.inbuf->getReadSpace() > 0) {

//        if (cd.inputSize != -1) {
//            cerr << "inputSize == " << cd.inputSize
//                 << ", readSpace == " << cd.inbuf->getReadSpace() << endl;
//        }
        
        bool any = false, last = false;
        m_s->processChunks(m_channel, any, last);

        if (last) break;

        if (any) m_s->m_spaceAvailable.signal();

        m_dataAvailable.lock();
        if (!m_s->testInbufReadSpace(m_channel) && !m_abandoning) {
            m_dataAvailable.wait();
        } else {
            m_dataAvailable.unlock();
        }

        if (m_abandoning) {
            if (m_s->m_debugLevel > 1) {
                cerr << "thread " << m_channel << " abandoning" << endl;
            }
            return;
        }
    }

    bool any = false, last = false;
    m_s->processChunks(m_channel, any, last);
    m_s->m_spaceAvailable.signal();
    
    if (m_s->m_debugLevel > 1) {
        cerr << "thread " << m_channel << " done" << endl;
    }
}

void
RubberBandStretcher::Impl::ProcessThread::signalDataAvailable()
{
    m_dataAvailable.signal();
}

void
RubberBandStretcher::Impl::ProcessThread::abandon()
{
    m_abandoning = true;
}

bool
RubberBandStretcher::Impl::resampleBeforeStretching() const
{
    // We can't resample before stretching in offline mode, because
    // the stretch calculation is based on doing it the other way
    // around.  It would take more work (and testing) to enable this.
    if (!m_realtime) return false;

    if (m_options & OptionPitchHighQuality) {
        return (m_pitchScale < 1.0); // better sound
    } else {
        return (m_pitchScale > 1.0); // better performance
    }
}
    
size_t
RubberBandStretcher::Impl::consumeChannel(size_t c, const float *input,
                                          size_t samples, bool final)
{
    Profiler profiler("RubberBandStretcher::Impl::consumeChannel");

    ChannelData &cd = *m_channelData[c];
    RingBuffer<float> &inbuf = *cd.inbuf;

    size_t toWrite = samples;
    size_t writable = inbuf.getWriteSpace();

    bool resampling = resampleBeforeStretching();

    if (resampling) {

        toWrite = int(ceil(samples / m_pitchScale));
        if (writable < toWrite) {
            samples = int(floor(writable * m_pitchScale));
            if (samples == 0) return 0;
        }

        size_t reqSize = int(ceil(samples / m_pitchScale));
        if (reqSize > cd.resamplebufSize) {
            cerr << "WARNING: RubberBandStretcher::Impl::consumeChannel: resizing resampler buffer from "
                 << cd.resamplebufSize << " to " << reqSize << endl;
            cd.setResampleBufSize(reqSize);
        }


        toWrite = cd.resampler->resample(&input,
                                         &cd.resamplebuf,
                                         samples,
                                         1.0 / m_pitchScale,
                                         final);

    }

    if (writable < toWrite) {
        if (resampling) {
            return 0;
        }
        toWrite = writable;
    }

    if (resampling) {
        inbuf.write(cd.resamplebuf, toWrite);
        cd.inCount += samples;
        return samples;
    } else {
        inbuf.write(input, toWrite);
        cd.inCount += toWrite;
        return toWrite;
    }
}

void
RubberBandStretcher::Impl::processChunks(size_t c, bool &any, bool &last)
{
    Profiler profiler("RubberBandStretcher::Impl::processChunks");

    // Process as many chunks as there are available on the input
    // buffer for channel c.  This requires that the increments have
    // already been calculated.

    ChannelData &cd = *m_channelData[c];

    last = false;
    any = false;

    while (!last) {

        if (!testInbufReadSpace(c)) {
//            cerr << "not enough input" << endl;
            break;
        }

        any = true;

        if (!cd.draining) {
            size_t got = cd.inbuf->peek(cd.fltbuf, m_windowSize);
            assert(got == m_windowSize || cd.inputSize >= 0);
            cd.inbuf->skip(m_increment);
            analyseChunk(c);
        }

        bool phaseReset = false;
        size_t phaseIncrement, shiftIncrement;
        getIncrements(c, phaseIncrement, shiftIncrement, phaseReset);

        last = processChunkForChannel(c, phaseIncrement, shiftIncrement, phaseReset);
        cd.chunkCount++;
        if (m_debugLevel > 2) {
            cerr << "channel " << c << ": last = " << last << ", chunkCount = " << cd.chunkCount << endl;
        }
    }
}

bool
RubberBandStretcher::Impl::processOneChunk()
{
    Profiler profiler("RubberBandStretcher::Impl::processOneChunk");

    // Process a single chunk for all channels, provided there is
    // enough data on each channel for at least one chunk.  This is
    // able to calculate increments as it goes along.

    for (size_t c = 0; c < m_channels; ++c) {
        if (!testInbufReadSpace(c)) return false;
        ChannelData &cd = *m_channelData[c];
        if (!cd.draining) {
            size_t got = cd.inbuf->peek(cd.fltbuf, m_windowSize);
            assert(got == m_windowSize || cd.inputSize >= 0);
            cd.inbuf->skip(m_increment);
            analyseChunk(c);
        }
    }
    
    bool phaseReset = false;
    size_t phaseIncrement, shiftIncrement;
    if (!getIncrements(0, phaseIncrement, shiftIncrement, phaseReset)) {
        calculateIncrements(phaseIncrement, shiftIncrement, phaseReset);
    }

    bool last = false;
    for (size_t c = 0; c < m_channels; ++c) {
        last = processChunkForChannel(c, phaseIncrement, shiftIncrement, phaseReset);
        m_channelData[c]->chunkCount++;
    }

    return last;
}

bool
RubberBandStretcher::Impl::testInbufReadSpace(size_t c)
{
    Profiler profiler("RubberBandStretcher::Impl::testInbufReadSpace");

    ChannelData &cd = *m_channelData[c];
    RingBuffer<float> &inbuf = *cd.inbuf;

    size_t rs = inbuf.getReadSpace();

    if (rs < m_windowSize && !cd.draining) {
            
        if (cd.inputSize == -1) {

            // Not all the input data has been written to the inbuf
            // (that's why the input size is not yet set).  We can't
            // process, because we don't have a full chunk of data, so
            // our process chunk would contain some empty padding in
            // its input -- and that would give incorrect output, as
            // we know there is more input to come.

            if (!m_threaded) {
//                cerr << "WARNING: RubberBandStretcher: read space < chunk size ("
//                          << inbuf.getReadSpace() << " < " << m_windowSize
//                          << ") when not all input written, on processChunks for channel " << c << endl;
            }
            return false;
        }
        
        if (rs == 0) {

            if (m_debugLevel > 1) {
                cerr << "read space = 0, giving up" << endl;
            }
            return false;

        } else if (rs < m_windowSize/2) {

            if (m_debugLevel > 1) {
                cerr << "read space = " << rs << ", setting draining true" << endl;
            }
            
            cd.draining = true;
        }
    }

    return true;
}

bool 
RubberBandStretcher::Impl::processChunkForChannel(size_t c,
                                                  size_t phaseIncrement,
                                                  size_t shiftIncrement,
                                                  bool phaseReset)
{
    Profiler profiler("RubberBandStretcher::Impl::processChunkForChannel");

    // Process a single chunk on a single channel.  This assumes
    // enough input data is available; caller must have tested this
    // using e.g. testInbufReadSpace first.  Return true if this is
    // the last chunk on the channel.

    if (phaseReset && (m_debugLevel > 1)) {
        cerr << "processChunkForChannel: phase reset found, incrs "
                  << phaseIncrement << ":" << shiftIncrement << endl;
    }

    ChannelData &cd = *m_channelData[c];

    if (!cd.draining) {
        
        // This is the normal processing case -- draining is only
        // set when all the input has been used and we only need
        // to write from the existing accumulator into the output.
        
        // We know we have enough samples available in m_inbuf --
        // this is usually m_windowSize, but we know that if fewer
        // are available, it's OK to use zeroes for the rest
        // (which the ring buffer will provide) because we've
        // reached the true end of the data.
        
        // We need to peek m_windowSize samples for processing, and
        // then skip m_increment to advance the read pointer.
        
        modifyChunk(c, phaseIncrement, phaseReset);
        synthesiseChunk(c); // reads from cd.mag, cd.phase

        if (m_debugLevel > 2) {
            if (phaseReset) {
                for (int i = 0; i < 10; ++i) {
                    cd.accumulator[i] = 1.2f - (i % 3) * 1.2f;
                }
            }
        }
    }

    bool last = false;

    if (cd.draining) {
        if (m_debugLevel > 1) {
            cerr << "draining: accumulator fill = " << cd.accumulatorFill << " (shiftIncrement = " << shiftIncrement << ")" <<  endl;
        }
        if (shiftIncrement == 0) {
            cerr << "WARNING: draining: shiftIncrement == 0, can't handle that in this context: setting to " << m_increment << endl;
            shiftIncrement = m_increment;
        }
        if (cd.accumulatorFill <= shiftIncrement) {
            if (m_debugLevel > 1) {
                cerr << "reducing shift increment from " << shiftIncrement
                          << " to " << cd.accumulatorFill
                          << " and marking as last" << endl;
            }
            shiftIncrement = cd.accumulatorFill;
            last = true;
        }
    }
        
    if (m_threaded) {

        int required = shiftIncrement;

        if (m_pitchScale != 1.0) {
            required = int(required / m_pitchScale) + 1;
        }
        
        if (cd.outbuf->getWriteSpace() < required) {
            if (m_debugLevel > 0) {
                cerr << "Buffer overrun on output for channel " << c << endl;
            }

            //!!! The only correct thing we can do here is resize the
            // buffer.  We can't wait for the client thread to read
            // some data out from the buffer so as to make more space,
            // because the client thread is probably stuck in a
            // process() call waiting for us to stow away enough input
            // increments to allow the process() call to complete.

        }
    }
    
    writeChunk(c, shiftIncrement, last);
    return last;
}

void
RubberBandStretcher::Impl::calculateIncrements(size_t &phaseIncrementRtn,
                                               size_t &shiftIncrementRtn,
                                               bool &phaseReset)
{
    Profiler profiler("RubberBandStretcher::Impl::calculateIncrements");

//    cerr << "calculateIncrements" << endl;
    
    // Calculate the next upcoming phase and shift increment, on the
    // basis that both channels are in sync.  This is in contrast to
    // getIncrements, which requires that all the increments have been
    // calculated in advance but can then return increments
    // corresponding to different chunks in different channels.

    // Requires frequency domain representations of channel data in
    // the mag and phase buffers in the channel.

    // This function is only used in real-time mode.

    phaseIncrementRtn = m_increment;
    shiftIncrementRtn = m_increment;
    phaseReset = false;

    if (m_channels == 0) return;

    ChannelData &cd = *m_channelData[0];

    size_t bc = cd.chunkCount;
    for (size_t c = 1; c < m_channels; ++c) {
        if (m_channelData[c]->chunkCount != bc) {
            cerr << "ERROR: RubberBandStretcher::Impl::calculateIncrements: Channels are not in sync" << endl;
            return;
        }
    }

    const int hs = m_windowSize/2 + 1;

    // Normally we would mix down the time-domain signal and apply a
    // single FFT, or else mix down the Cartesian form of the
    // frequency-domain signal.  Both of those would be inefficient
    // from this position.  Fortunately, the onset detectors should
    // work reasonably well (maybe even better?) if we just sum the
    // magnitudes of the frequency-domain channel signals and forget
    // about phase entirely.  Normally we don't expect the channel
    // phases to cancel each other, and broadband effects will still
    // be apparent.

    float df = 0.f;

    if (m_channels == 1) {

        df = m_phaseResetAudioCurve->process(cd.mag, m_increment);

    } else {

        double *tmp = (double *)alloca(hs * sizeof(double));

        for (int i = 0; i < hs; ++i) {
            tmp[i] = 0.0;
        }
        for (size_t c = 0; c < m_channels; ++c) {
            for (int i = 0; i < hs; ++i) {
                tmp[i] += m_channelData[c]->mag[i];
            }
        }
    
        df = m_phaseResetAudioCurve->process(tmp, m_increment);
	}

    int incr = m_stretchCalculator->calculateSingle
            (getEffectiveRatio(), df, m_increment);

    m_lastProcessPhaseResetDf.write(&df, 1);
    m_lastProcessOutputIncrements.write(&incr, 1);

    if (incr < 0) {
        phaseReset = true;
        incr = -incr;
    }
    
    // The returned increment is the phase increment.  The shift
    // increment for one chunk is the same as the phase increment for
    // the following chunk (see comment below).  This means we don't
    // actually know the shift increment until we see the following
    // phase increment... which is a bit of a problem.

    // This implies we should use this increment for the shift
    // increment, and make the following phase increment the same as
    // it.  This means in RT mode we'll be one chunk later with our
    // phase reset than we would be in non-RT mode.  The sensitivity
    // of the broadband onset detector may mean that this isn't a
    // problem -- test it and see.

    shiftIncrementRtn = incr;

    if (cd.prevIncrement == 0) {
        phaseIncrementRtn = shiftIncrementRtn;
    } else {
        phaseIncrementRtn = cd.prevIncrement;
    }

    cd.prevIncrement = shiftIncrementRtn;
}

bool
RubberBandStretcher::Impl::getIncrements(size_t channel,
                                         size_t &phaseIncrementRtn,
                                         size_t &shiftIncrementRtn,
                                         bool &phaseReset)
{
    Profiler profiler("RubberBandStretcher::Impl::getIncrements");

    if (channel >= m_channels) {
        phaseIncrementRtn = m_increment;
        shiftIncrementRtn = m_increment;
        phaseReset = false;
        return false;
    }

    // There are two relevant output increments here.  The first is
    // the phase increment which we use when recalculating the phases
    // for the current chunk; the second is the shift increment used
    // to determine how far to shift the processing buffer after
    // writing the chunk.  The shift increment for one chunk is the
    // same as the phase increment for the following chunk.
    
    // When an onset occurs for which we need to reset phases, the
    // increment given will be negative.
    
    // When we reset phases, the previous shift increment (and so
    // current phase increments) must have been m_increment to ensure
    // consistency.
    
    // m_outputIncrements stores phase increments.

    ChannelData &cd = *m_channelData[channel];
    bool gotData = true;

    if (cd.chunkCount >= m_outputIncrements.size()) {
//        cerr << "WARNING: RubberBandStretcher::Impl::getIncrements:"
//             << " chunk count " << cd.chunkCount << " >= "
//             << m_outputIncrements.size() << endl;
        if (m_outputIncrements.size() == 0) {
            phaseIncrementRtn = m_increment;
            shiftIncrementRtn = m_increment;
            phaseReset = false;
            return false;
        } else {
            cd.chunkCount = m_outputIncrements.size()-1;
            gotData = false;
        }
    }
    
    int phaseIncrement = m_outputIncrements[cd.chunkCount];
    
    int shiftIncrement = phaseIncrement;
    if (cd.chunkCount + 1 < m_outputIncrements.size()) {
        shiftIncrement = m_outputIncrements[cd.chunkCount + 1];
    }
    
    if (phaseIncrement < 0) {
        phaseIncrement = -phaseIncrement;
        phaseReset = true;
    }
    
    if (shiftIncrement < 0) {
        shiftIncrement = -shiftIncrement;
    }
    
    if (shiftIncrement >= int(m_windowSize)) {
        cerr << "*** ERROR: RubberBandStretcher::Impl::processChunks: shiftIncrement " << shiftIncrement << " >= windowSize " << m_windowSize << " at " << cd.chunkCount << " (of " << m_outputIncrements.size() << ")" << endl;
        shiftIncrement = m_windowSize;
    }

    phaseIncrementRtn = phaseIncrement;
    shiftIncrementRtn = shiftIncrement;
    if (cd.chunkCount == 0) phaseReset = true; // don't mess with the first chunk
    return gotData;
}

void
RubberBandStretcher::Impl::analyseChunk(size_t channel)
{
    Profiler profiler("RubberBandStretcher::Impl::analyseChunk");

    int i;

    ChannelData &cd = *m_channelData[channel];

    double *const R__ dblbuf = cd.dblbuf;
    float *const R__ fltbuf = cd.fltbuf;

    int sz = m_windowSize;
    int hs = m_windowSize/2;

    // cd.fltbuf is known to contain m_windowSize samples

    m_window->cut(fltbuf);

    if (cd.oversample > 1) {

        int bufsiz = sz * cd.oversample;
        int offset = (bufsiz - sz) / 2;

        // eek

        for (i = 0; i < offset; ++i) {
            dblbuf[i] = 0.0;
        }
        for (i = 0; i < offset; ++i) {
            dblbuf[bufsiz - i - 1] = 0.0;
        }
        for (i = 0; i < sz; ++i) {
            dblbuf[offset + i] = fltbuf[i];
        }
        for (i = 0; i < bufsiz / 2; ++i) {
            double tmp = dblbuf[i];
            dblbuf[i] = dblbuf[i + bufsiz/2];
            dblbuf[i + bufsiz/2] = tmp;
        }
    } else {
        for (i = 0; i < hs; ++i) {
            dblbuf[i] = fltbuf[i + hs];
            dblbuf[i + hs] = fltbuf[i];
        }
    }

    cd.fft->forwardPolar(dblbuf, cd.mag, cd.phase);
}

static inline double mod(double x, double y) { return x - (y * floor(x / y)); }
static inline double princarg(double a) { return mod(a + M_PI, -2.0 * M_PI) + M_PI; }

void
RubberBandStretcher::Impl::modifyChunk(size_t channel, size_t outputIncrement,
                                       bool phaseReset)
{
    Profiler profiler("RubberBandStretcher::Impl::modifyChunk");

    ChannelData &cd = *m_channelData[channel];

    if (phaseReset && m_debugLevel > 1) {
        cerr << "phase reset: leaving phases unmodified" << endl;
    }

    int pfp = 0;
    double rate = m_sampleRate;

    int sz = m_windowSize;

    int count = (sz * cd.oversample) / 2;

    bool unchanged = cd.unchanged && (outputIncrement == m_increment);
    bool fullReset = phaseReset;

    if (!(m_options & OptionPhaseIndependent)) {

        cd.freqPeak[0] = 0;

        float freq0 = m_freq0;
        float freq1 = m_freq1;
        float freq2 = m_freq2;

        // As the stretch ratio increases, so the frequency thresholds
        // for phase lamination should increase.  Beyond a ratio of
        // about 1.5, the threshold should be about 1200Hz; beyond a
        // ratio of 2, we probably want no lamination to happen at all
        // by default.  This calculation aims for more or less that.
        // We only do this if the phase option is OptionPhaseAdaptive
        // (the default), i.e. not Independent or PeakLocked.

        if (!(m_options & OptionPhasePeakLocked)) {
            float r = getEffectiveRatio();
            if (r > 1) {
                float rf0 = 600 + (600 * ((r-1)*(r-1)*(r-1)*2));
                float f1ratio = freq1 / freq0;
                float f2ratio = freq2 / freq0;
                freq0 = std::max(freq0, rf0);
                freq1 = freq0 * f1ratio;
                freq2 = freq0 * f2ratio;
            }
        }

        int limit0 = lrint((freq0 * sz * cd.oversample) / rate);
        int limit1 = lrint((freq1 * sz * cd.oversample) / rate);
        int limit2 = lrint((freq2 * sz * cd.oversample) / rate);

        int range = 0;

        if (limit1 < limit0) limit1 = limit0;
        if (limit2 < limit1) limit2 = limit1;
    
//        cerr << "limit0 = " << limit0 << " limit1 = " << limit1 << " limit2 = " << limit2 << endl;

        int peakCount = 0;

        for (int i = 0; i <= count; ++i) {

            double mag = cd.mag[i];
            bool isPeak = true;

            for (int j = 1; j <= range; ++j) {

                if (mag < cd.mag[i-j]) {
                    isPeak = false;
                    break;
                }

                if (mag < cd.mag[i+j]) {
                    isPeak = false;
                    break;
                }
            }        

            if (isPeak) {

                // i is a peak bin.

                // The previous peak bin was at pfp; make freqPeak entries
                // from pfp to half-way between pfp and i point at pfp, and
                // those from the half-way mark to i point at i.
            
                int halfway = (pfp + i) / 2;
                if (halfway == pfp) halfway = pfp + 1;

                for (int j = pfp + 1; j < halfway; ++j) {
                    cd.freqPeak[j] = pfp;
                }
                for (int j = halfway; j <= i; ++j) {
                    cd.freqPeak[j] = i;
                }

                pfp = i;

                ++peakCount;
            }

            if (i == limit0) range = 1;
            if (i == limit1) range = 2;
            if (i >= limit2) {
                range = 3;
                if (i + range + 1 > count) range = count - i - 1;
            }
        }

//        cerr << "peakCount = " << peakCount << endl;
        
        cd.freqPeak[count-1] = count-1;
        cd.freqPeak[count] = count;
    }

    Profiler profiler2("RubberBandStretcher::Impl::modifyChunk part 2");

    double peakInPhase = 0.0;
    double peakOutPhase = 0.0;
    int p = -1, pp = -1;


    for (int i = 0; i <= count; ++i) {
        
        if (m_options & OptionPhaseIndependent) {
            p = i;
            pp = i-1;
        } else {
            p = cd.freqPeak[i];
            if (i > 0) pp = cd.freqPeak[i-1];
        }

        bool resetThis = phaseReset;
        
        if (m_options & OptionTransientsMixed) {
            int low = lrint((150 * sz * cd.oversample) / rate);
            int high = lrint((1000 * sz * cd.oversample) / rate);
            if (resetThis) {
                if (i > low && i < high) {
                    resetThis = false;
                    fullReset = false;
                }
            }
        }

        if (!resetThis) {

            if (i == 0 || p != pp) {
	

                double omega = (2 * M_PI * m_increment * p) /
                    (sz * cd.oversample);
                double expectedPhase = cd.prevPhase[p] + omega;
                double phaseError = princarg(cd.phase[p] - expectedPhase);
                double phaseIncrement = (omega + phaseError) / m_increment;
            
                double unwrappedPhase = cd.unwrappedPhase[p] +
                    outputIncrement * phaseIncrement;

                cd.prevPhase[p] = cd.phase[p];
                cd.phase[p] = unwrappedPhase;
                cd.unwrappedPhase[p] = unwrappedPhase;

                peakInPhase = cd.prevPhase[p];
                peakOutPhase = unwrappedPhase;

            }

            if (i != p) {

//                cd.phase[i] = atan2(cd.imag[i], cd.real[i]);//!!!

                double diffToPeak = peakInPhase - cd.phase[i];
                double unwrappedPhase = peakOutPhase - diffToPeak;
                
                cd.prevPhase[i] = cd.phase[i];
                cd.phase[i] = unwrappedPhase;
                cd.unwrappedPhase[i] = unwrappedPhase;
            }

        } else {
            // resetThis == true
            cd.prevPhase[i] = cd.phase[i];
            cd.unwrappedPhase[i] = cd.phase[i];
        }
    }


    if (fullReset) unchanged = true;
    cd.unchanged = unchanged;

    if (unchanged && m_debugLevel > 1) {
        cerr << "frame unchanged on channel " << channel << endl;
    }
}


void
RubberBandStretcher::Impl::formantShiftChunk(size_t channel)
{
    Profiler profiler("RubberBandStretcher::Impl::formantShiftChunk");

    ChannelData &cd = *m_channelData[channel];

    double *const R__ mag = cd.mag;
    double *const R__ envelope = cd.envelope;
    double *const R__ dblbuf = cd.dblbuf;

    const int sz = m_windowSize;
    const int hs = m_windowSize/2;
    const double denom = sz;

    
    cd.fft->inverseCepstral(mag, dblbuf);

    for (int i = 0; i < sz; ++i) {
        dblbuf[i] /= denom;
    }

    //!!! calculate this value -- the divisor should be the highest fundamental frequency we expect to find, plus a bit
    const int cutoff = m_sampleRate / 700;
//    const int cutoff = 1000;
//    const int cutoff = 20;

//    cerr <<"cutoff = "<< cutoff << ", m_sampleRate/cutoff = " << m_sampleRate/cutoff << ", m_sampleRate/700 = " << m_sampleRate/700 << endl;

    dblbuf[0] /= 2;
    dblbuf[cutoff-1] /= 2;

    for (int i = cutoff; i < sz; ++i) {
        dblbuf[i] = 0.0;
    }

    cd.fft->forward(dblbuf, envelope, 0);


    for (int i = 0; i <= hs; ++i) {
        envelope[i] = exp(envelope[i]);
    }
    for (int i = 0; i <= hs; ++i) {
        mag[i] /= envelope[i];
    }

    if (m_pitchScale > 1.0) {
        // scaling up, we want a new envelope that is lower by the pitch factor
        for (int target = 0; target <= hs; ++target) {
            int source = lrint(target * m_pitchScale);
            if (source > m_windowSize) {
                envelope[target] = 0.0;
            } else {
                envelope[target] = envelope[source];
            }
        }
    } else {
        // scaling down, we want a new envelope that is higher by the pitch factor
        for (int target = hs; target > 0; ) {
            --target;
            int source = lrint(target * m_pitchScale);
            envelope[target] = envelope[source];
        }
    }

    for (int i = 0; i <= hs; ++i) {
        mag[i] *= envelope[i];
//        mag[i] = envelope[i];
    }

    cd.unchanged = false;
}

void
RubberBandStretcher::Impl::synthesiseChunk(size_t channel)
{
    Profiler profiler("RubberBandStretcher::Impl::synthesiseChunk");


    if ((m_options & OptionFormantPreserved) &&
        (m_pitchScale != 1.0)) {
        formantShiftChunk(channel);
    }

    ChannelData &cd = *m_channelData[channel];

    double *const R__ dblbuf = cd.dblbuf;
    float *const R__ fltbuf = cd.fltbuf;
    float *const R__ accumulator = cd.accumulator;
    float *const R__ windowAccumulator = cd.windowAccumulator;
    
    int sz = m_windowSize;
    int hs = m_windowSize/2;
    int i;


    if (!cd.unchanged) {

        cd.fft->inversePolar(cd.mag, cd.phase, cd.dblbuf);

        if (cd.oversample > 1) {

            int bufsiz = sz * cd.oversample;
            int hbs = hs * cd.oversample;
            int offset = (bufsiz - sz) / 2;

            for (i = 0; i < hbs; ++i) {
                double tmp = dblbuf[i];
                dblbuf[i] = dblbuf[i + hbs];
                dblbuf[i + hbs] = tmp;
            }
            for (i = 0; i < sz; ++i) {
                fltbuf[i] = float(dblbuf[i + offset]);
            }
        } else {
            for (i = 0; i < hs; ++i) {
                fltbuf[i] = float(dblbuf[i + hs]);
            }
            for (i = 0; i < hs; ++i) {
                fltbuf[i + hs] = float(dblbuf[i]);
            }
        }

        float denom = float(sz * cd.oversample);

        // our ffts produced unscaled results
        for (i = 0; i < sz; ++i) {
            fltbuf[i] = fltbuf[i] / float(sz * cd.oversample);
        }
//    } else {
//        cerr << "unchanged on channel " << channel << endl;
    }

    m_window->cut(fltbuf);

    for (i = 0; i < sz; ++i) {
        accumulator[i] += fltbuf[i];
    }

    cd.accumulatorFill = m_windowSize;

    float fixed = m_window->getArea() * 1.5f;

    for (i = 0; i < sz; ++i) {
        float val = m_window->getValue(i);
        windowAccumulator[i] += val * fixed;
    }
}

void
RubberBandStretcher::Impl::writeChunk(size_t channel, size_t shiftIncrement, bool last)
{
    Profiler profiler("RubberBandStretcher::Impl::writeChunk");

    ChannelData &cd = *m_channelData[channel];
    
    float *const R__ accumulator = cd.accumulator;
    float *const R__ windowAccumulator = cd.windowAccumulator;

    const int sz = m_windowSize;
    const int si = shiftIncrement;

    int i;
    
    if (m_debugLevel > 2) {
        cerr << "writeChunk(" << channel << ", " << shiftIncrement << ", " << last << ")" << endl;
    }

    for (i = 0; i < si; ++i) {
        if (windowAccumulator[i] > 0.f) {
            accumulator[i] /= windowAccumulator[i];
        }
    }

    // for exact sample scaling (probably not meaningful if we
    // were running in RT mode)
    size_t theoreticalOut = 0;
    if (cd.inputSize >= 0) {
        theoreticalOut = lrint(cd.inputSize * m_timeRatio);
    }

    bool resampledAlready = resampleBeforeStretching();

    if (!resampledAlready && m_pitchScale != 1.0 && cd.resampler) {

        size_t reqSize = int(ceil(si / m_pitchScale));
        if (reqSize > cd.resamplebufSize) {
            // This shouldn't normally happen -- the buffer is
            // supposed to be initialised with enough space in the
            // first place.  But we retain this check in case the
            // pitch scale has changed since then, or the stretch
            // calculator has gone mad, or something.
            cerr << "WARNING: RubberBandStretcher::Impl::writeChunk: resizing resampler buffer from "
                      << cd.resamplebufSize << " to " << reqSize << endl;
            cd.setResampleBufSize(reqSize);
        }


        size_t outframes = cd.resampler->resample(&cd.accumulator,
                                                  &cd.resamplebuf,
                                                  si,
                                                  1.0 / m_pitchScale,
                                                  last);


        writeOutput(*cd.outbuf, cd.resamplebuf,
                    outframes, cd.outCount, theoreticalOut);

    } else {
        writeOutput(*cd.outbuf, accumulator,
                    si, cd.outCount, theoreticalOut);
    }
    
    for (i = 0; i < sz - si; ++i) {
        accumulator[i] = accumulator[i + si];
    }
    
    for (i = sz - si; i < sz; ++i) {
        accumulator[i] = 0.0f;
    }
    
    for (i = 0; i < sz - si; ++i) {
        windowAccumulator[i] = windowAccumulator[i + si];
    }
    
    for (i = sz - si; i < sz; ++i) {
        windowAccumulator[i] = 0.0f;
    }
    
    if (cd.accumulatorFill > si) {
        cd.accumulatorFill -= si;
    } else {
        cd.accumulatorFill = 0;
        if (cd.draining) {
            if (m_debugLevel > 1) {
                cerr << "RubberBandStretcher::Impl::processChunks: setting outputComplete to true" << endl;
            }
            cd.outputComplete = true;
        }
    }
}

void
RubberBandStretcher::Impl::writeOutput(RingBuffer<float> &to, float *from, size_t qty, size_t &outCount, size_t theoreticalOut)
{
    Profiler profiler("RubberBandStretcher::Impl::writeOutput");

    // In non-RT mode, we don't want to write the first startSkip
    // samples, because the first chunk is centred on the start of the
    // output.  In RT mode we didn't apply any pre-padding in
    // configure(), so we don't want to remove any here.

    size_t startSkip = 0;
    if (!m_realtime) {
        startSkip = lrintf((m_windowSize/2) / m_pitchScale);
    }

    if (outCount > startSkip) {
        
        // this is the normal case

        if (theoreticalOut > 0) {
            if (m_debugLevel > 1) {
                cerr << "theoreticalOut = " << theoreticalOut
                     << ", outCount = " << outCount
                     << ", startSkip = " << startSkip
                     << ", qty = " << qty << endl;
            }
            if (outCount - startSkip <= theoreticalOut &&
                outCount - startSkip + qty > theoreticalOut) {
                qty = theoreticalOut - (outCount - startSkip);
                if (m_debugLevel > 1) {
                    cerr << "reduce qty to " << qty << endl;
                }
            }
        }

        if (m_debugLevel > 2) {
            cerr << "writing " << qty << endl;
        }

        size_t written = to.write(from, qty);

        if (written < qty) {
            cerr << "WARNING: RubberBandStretcher::Impl::writeOutput: "
                 << "Buffer overrun on output: wrote " << written
                 << " of " << qty << " samples" << endl;
        }

        outCount += written;
        return;
    }

    // the rest of this is only used during the first startSkip samples

    if (outCount + qty <= startSkip) {
        if (m_debugLevel > 1) {
            cerr << "qty = " << qty << ", startSkip = "
                 << startSkip << ", outCount = " << outCount
                 << ", discarding" << endl;
        }
        outCount += qty;
        return;
    }

    size_t off = startSkip - outCount;
    if (m_debugLevel > 1) {
        cerr << "qty = " << qty << ", startSkip = "
             << startSkip << ", outCount = " << outCount
             << ", writing " << qty - off
             << " from start offset " << off << endl;
    }
    to.write(from + off, qty - off);
    outCount += qty;
}

int
RubberBandStretcher::Impl::available() const
{
    Profiler profiler("RubberBandStretcher::Impl::available");

    if (m_threaded) {
        MutexLocker locker(&m_threadSetMutex);
        if (m_channelData.empty()) return 0;
    } else {
        if (m_channelData.empty()) return 0;
    }

    if (!m_threaded) {
        for (size_t c = 0; c < m_channels; ++c) {
            if (m_channelData[c]->inputSize >= 0) {
//                cerr << "available: m_done true" << endl;
                if (m_channelData[c]->inbuf->getReadSpace() > 0) {
//                    cerr << "calling processChunks(" << c << ") from available" << endl;
                    //!!! do we ever actually do this? if so, this method should not be const
                    // ^^^ yes, we do sometimes -- e.g. when fed a very short file
                    bool any = false, last = false;
                    ((RubberBandStretcher::Impl *)this)->processChunks(c, any, last);
                }
            }
        }
    }

    size_t min = 0;
    bool consumed = true;
    bool haveResamplers = false;

    for (size_t i = 0; i < m_channels; ++i) {
        size_t availIn = m_channelData[i]->inbuf->getReadSpace();
        size_t availOut = m_channelData[i]->outbuf->getReadSpace();
        if (m_debugLevel > 2) {
            cerr << "available on channel " << i << ": " << availOut << " (waiting: " << availIn << ")" << endl;
        }
        if (i == 0 || availOut < min) min = availOut;
        if (!m_channelData[i]->outputComplete) consumed = false;
        if (m_channelData[i]->resampler) haveResamplers = true;
    }

    if (min == 0 && consumed) return -1;
    if (m_pitchScale == 1.0) return min;

    if (haveResamplers) return min; // resampling has already happened
    return int(floor(min / m_pitchScale));
}

size_t
RubberBandStretcher::Impl::retrieve(float *const *output, size_t samples) const
{
    Profiler profiler("RubberBandStretcher::Impl::retrieve");

    size_t got = samples;

    for (size_t c = 0; c < m_channels; ++c) {
        size_t gotHere = m_channelData[c]->outbuf->read(output[c], got);
        if (gotHere < got) {
            if (c > 0) {
                if (m_debugLevel > 0) {
                    cerr << "RubberBandStretcher::Impl::retrieve: WARNING: channel imbalance detected" << endl;
                }
            }
            got = gotHere;
        }
    }

    return got;
}

}

