/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band Library
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2022 Particular Programs Ltd.

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

#include "R3StretcherImpl.h"

#include <array>

namespace RubberBand {

void
R3StretcherImpl::setTimeRatio(double ratio)
{
    m_timeRatio = ratio;
}

void
R3StretcherImpl::setPitchScale(double scale)
{
    m_pitchScale = scale;
}

double
R3StretcherImpl::getTimeRatio() const
{
    return m_timeRatio;
}

double
R3StretcherImpl::getPitchScale() const
{
    return m_pitchScale;
}

size_t
R3StretcherImpl::getLatency() const
{
    return 0; //!!!
}

size_t
R3StretcherImpl::getChannelCount() const
{
    return m_parameters.channels;
}

void
R3StretcherImpl::reset()
{
    //!!!
}

size_t
R3StretcherImpl::getSamplesRequired() const
{
    int longest = m_guideConfiguration.longestFftSize;
    size_t rs = m_channelData[0]->inbuf->getReadSpace();
    if (rs < longest) {
        return longest - rs;
    } else {
        return 0;
    }
}

void
R3StretcherImpl::process(const float *const *input, size_t samples, bool final)
{
    //!!! todo: final

    m_parameters.logger("process called");
    if (final) {
        m_parameters.logger("final = true");
        m_draining = true;
    }
    
    bool allConsumed = false;

    size_t ws = m_channelData[0]->inbuf->getWriteSpace();
    if (samples > ws) {
        //!!! check this
        m_parameters.logger("R3StretcherImpl::process: WARNING: Forced to increase input buffer size. Either setMaxProcessSize was not properly called or process is being called repeatedly without retrieve.");
        size_t newSize = m_channelData[0]->inbuf->getSize() - ws + samples;
        for (int c = 0; c < m_parameters.channels; ++c) {
            m_channelData[c]->inbuf =
                std::unique_ptr<RingBuffer<float>>
                (m_channelData[c]->inbuf->resized(newSize));
        }
    }

    for (int c = 0; c < m_parameters.channels; ++c) {
        m_channelData[c]->inbuf->write(input[c], samples);
    }

    consume();
}

int
R3StretcherImpl::available() const
{
    m_parameters.logger("available called");
    int av = int(m_channelData[0]->outbuf->getReadSpace());
    if (av == 0 && m_draining) return -1;
    else return av;
}

size_t
R3StretcherImpl::retrieve(float *const *output, size_t samples) const
{
    m_parameters.logger("retrieve called");
    size_t got = samples;
    
    for (size_t c = 0; c < m_parameters.channels; ++c) {
        size_t gotHere = m_channelData[c]->outbuf->read(output[c], got);
        if (gotHere < got) {
            if (c > 0) {
                m_parameters.logger("R3StretcherImpl::retrieve: WARNING: channel imbalance detected");
            }
            got = gotHere;
        }
    }

    return got;
}

void
R3StretcherImpl::consume()
{
    int inhop = 171, outhop = 256; //!!!
    double ratio = double(outhop) / double(inhop);

    int longest = m_guideConfiguration.longestFftSize;
    int classify = m_guideConfiguration.classificationFftSize;
    
    while ((m_draining || m_channelData[0]->inbuf->getReadSpace() >= longest) &&
           m_channelData[0]->outbuf->getWriteSpace() >= outhop) {

        m_parameters.logger("consume looping");
        
        for (int c = 0; c < m_parameters.channels; ++c) {

            auto cd = m_channelData[c];
            auto longestScale = cd->scales.at(longest);
            
            cd->inbuf->peek(longestScale->timeDomainFrame.data(), longest);

            for (auto it: cd->scales) {
                int fftSize = it.first;
                auto scale = it.second;
                if (fftSize == longest) continue;
                int offset = (longest - fftSize) / 2;
                m_scaleData.at(fftSize)->analysisWindow.cut
                    (longestScale->timeDomainFrame.data() + offset,
                     scale->timeDomainFrame.data());
            }

            m_scaleData.at(longest)->analysisWindow.cut
                (longestScale->timeDomainFrame.data());
        }

        for (int c = 0; c < m_parameters.channels; ++c) {

            auto cd = m_channelData[c];

            //!!! There are some aspects of scaling etc handled in bsq
            //!!! that are not yet here

            for (auto it: cd->scales) {
                int fftSize = it.first;
                auto scale = it.second;
                m_scaleData.at(fftSize)->fft.forwardPolar
                    (scale->timeDomainFrame.data(),
                     scale->mag.data(),
                     scale->phase.data());
            }
        }

        for (int c = 0; c < m_parameters.channels; ++c) {
            auto cd = m_channelData[c];
            auto classifyScale = cd->scales.at(classify);
            cd->prevSegmentation = cd->segmentation;
            cd->segmentation = cd->segmenter->segment(classifyScale->mag.data());
            m_troughPicker.findNearestAndNextPeaks
                (classifyScale->mag.data(), 3, nullptr,
                 classifyScale->nextTroughs.data());
            m_guide.calculate(ratio, classifyScale->mag.data(),
                              classifyScale->nextTroughs.data(),
                              classifyScale->prevMag.data(),
                              cd->segmentation,
                              cd->prevSegmentation,
                              BinSegmenter::Segmentation(), //!!!
                              cd->guidance);
        }

        for (auto it : m_channelData[0]->scales) {
            int fftSize = it.first;
            for (int c = 0; c < m_parameters.channels; ++c) {
                auto cd = m_channelData[c];
                auto classifyScale = cd->scales.at(fftSize);
                m_channelAssembly.mag[c] = classifyScale->mag.data();
                m_channelAssembly.phase[c] = classifyScale->phase.data();
                m_channelAssembly.guidance[c] = &cd->guidance;
                m_channelAssembly.outPhase[c] = classifyScale->outPhase.data();
            }
            m_scaleData.at(fftSize)->guided.advance
                (m_channelAssembly.outPhase.data(),
                 m_channelAssembly.mag.data(),
                 m_channelAssembly.phase.data(),
                 m_guideConfiguration,
                 m_channelAssembly.guidance.data(),
                 inhop,
                 outhop);
        }

        for (int c = 0; c < m_parameters.channels; ++c) {
            for (auto it : m_channelData[c]->scales) {
                auto scale = it.second;
                int bufSize = scale->bufSize;
                // copy to prevMag before filtering
                v_copy(scale->prevMag.data(), scale->mag.data(), bufSize);
                v_copy(scale->prevOutPhase.data(), scale->outPhase.data(), bufSize);
                //!!! seems wasteful
                for (int i = 0; i < bufSize; ++i) {
                    scale->phase[i] = princarg(scale->outPhase[i]);
                }
            }
        }
        
        //!!! + filter here
        
        for (int c = 0; c < m_parameters.channels; ++c) {
            for (auto it : m_channelData[c]->scales) {
                int fftSize = it.first;
                auto scale = it.second;
                auto scaleData = m_scaleData.at(fftSize);
                int bufSize = scale->bufSize;
                scaleData->fft.inversePolar(scale->mag.data(),
                                            scale->phase.data(),
                                            scale->timeDomainFrame.data());
                /*
                int synthesisWindowSize = scaleData->synthesisWindow.getSize();
                int fromOffset = (fftSize - synthesisWindowSize) / 2;
                int toOffset = (longest - synthesisWindowSize) / 2;
                //!!! not right - accumulator is of scale data size, not full longest size - we need offset when mixing into mixdown buffer below as well
                scaleData->synthesisWindow.cutAndAdd
                    (scale->timeDomainFrame.data() + fromOffset,
                     scale->accumulator.data() + toOffset);
                */
                int synthesisWindowSize = scaleData->synthesisWindow.getSize();
                int offset = (fftSize - synthesisWindowSize) / 2;
                scaleData->synthesisWindow.cutAndAdd
                    (scale->timeDomainFrame.data() + offset,
                     scale->accumulator.data());
            }
        }

        for (int c = 0; c < m_parameters.channels; ++c) {
            auto cd = m_channelData[c];
            v_zero(cd->mixdown.data(), outhop);
            for (auto it : cd->scales) {
                auto scale = it.second;
                auto &acc = scale->accumulator;
                v_add(cd->mixdown.data(), acc.data(), outhop);
                int n = acc.size() - outhop;
                v_move(acc.data(), acc.data() + outhop, n);
                v_zero(acc.data() + n, outhop);
            }
            m_channelData[c]->outbuf->write(cd->mixdown.data(), outhop);
            m_channelData[c]->inbuf->skip(inhop);
        }
    }
    
}


}

