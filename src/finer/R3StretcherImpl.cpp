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

#include "../common/VectorOpsComplex.h"

#include <array>

namespace RubberBand {

void
R3StretcherImpl::setTimeRatio(double ratio)
{
    m_timeRatio = ratio;
    calculateHop();
}

void
R3StretcherImpl::setPitchScale(double scale)
{
    m_pitchScale = scale;
    calculateHop();
}

void
R3StretcherImpl::calculateHop()
{
    double ratio = getEffectiveRatio();
    double proposedOuthop = 256;
    
    if (ratio > 1.0) {
        double inhop = proposedOuthop / ratio;
        if (inhop < 1.0) {
            m_parameters.logger("WARNING: Extreme ratio yields ideal inhop < 1, results may be suspect");
            m_inhop = 1;
        } else {
            m_inhop = int(round(inhop));
        }
    } else {
        double inhop = std::min(proposedOuthop / ratio, 340.0);
        m_inhop = int(round(inhop));
    }

    m_prevOuthop = int(round(m_inhop * ratio));
    
    std::ostringstream str;
    str << "R3StretcherImpl::calculateHop: for effective ratio " << ratio
        << " calculated (typical) inhop of " << m_inhop << std::endl;
    m_parameters.logger(str.str());
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

//!!! __attribute__((annotate("realtime")))
void
R3StretcherImpl::process(const float *const *input, size_t samples, bool final)
{
    //!!! todo: final

//!!!    m_parameters.logger("process called");
    if (final) {
//        m_parameters.logger("final = true");
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

//!!! __attribute__((annotate("realtime")))
int
R3StretcherImpl::available() const
{
//!!!    m_parameters.logger("available called");
    int av = int(m_channelData[0]->outbuf->getReadSpace());
    if (av == 0 && m_draining) return -1;
    else return av;
}

//!!! __attribute__((annotate("realtime")))
size_t
R3StretcherImpl::retrieve(float *const *output, size_t samples) const
{
//!!!    m_parameters.logger("retrieve called");
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
    double ratio = getEffectiveRatio();

    int longest = m_guideConfiguration.longestFftSize;
    int classify = m_guideConfiguration.classificationFftSize;

    m_calculator->setDebugLevel(3);
    
    int outhop = m_calculator->calculateSingle(ratio,
                                               1.0 / m_pitchScale,
                                               1.f,
                                               m_inhop,
                                               longest,
                                               longest);

    std::cout << "outhop = " << outhop << std::endl;

    double instantaneousRatio = double(m_prevOuthop) / double(m_inhop);
    m_prevOuthop = outhop;

    while (m_channelData.at(0)->outbuf->getWriteSpace() >= outhop) {

        int readSpace = m_channelData.at(0)->inbuf->getReadSpace();
        if (readSpace < longest) {
            if (m_draining) {
                if (readSpace == 0) {
                    break;
                }
            } else {
                break;
            }
        }

        // Analysis. This is per channel
        
        for (int c = 0; c < m_parameters.channels; ++c) {

            // Our ChannelData, ScaleData, and ChannelScaleData maps
            // contain shared_ptrs; whenever we retain one of them in
            // a variable here, we do so by reference to avoid copying
            // the shared_ptr (as that is not realtime safe). Same
            // goes for the map iterators
            
            auto &cd = m_channelData.at(c);
            auto &longestScale = cd->scales.at(longest);
            double *buf = longestScale->timeDomain.data();

            if (readSpace < longest) {
                cd->inbuf->peek(buf, readSpace);
                v_zero(buf + readSpace, longest - readSpace);
            } else {
                cd->inbuf->peek(buf, longest);
            }

            // We have a single unwindowed frame at the longest FFT
            // size ("scale"). Populate the shorter FFT sizes from the
            // centre of it, windowing as we copy
            for (auto &it: cd->scales) {
                int fftSize = it.first;
                auto &scale = it.second;
                if (fftSize == longest) continue;
                int offset = (longest - fftSize) / 2;
                m_scaleData.at(fftSize)->analysisWindow.cut
                    (buf + offset, scale->timeDomain.data());
            }

            // Then window the longest one
            m_scaleData.at(longest)->analysisWindow.cut(buf);

            // FFT shift, forward FFT, and carry out cartesian-polar
            // conversion for each FFT size
            for (auto &it: cd->scales) {
                int fftSize = it.first;
                auto &scale = it.second;
                
                v_fftshift(scale->timeDomain.data(), fftSize);

                if (fftSize == m_guideConfiguration.classificationFftSize) {

                    // For the classification scale we need the full range
                    m_scaleData.at(fftSize)->fft.forwardPolar
                        (scale->timeDomain.data(),
                         scale->mag.data(),
                         scale->phase.data());

                } else {

                    // For other scales we only need do
                    // cartesian-polar conversion for the necessary
                    // frequency subset
                    m_scaleData.at(fftSize)->fft.forward
                    (scale->timeDomain.data(),
                     scale->real.data(),
                     scale->imag.data());
                
                    for (const auto &b : m_guideConfiguration.fftBandLimits) {
                        if (b.fftSize == fftSize) {
                            int offset = b.b0min;
                            v_cartesian_to_polar
                                (scale->mag.data() + offset,
                                 scale->phase.data() + offset,
                                 scale->real.data() + offset,
                                 scale->imag.data() + offset,
                                 b.b1max - offset);
                            break;
                        }
                    }
                }
                
                v_scale(scale->mag.data(), 1.0 / double(fftSize),
                        scale->mag.size());
            }

            // Use the classification scale to get a bin segmentation
            // and calculate the adaptive frequency guide for this
            // channel
            auto &classifyScale = cd->scales.at(classify);
            cd->prevSegmentation = cd->segmentation;
            cd->segmentation =
                cd->segmenter->segment(classifyScale->mag.data());
            m_troughPicker.findNearestAndNextPeaks
                (classifyScale->mag.data(), 3, nullptr,
                 classifyScale->nextTroughs.data());
            m_guide.calculate(instantaneousRatio,
                              classifyScale->mag.data(),
                              classifyScale->nextTroughs.data(),
                              classifyScale->prevMag.data(),
                              cd->segmentation,
                              cd->prevSegmentation,
                              BinSegmenter::Segmentation(), //!!!
                              cd->guidance);
        }

        // Phase update. This is synchronised across all channels
        
        for (auto &it : m_channelData[0]->scales) {
            int fftSize = it.first;
            for (int c = 0; c < m_parameters.channels; ++c) {
                auto &cd = m_channelData.at(c);
                auto &classifyScale = cd->scales.at(fftSize);
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
                 m_inhop,
                 outhop);
        }

        // Resynthesis. This is per channel
        
        for (int c = 0; c < m_parameters.channels; ++c) {

            auto &cd = m_channelData.at(c);

            for (auto &it : cd->scales) {
                auto &scale = it.second;
                int bufSize = scale->bufSize;
                // copy to prevMag before filtering
                v_copy(scale->prevMag.data(), scale->mag.data(), bufSize);
                v_copy(scale->prevOutPhase.data(), scale->outPhase.data(), bufSize);
            }
        
            for (const auto &band : cd->guidance.fftBands) {
                int fftSize = band.fftSize;
                auto &scale = cd->scales.at(fftSize);
                auto &scaleData = m_scaleData.at(fftSize);

                //!!! messy and slow, but leave it until we've
                //!!! discovered whether we need a window accumulator
                //!!! (we probably do)
                int analysisWindowSize = scaleData->analysisWindow.getSize();
                int synthesisWindowSize = scaleData->synthesisWindow.getSize();
                int offset = (analysisWindowSize - synthesisWindowSize) / 2;
                double winscale = 0.0;
                for (int i = 0; i < synthesisWindowSize; ++i) {
                    winscale += scaleData->analysisWindow.getValue(i + offset) *
                        scaleData->synthesisWindow.getValue(i);
                }
                winscale = double(outhop) / winscale;

                // The frequency filter is applied naively in the
                // frequency domain. Aliasing is reduced by the
                // shorter resynthesis window
                
                double factor = m_parameters.sampleRate / double(fftSize);
                for (int i = 0; i < fftSize/2 + 1; ++i) {
                    double f = double(i) * factor;
                    if (f >= band.f0 && f < band.f1) {
                        //!!! check the mod 2 bit from stretch-fn
                        scale->mag[i] *= winscale;
                    } else {
                        scale->mag[i] = 0.f;
                    }
                }
            }

            // Resynthesise each FFT size (scale) individually, then
            // sum. This is easier to manage scaling for in situations
            // with a varying resynthesis hop
            
            for (auto &it : cd->scales) {
                int fftSize = it.first;
                auto &scale = it.second;
                auto &scaleData = m_scaleData.at(fftSize);
                
                for (const auto &b : m_guideConfiguration.fftBandLimits) {
                    if (b.fftSize == fftSize) {
                        int offset = b.b0min;
                        v_zero(scale->real.data(), fftSize/2 + 1);
                        v_zero(scale->imag.data(), fftSize/2 + 1);
                        v_polar_to_cartesian
                            (scale->real.data() + offset,
                             scale->imag.data() + offset,
                             scale->mag.data() + offset,
                             scale->outPhase.data() + offset,
                             b.b1max - offset);
                        break;
                    }
                }

                scaleData->fft.inverse(scale->real.data(),
                                       scale->imag.data(),
                                       scale->timeDomain.data());

                v_fftshift(scale->timeDomain.data(), fftSize);

                // Synthesis window is shorter than analysis window,
                // so copy and cut only from the middle of the
                // time-domain frame; and the accumulator length
                // always matches the longest FFT size, so as to make
                // mixing straightforward, so there is an additional
                // offset needed for the target
                
                int synthesisWindowSize = scaleData->synthesisWindow.getSize();
                int fromOffset = (fftSize - synthesisWindowSize) / 2;
                int toOffset = (m_guideConfiguration.longestFftSize -
                                synthesisWindowSize) / 2;

                scaleData->synthesisWindow.cutAndAdd
                    (scale->timeDomain.data() + fromOffset,
                     scale->accumulator.data() + toOffset);
            }

            // Mix and emit this channel
            
            double *mixptr = cd->mixdown.data();
            v_zero(mixptr, outhop);

            for (auto &it : cd->scales) {
                auto &scale = it.second;
                v_add(mixptr, scale->accumulator.data(), outhop);
            }

            cd->outbuf->write(mixptr, outhop);

            for (auto &it : cd->scales) {
                int fftSize = it.first;
                auto &scale = it.second;
                double *accptr = scale->accumulator.data();

                int n = scale->accumulator.size() - outhop;
                v_move(accptr, accptr + outhop, n);
                v_zero(accptr + n, outhop);
            }
            
            if (readSpace < m_inhop) {
                // This should happen only when draining
                cd->inbuf->skip(readSpace);
            } else {
                cd->inbuf->skip(m_inhop);
            }
        }
    }
}


}

