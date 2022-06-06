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

R3StretcherImpl::R3StretcherImpl(Parameters parameters,
                                 double initialTimeRatio,
                                 double initialPitchScale) :
    m_parameters(parameters),
    m_timeRatio(initialTimeRatio),
    m_pitchScale(initialPitchScale),
    m_guide(Guide::Parameters(m_parameters.sampleRate, parameters.logger)),
    m_guideConfiguration(m_guide.getConfiguration()),
    m_channelAssembly(m_parameters.channels),
    m_troughPicker(m_guideConfiguration.classificationFftSize / 2 + 1),
    m_inhop(1),
    m_prevOuthop(1),
    m_draining(false)
{
    BinSegmenter::Parameters segmenterParameters
        (m_guideConfiguration.classificationFftSize,
         m_parameters.sampleRate);
    BinClassifier::Parameters classifierParameters
        (m_guideConfiguration.classificationFftSize / 2 + 1,
         9, 1, 10, 2.0, 2.0, 1.0e-7);

    int inRingBufferSize = m_guideConfiguration.longestFftSize * 2;
    int outRingBufferSize = m_guideConfiguration.longestFftSize * 16;

    for (int c = 0; c < m_parameters.channels; ++c) {
        m_channelData.push_back(std::make_shared<ChannelData>
                                (segmenterParameters,
                                 classifierParameters,
                                 m_guideConfiguration.longestFftSize,
                                 inRingBufferSize,
                                 outRingBufferSize));
        for (auto band: m_guideConfiguration.fftBandLimits) {
            int fftSize = band.fftSize;
            m_channelData[c]->scales[fftSize] =
                std::make_shared<ChannelScaleData>
                (fftSize, m_guideConfiguration.longestFftSize);
        }
    }
        
    for (auto band: m_guideConfiguration.fftBandLimits) {
        int fftSize = band.fftSize;
        GuidedPhaseAdvance::Parameters guidedParameters
            (fftSize, m_parameters.sampleRate, m_parameters.channels,
             m_parameters.logger);
        m_scaleData[fftSize] = std::make_shared<ScaleData>(guidedParameters);
    }

    m_calculator = std::unique_ptr<StretchCalculator>
        (new StretchCalculator(int(round(m_parameters.sampleRate)), //!!! which is a double...
                               1, false)); // no fixed inputIncrement

    Resampler::Parameters resamplerParameters;
    resamplerParameters.quality = Resampler::FastestTolerable;
    resamplerParameters.dynamism = Resampler::RatioOftenChanging;
    resamplerParameters.ratioChange = Resampler::SmoothRatioChange;
    resamplerParameters.initialSampleRate = m_parameters.sampleRate;
    resamplerParameters.maxBufferSize = m_guideConfiguration.longestFftSize; //!!!???
    m_resampler = std::unique_ptr<Resampler>
        (new Resampler(resamplerParameters, m_parameters.channels));

    m_formant = std::unique_ptr<FormantData>
        (new FormantData(m_guideConfiguration.classificationFftSize));
    
    calculateHop();

    m_prevInhop = m_inhop;
    m_prevOuthop = int(round(m_inhop * getEffectiveRatio()));

    if (!m_inhop.is_lock_free()) {
        m_parameters.logger("WARNING: std::atomic<int> is not lock-free");
    }
    if (!m_timeRatio.is_lock_free()) {
        m_parameters.logger("WARNING: std::atomic<double> is not lock-free");
    }
}

WindowType
R3StretcherImpl::ScaleData::analysisWindowShape(int fftSize)
{
//!!!    return HannWindow;
    if (fftSize == 4096) return HannWindow;
    else return NiemitaloForwardWindow;
}

int
R3StretcherImpl::ScaleData::analysisWindowLength(int fftSize)
{
    return fftSize;
}

WindowType
R3StretcherImpl::ScaleData::synthesisWindowShape(int fftSize)
{
//!!!    return HannWindow;
    if (fftSize == 4096) return HannWindow;
    else return NiemitaloReverseWindow;
}

int
R3StretcherImpl::ScaleData::synthesisWindowLength(int fftSize)
{
//!!!    return fftSize/2;
    if (fftSize == 4096) return fftSize/2;
    else return fftSize;
}

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
R3StretcherImpl::setFormantOption(RubberBandStretcher::Options options)
{
    int mask = (RubberBandStretcher::OptionFormantShifted |
                RubberBandStretcher::OptionFormantPreserved);
    m_parameters.options &= ~mask;
    options &= mask;
    m_parameters.options |= options;
}

void
R3StretcherImpl::setPitchOption(RubberBandStretcher::Options options)
{
    int mask = (RubberBandStretcher::OptionPitchHighQuality |
                RubberBandStretcher::OptionPitchHighSpeed |
                RubberBandStretcher::OptionPitchHighConsistency);
    m_parameters.options &= ~mask;
    options &= mask;
    m_parameters.options |= options;
}

void
R3StretcherImpl::calculateHop()
{
    double ratio = getEffectiveRatio();
    double proposedOuthop = 256.0;
//!!!    if (proposedOuthop * m_pitchScale > 2048.0) {
//        proposedOuthop = 2048.0 / m_pitchScale;
//    }
    double inhop = 1.0;
    
    if (ratio > 1.0) {
        inhop = proposedOuthop / ratio;
        if (inhop < 1.0) {
            m_parameters.logger("WARNING: Extreme ratio yields ideal inhop < 1, results may be suspect");
            inhop = 1.0;
        }
    } else {
        inhop = std::min(proposedOuthop / ratio, 340.0);
    }

    m_inhop = int(round(inhop));

    //!!! but if we now have outhop > 4096 ever, we will crash, so we must check

//    std::cout << "R3StretcherImpl::calculateHop: inhop = " << m_inhop << ", proposed outhop = " << proposedOuthop << ", mean outhop = " << m_inhop * ratio << std::endl;
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
    m_calculator->reset();
    m_resampler->reset();

    for (auto &it : m_scaleData) {
        it.second->guided.reset();
    }

    for (auto &cd : m_channelData) {
        cd->reset();
    }

    m_prevInhop = m_inhop;
    m_prevOuthop = int(round(m_inhop * getEffectiveRatio()));
}

size_t
R3StretcherImpl::getSamplesRequired() const
{
    if (available() != 0) return 0;
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
            auto newBuf = m_channelData[c]->inbuf->resized(newSize);
            m_channelData[c]->inbuf = std::unique_ptr<RingBuffer<float>>(newBuf);
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
    int longest = m_guideConfiguration.longestFftSize;
    int channels = m_parameters.channels;

//    m_calculator->setDebugLevel(3);

    int inhop = m_inhop;

    double effectivePitchRatio = 1.0 / m_pitchScale;
    if (m_resampler) {
        effectivePitchRatio = m_resampler->getEffectiveRatio(effectivePitchRatio);
    }
    
    int outhop = m_calculator->calculateSingle(m_timeRatio,
                                               effectivePitchRatio,
                                               1.f,
                                               inhop,
                                               longest,
                                               longest);

    // Now inhop is the distance by which the input stream will be
    // advanced after our current frame has been read, and outhop is
    // the distance by which the output will be advanced after it has
    // been emitted; m_prevInhop and m_prevOuthop are the
    // corresponding values the last time a frame was processed (*not*
    // just the last time this function was called, since we can
    // return without doing anything if the output buffer is full).
    //
    // Our phase adjustments need to be based on the distances we have
    // advanced the input and output since the previous frame, not the
    // distances we are about to advance them, so they use the m_prev
    // values.
/*
    if (inhop != m_prevInhop) {
        std::cout << "Note: inhop has changed from " << m_prevInhop
                  << " to " << inhop << std::endl;
    }
    if (outhop != m_prevOuthop) {
        std::cout << "Note: outhop has changed from " << m_prevOuthop
                  << " to " << outhop << std::endl;
    }
*/    
    while (m_channelData.at(0)->outbuf->getWriteSpace() >= outhop) {

        // NB our ChannelData, ScaleData, and ChannelScaleData maps
        // contain shared_ptrs; whenever we retain one of them in a
        // variable, we do so by reference to avoid copying the
        // shared_ptr (as that is not realtime safe). Same goes for
        // the map iterators

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

        // Analysis
        
        for (int c = 0; c < channels; ++c) {
            analyseChannel(c, inhop, m_prevInhop, m_prevOuthop);
        }

        if (m_parameters.options & RubberBandStretcher::OptionFormantPreserved) {
            m_formant->enabled = true;
            analyseFormant();
        } else {
            m_formant->enabled = false;
        }
        
        // Phase update. This is synchronised across all channels
        
        for (auto &it : m_channelData[0]->scales) {
            int fftSize = it.first;
            for (int c = 0; c < channels; ++c) {
                auto &cd = m_channelData.at(c);
                auto &scale = cd->scales.at(fftSize);
                m_channelAssembly.mag[c] = scale->mag.data();
                m_channelAssembly.phase[c] = scale->phase.data();
                m_channelAssembly.guidance[c] = &cd->guidance;
                m_channelAssembly.outPhase[c] = scale->advancedPhase.data();
            }
            m_scaleData.at(fftSize)->guided.advance
                (m_channelAssembly.outPhase.data(),
                 m_channelAssembly.mag.data(),
                 m_channelAssembly.phase.data(),
                 m_guideConfiguration,
                 m_channelAssembly.guidance.data(),
                 m_prevInhop,
                 m_prevOuthop);
        }

        // Resynthesis
        
        for (int c = 0; c < channels; ++c) {
            synthesiseChannel(c, outhop);
        }

        // Resample

        int resampledCount = 0;
        if (m_resampler) {
            for (int c = 0; c < channels; ++c) {
                auto &cd = m_channelData.at(c);
                m_channelAssembly.mixdown[c] = cd->mixdown.data();
                m_channelAssembly.resampled[c] = cd->resampled.data();
            }
            resampledCount = m_resampler->resample
                (m_channelAssembly.resampled.data(),
                 m_channelData[0]->resampled.size(),
                 m_channelAssembly.mixdown.data(),
                 outhop,
                 1.0 / m_pitchScale,
                 m_draining && readSpace < longest);
        }

        // Emit

        for (int c = 0; c < channels; ++c) {
            auto &cd = m_channelData.at(c);
            if (m_resampler) {
                cd->outbuf->write(cd->resampled.data(), resampledCount);
            } else {
                cd->outbuf->write(cd->mixdown.data(), outhop);
            }
        
            int readSpace = cd->inbuf->getReadSpace();
            if (readSpace < inhop) {
                // This should happen only when draining
                cd->inbuf->skip(readSpace);
            } else {
                cd->inbuf->skip(inhop);
            }
        }

        m_prevInhop = inhop;
        m_prevOuthop = outhop;
    }
}

void
R3StretcherImpl::analyseChannel(int c, int inhop, int prevInhop, int prevOuthop)
{
    int longest = m_guideConfiguration.longestFftSize;
    int classify = m_guideConfiguration.classificationFftSize;

    auto &cd = m_channelData.at(c);
    double *buf = cd->scales.at(longest)->timeDomain.data();

    int readSpace = cd->inbuf->getReadSpace();
    if (readSpace < longest) {
        cd->inbuf->peek(buf, readSpace);
        v_zero(buf + readSpace, longest - readSpace);
    } else {
        cd->inbuf->peek(buf, longest);
    }
    
    // We have a single unwindowed frame at the longest FFT size
    // ("scale"). Populate the shorter FFT sizes from the centre of
    // it, windowing as we copy. The classification scale is handled
    // separately because it has readahead, so skip it here as well as
    // the longest. (In practice this means we are probably only
    // populating one scale)

    for (auto &it: cd->scales) {
        int fftSize = it.first;
        if (fftSize == classify || fftSize == longest) continue;
        int offset = (longest - fftSize) / 2;
        m_scaleData.at(fftSize)->analysisWindow.cut
            (buf + offset, it.second->timeDomain.data());
    }

    // The classification scale has a one-hop readahead, so populate
    // the readahead from further down the long unwindowed frame.

    auto &classifyScale = cd->scales.at(classify);
    ClassificationReadaheadData &readahead = cd->readahead;

    m_scaleData.at(classify)->analysisWindow.cut
        (buf + (longest - classify) / 2 + inhop,
         readahead.timeDomain.data());

    // If inhop has changed since the previous frame, we'll have to
    // populate the classification scale (but for analysis/resynthesis
    // rather than classification) anew rather than reuse the previous
    // readahead. Pity...

    bool haveValidReadahead = cd->haveReadahead;
    if (inhop != prevInhop) haveValidReadahead = false;

    if (!haveValidReadahead) {
        m_scaleData.at(classify)->analysisWindow.cut
            (buf + (longest - classify) / 2,
             classifyScale->timeDomain.data());
    }
            
    // Finally window the longest scale
    m_scaleData.at(longest)->analysisWindow.cut(buf);

    // FFT shift, forward FFT, and carry out cartesian-polar
    // conversion for each FFT size.

    // For the classification scale we need magnitudes for the full
    // range (polar only in a subset) and we operate in the readahead,
    // pulling current values from the existing readahead (except
    // where the inhop has changed as above, in which case we need to
    // do both readahead and current)

    v_fftshift(readahead.timeDomain.data(), classify);

    if (haveValidReadahead) {
        v_copy(classifyScale->mag.data(),
               readahead.mag.data(),
               classifyScale->bufSize);
        v_copy(classifyScale->phase.data(),
               readahead.phase.data(),
               classifyScale->bufSize);
    }

    m_scaleData.at(classify)->fft.forward(readahead.timeDomain.data(),
                                          classifyScale->real.data(),
                                          classifyScale->imag.data());

    for (const auto &b : m_guideConfiguration.fftBandLimits) {
        if (b.fftSize == classify) {
            if (b.b0min > 0) {
                v_cartesian_to_magnitudes(readahead.mag.data(),
                                          classifyScale->real.data(),
                                          classifyScale->imag.data(),
                                          b.b0min);
            }
                    
            v_cartesian_to_polar(readahead.mag.data() + b.b0min,
                                 readahead.phase.data() + b.b0min,
                                 classifyScale->real.data() + b.b0min,
                                 classifyScale->imag.data() + b.b0min,
                                 b.b1max - b.b0min);
                    
            if (b.b1max < classify/2 + 1) {
                v_cartesian_to_magnitudes
                    (readahead.mag.data() + b.b1max,
                     classifyScale->real.data() + b.b1max,
                     classifyScale->imag.data() + b.b1max,
                     classify/2 + 1 - b.b1max);
            }
                    
            v_scale(classifyScale->mag.data(),
                    1.0 / double(classify),
                    classifyScale->mag.size());
            break;
        }
    }

    cd->haveReadahead = true;

    // For the others (and the classify as well, if the inhop has
    // changed or we haven't filled the readahead yet) we operate
    // directly in the scale data and restrict the range for
    // cartesian-polar conversion
            
    for (auto &it: cd->scales) {
        int fftSize = it.first;
        if (fftSize == classify && haveValidReadahead) continue;
        auto &scale = it.second;

        v_fftshift(scale->timeDomain.data(), fftSize);

        m_scaleData.at(fftSize)->fft.forward(scale->timeDomain.data(),
                                             scale->real.data(),
                                             scale->imag.data());

        //!!! This should be a map
        for (const auto &b : m_guideConfiguration.fftBandLimits) {
            if (b.fftSize == fftSize) {
                v_cartesian_to_polar(scale->mag.data() + b.b0min,
                                     scale->phase.data() + b.b0min,
                                     scale->real.data() + b.b0min,
                                     scale->imag.data() + b.b0min,
                                     b.b1max - b.b0min);
                v_scale(scale->mag.data() + b.b0min,
                        1.0 / double(fftSize),
                        b.b1max - b.b0min);
                break;
            }
        }
    }

    // Use the classification scale to get a bin segmentation and
    // calculate the adaptive frequency guide for this channel

    v_copy(cd->classification.data(), cd->nextClassification.data(),
           cd->classification.size());
    cd->classifier->classify(readahead.mag.data(),
                             cd->nextClassification.data());

    cd->prevSegmentation = cd->segmentation;
    cd->segmentation = cd->nextSegmentation;
    cd->nextSegmentation = cd->segmenter->segment(cd->nextClassification.data());

    m_troughPicker.findNearestAndNextPeaks
        (classifyScale->mag.data(), 3, nullptr,
         classifyScale->troughs.data());
            
    double instantaneousRatio = double(prevOuthop) / double(prevInhop);
//!!!???    bool specialCaseUnity = !(m_parameters.options &
//                              RubberBandStretcher::OptionPitchHighConsistency);

    bool specialCaseUnity = true;
        
    m_guide.calculate(instantaneousRatio,
                      classifyScale->mag.data(),
                      classifyScale->troughs.data(),
                      classifyScale->prevMag.data(),
                      cd->segmentation,
                      cd->prevSegmentation,
                      cd->nextSegmentation,
                      specialCaseUnity,
                      cd->guidance);
}

void
R3StretcherImpl::analyseFormant()
{
    int classify = m_guideConfiguration.classificationFftSize;
    int binCount = classify/2 + 1;
    int channels = m_parameters.channels;

    auto &f = *m_formant;

    v_zero(f.envelope.data(), binCount);
    
    for (int c = 0; c < channels; ++c) {
        auto &cd = m_channelData.at(c);
        auto &scale = cd->scales.at(classify);
        for (int i = 0; i < binCount; ++i) {
            f.envelope.at(i) += scale->mag.at(i) / double(channels);
        }
    }

    m_scaleData.at(classify)->fft.inverseCepstral
        (f.envelope.data(), f.cepstra.data());
    
    int cutoff = int(floor(m_parameters.sampleRate / 650.0));
    if (cutoff < 1) cutoff = 1;

    f.cepstra[0] /= 2.0;
    f.cepstra[cutoff-1] /= 2.0;
    for (int i = cutoff; i < classify; ++i) {
        f.cepstra[i] = 0.0;
    }
    v_scale(f.cepstra.data(), 1.0 / double(classify), cutoff);

    m_scaleData.at(classify)->fft.forward
        (f.cepstra.data(), f.envelope.data(),
         f.spare.data()); // shifted is just a spare for this one

    v_exp(f.envelope.data(), binCount);
    v_square(f.envelope.data(), binCount);

    for (int i = 0; i < binCount; ++i) {
        if (f.envelope[i] > 1.0e10) f.envelope[i] = 1.0e10;
    }
}

void
R3StretcherImpl::synthesiseChannel(int c, int outhop)
{
    int longest = m_guideConfiguration.longestFftSize;

    auto &cd = m_channelData.at(c);

    for (auto &it : cd->scales) {

        auto &scale = it.second;
        int bufSize = scale->bufSize;

        // copy to prevMag before filtering
        v_copy(scale->prevMag.data(),
               scale->mag.data(),
               bufSize);
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

        // The frequency filter is applied naively in the frequency
        // domain. Aliasing is reduced by the shorter resynthesis
        // window

        //!!! I don't think we have binForFrequency etc available in
        //!!! this class - really that's ridiculous
        
        int lowBin = int(floor(fftSize * band.f0 / m_parameters.sampleRate));
        int highBin = int(floor(fftSize * band.f1 / m_parameters.sampleRate));
        if (highBin % 2 == 1) --highBin;

        int formantHigh = int(floor(fftSize * 10000.0 / m_parameters.sampleRate));
        for (int i = 0; i < lowBin; ++i) {
            scale->mag[i] = 0.0;
        }
        if (m_formant->enabled) {
            double targetFactor = double(m_formant->fftSize) / double(fftSize);
            double sourceFactor = targetFactor * m_pitchScale;
            double maxRatio = 60.0;
            double minRatio = 1.0 / maxRatio;
            for (int i = lowBin; i < highBin && i < formantHigh; ++i) {
                double source = m_formant->envelopeAt(i * sourceFactor);
                double target = m_formant->envelopeAt(i * targetFactor);
                if (target > 0.0) {
                    double ratio = source / target;
                    if (ratio < minRatio) ratio = minRatio;
                    if (ratio > maxRatio) ratio = maxRatio;
                    scale->mag[i] *= ratio;
                }
            }
        }
        for (int i = lowBin; i < highBin; ++i) {
            scale->mag[i] *= winscale;
        }
        for (int i = highBin; i < fftSize/2 + 1; ++i) {
            scale->mag[i] = 0.0;
        }
    }

    // Resynthesise each FFT size (scale) individually, then sum. This
    // is easier to manage scaling for in situations with a varying
    // resynthesis hop
            
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
                     scale->advancedPhase.data() + offset,
                     b.b1max - offset);
                break;
            }
        }

        scaleData->fft.inverse(scale->real.data(),
                               scale->imag.data(),
                               scale->timeDomain.data());

        v_fftshift(scale->timeDomain.data(), fftSize);

        // Synthesis window may be shorter than analysis window, so
        // copy and cut only from the middle of the time-domain frame;
        // and the accumulator length always matches the longest FFT
        // size, so as to make mixing straightforward, so there is an
        // additional offset needed for the target
                
        int synthesisWindowSize = scaleData->synthesisWindow.getSize();
        int fromOffset = (fftSize - synthesisWindowSize) / 2;
        int toOffset = (longest - synthesisWindowSize) / 2;

        scaleData->synthesisWindow.cutAndAdd
            (scale->timeDomain.data() + fromOffset,
             scale->accumulator.data() + toOffset);
    }

    // Mix this channel and move the accumulator along
            
    float *mixptr = cd->mixdown.data();
    v_zero(mixptr, outhop);

    for (auto &it : cd->scales) {
        auto &scale = it.second;

        double *accptr = scale->accumulator.data();
        for (int i = 0; i < outhop; ++i) {
            mixptr[i] += float(accptr[i]);
        }

        int n = scale->accumulator.size() - outhop;
        v_move(accptr, accptr + outhop, n);
        v_zero(accptr + n, outhop);
    }
}

}

