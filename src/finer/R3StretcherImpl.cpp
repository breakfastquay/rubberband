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
    m_formantScale(0.0),
    m_guide(Guide::Parameters(m_parameters.sampleRate, parameters.logger)),
    m_guideConfiguration(m_guide.getConfiguration()),
    m_channelAssembly(m_parameters.channels),
    m_inhop(1),
    m_prevOuthop(1),
    m_startSkip(0),
    m_studyInputDuration(0),
    m_totalTargetDuration(0),
    m_processInputDuration(0),
    m_totalOutputDuration(0),
    m_mode(ProcessMode::JustCreated)
{
    double maxClassifierFrequency = 16000.0;
    if (maxClassifierFrequency > m_parameters.sampleRate/2) {
        maxClassifierFrequency = m_parameters.sampleRate/2;
    }
    int classificationBins =
        int(floor(m_guideConfiguration.classificationFftSize *
                  maxClassifierFrequency / m_parameters.sampleRate));
    
    BinSegmenter::Parameters segmenterParameters
        (m_guideConfiguration.classificationFftSize,
         classificationBins, m_parameters.sampleRate, 18);
    
    BinClassifier::Parameters classifierParameters
        (classificationBins, 9, 1, 10, 2.0, 2.0);

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

    if (isRealTime()) {
        resamplerParameters.dynamism = Resampler::RatioOftenChanging;
        resamplerParameters.ratioChange = Resampler::SmoothRatioChange;
    } else {
        // ratio can't be changed in offline mode
        resamplerParameters.dynamism = Resampler::RatioMostlyFixed;
        resamplerParameters.ratioChange = Resampler::SuddenRatioChange;
    }
        
    resamplerParameters.initialSampleRate = m_parameters.sampleRate;
    resamplerParameters.maxBufferSize = m_guideConfiguration.longestFftSize; //!!!???
    m_resampler = std::unique_ptr<Resampler>
        (new Resampler(resamplerParameters, m_parameters.channels));

    calculateHop();

    m_prevInhop = m_inhop;
    m_prevOuthop = int(round(m_inhop * getEffectiveRatio()));

    if (!m_inhop.is_lock_free()) {
        m_parameters.logger("WARNING: std::atomic<int> is not lock-free");
    }
    if (!m_timeRatio.is_lock_free()) {
        m_parameters.logger("WARNING: std::atomic<double> is not lock-free");
    }

    // Pad to half of the longest frame. As with R2, in real-time mode
    // we don't do this -- it's better to start with a swoosh than
    // introduce more latency, and we don't want gaps when the ratio
    // changes.
    
    if (!isRealTime()) {
        m_parameters.logger("Offline mode: pre-padding");
        int pad = m_guideConfiguration.longestFftSize / 2;
        for (int c = 0; c < m_parameters.channels; ++c) {
            m_channelData[c]->inbuf->zero(pad);
        }
        // By the time we skip this later we will have resampled
        m_startSkip = int(round(pad / m_pitchScale));
    } else {
        m_parameters.logger("RT mode: no internal pre-pad");
    }
}

WindowType
R3StretcherImpl::ScaleData::analysisWindowShape(int fftSize)
{
    if (fftSize > 2048) return HannWindow;
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
    if (fftSize > 2048) return HannWindow;
    else return NiemitaloReverseWindow;
}

int
R3StretcherImpl::ScaleData::synthesisWindowLength(int fftSize)
{
    if (fftSize > 2048) return fftSize/2;
    else return fftSize;
}

void
R3StretcherImpl::setTimeRatio(double ratio)
{
    if (!isRealTime()) {
        if (m_mode == ProcessMode::Studying ||
            m_mode == ProcessMode::Processing) {
            m_parameters.logger("R3StretcherImpl::setTimeRatio: Cannot set time ratio while studying or processing in non-RT mode");
            return;
        }
    }

    if (ratio == m_timeRatio) return;
    m_timeRatio = ratio;
    calculateHop();
}

void
R3StretcherImpl::setPitchScale(double scale)
{
    if (!isRealTime()) {
        if (m_mode == ProcessMode::Studying ||
            m_mode == ProcessMode::Processing) {
            m_parameters.logger("R3StretcherImpl::setTimeRatio: Cannot set pitch scale while studying or processing in non-RT mode");
            return;
        }
    }

    if (scale == m_pitchScale) return;
    m_pitchScale = scale;
    calculateHop();
}

void
R3StretcherImpl::setFormantScale(double scale)
{
    if (!isRealTime()) {
        if (m_mode == ProcessMode::Studying ||
            m_mode == ProcessMode::Processing) {
            m_parameters.logger("R3StretcherImpl::setTimeRatio: Cannot set formant scale while studying or processing in non-RT mode");
            return;
        }
    }

    m_formantScale = scale;
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
R3StretcherImpl::setKeyFrameMap(const std::map<size_t, size_t> &mapping)
{
    if (isRealTime()) {
        m_parameters.logger("R3StretcherImpl::setKeyFrameMap: Cannot specify key frame map in RT mode");
        return;
    }
    if (m_mode == ProcessMode::Processing || m_mode == ProcessMode::Finished) {
        m_parameters.logger("R3StretcherImpl::setKeyFrameMap: Cannot specify key frame map after process() has begun");
        return;
    }

    m_keyFrameMap = mapping;
}

void
R3StretcherImpl::calculateHop()
{
    double ratio = getEffectiveRatio();

    // In R2 we generally targeted a certain inhop, and calculated
    // outhop from that. Here we are the other way around. We aim for
    // outhop = 256 at ratios around 1, reducing down to 128 for
    // ratios far below 1 and up to 512 for ratios far above. As soon
    // as outhop exceeds 256 we have to drop the 1024-bin FFT, as the
    // overlap will be inadequate for it. That's among the jobs of the
    // Guide class. (We can't go above 512 without changing the window
    // shape or dropping the 2048-bin FFT, and we can't do either of
    // those dynamically.)
    
    double proposedOuthop = pow(2.0, 8.0 + 2.0 * log10(ratio));
    if (proposedOuthop > 512.0) proposedOuthop = 512.0;
    if (proposedOuthop < 128.0) proposedOuthop = 128.0;

    std::cout << "calculateHop: for ratio " << ratio << " proposedOuthop = "
              << proposedOuthop << std::endl;
    
    double inhop = proposedOuthop / ratio;
    if (inhop < 1.0) {
        m_parameters.logger("WARNING: Extreme ratio yields ideal inhop < 1, results may be suspect");
        inhop = 1.0;
    }
    if (inhop > 768.0) {
        m_parameters.logger("WARNING: Extreme ratio yields ideal inhop > 768, results may be suspect");
        inhop = 768.0;
    }

    m_inhop = int(round(inhop));

    std::cout << "R3StretcherImpl::calculateHop: inhop = " << m_inhop << ", proposed outhop = " << proposedOuthop << ", mean outhop = " << m_inhop * ratio << std::endl;
}

void
R3StretcherImpl::updateRatioFromMap()
{
    if (m_keyFrameMap.empty()) return;
//!!!    auto itr = m_keyFrameMap.upper_bound(m_processInputDuration);
      
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

double
R3StretcherImpl::getFormantScale() const
{
    return m_formantScale;
}

size_t
R3StretcherImpl::getLatency() const
{
    if (!isRealTime()) {
        return 0;
    } else {
        double factor = m_pitchScale * 0.5;
        return size_t(ceil(m_guideConfiguration.longestFftSize * factor));
    }
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

    m_studyInputDuration = 0;
    m_totalTargetDuration = 0;
    m_processInputDuration = 0;
    m_totalOutputDuration = 0;
    m_keyFrameMap.clear();

    m_mode = ProcessMode::JustCreated;
}

void
R3StretcherImpl::study(const float *const *, size_t samples, bool)
{
    if (isRealTime()) {
        m_parameters.logger("R3StretcherImpl::study: Not meaningful in realtime mode");
        return;
    }

    if (m_mode == ProcessMode::Processing || m_mode == ProcessMode::Finished) {
        m_parameters.logger("R3StretcherImpl::study: Cannot study after processing");
        return;
    }
    
    if (m_mode == ProcessMode::JustCreated) {
        m_studyInputDuration = 0;
    }

    m_mode = ProcessMode::Studying;
    m_studyInputDuration += samples;
}

size_t
R3StretcherImpl::getSamplesRequired() const
{
    if (available() != 0) return 0;
    int longest = m_guideConfiguration.longestFftSize;
    int rs = m_channelData[0]->inbuf->getReadSpace();
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

    if (m_mode == ProcessMode::Finished) {
        m_parameters.logger("R3StretcherImpl::process: Cannot process again after final chunk");
        return;
    }

    if (!isRealTime() && !m_keyFrameMap.empty()) {
        if (m_mode == ProcessMode::Studying) {
            m_totalTargetDuration =
                round(m_studyInputDuration * getEffectiveRatio());
        }
    }

    if (final) {
        // We don't distinguish between Finished and "draining, but
        // haven't yet delivered all the samples" because the
        // distinction is meaningless internally - it only affects
        // whether available() finds any samples in the buffer
        m_mode = ProcessMode::Finished;
    } else {
        m_mode = ProcessMode::Processing;
    }
    
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

    m_processInputDuration += samples;
    
    consume();
}

//!!! __attribute__((annotate("realtime")))
int
R3StretcherImpl::available() const
{
    int av = int(m_channelData[0]->outbuf->getReadSpace());
    if (av == 0 && m_mode == ProcessMode::Finished) {
        return -1;
    } else {
        return av;
    }
}

//!!! __attribute__((annotate("realtime")))
size_t
R3StretcherImpl::retrieve(float *const *output, size_t samples) const
{
    int got = samples;
    
    for (int c = 0; c < m_parameters.channels; ++c) {
        int gotHere = m_channelData[c]->outbuf->read(output[c], got);
        if (gotHere < got) {
            if (c > 0) {
                m_parameters.logger("R3StretcherImpl::retrieve: WARNING: channel imbalance detected");
            }
            got = std::min(got, std::max(gotHere, 0));
        }
    }

    return got;
}

void
R3StretcherImpl::consume()
{
    int longest = m_guideConfiguration.longestFftSize;
    int channels = m_parameters.channels;

    //!!! todo: wire debug level & logger throughout
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
                                               longest,
                                               true);
    
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
            if (m_mode == ProcessMode::Finished) {
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

        // Phase update. This is synchronised across all channels
        
        for (auto &it : m_channelData[0]->scales) {
            int fftSize = it.first;
            for (int c = 0; c < channels; ++c) {
                auto &cd = m_channelData.at(c);
                auto &scale = cd->scales.at(fftSize);
                m_channelAssembly.mag[c] = scale->mag.data();
                m_channelAssembly.phase[c] = scale->phase.data();
                m_channelAssembly.prevMag[c] = scale->prevMag.data();
                m_channelAssembly.guidance[c] = &cd->guidance;
                m_channelAssembly.outPhase[c] = scale->advancedPhase.data();
            }
            m_scaleData.at(fftSize)->guided.advance
                (m_channelAssembly.outPhase.data(),
                 m_channelAssembly.mag.data(),
                 m_channelAssembly.phase.data(),
                 m_channelAssembly.prevMag.data(),
                 m_guideConfiguration,
                 m_channelAssembly.guidance.data(),
                 m_prevInhop,
                 m_prevOuthop);
        }

        for (int c = 0; c < channels; ++c) {
            adjustPreKick(c);
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
                 m_mode == ProcessMode::Finished && readSpace < longest);
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
                // This should happen only when draining (Finished)
                cd->inbuf->skip(readSpace);
            } else {
                cd->inbuf->skip(inhop);
            }
        }

        if (m_startSkip > 0) {
            int toSkip = std::min
                (m_startSkip, m_channelData.at(0)->outbuf->getReadSpace());
            for (int c = 0; c < channels; ++c) {
                m_channelData.at(c)->outbuf->skip(toSkip);
            }
            m_startSkip -= toSkip;
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

    if (haveValidReadahead) {
        v_copy(classifyScale->mag.data(),
               readahead.mag.data(),
               classifyScale->bufSize);
        v_copy(classifyScale->phase.data(),
               readahead.phase.data(),
               classifyScale->bufSize);
    }

    v_fftshift(readahead.timeDomain.data(), classify);
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
        if (fftSize == classify && haveValidReadahead) {
            continue;
        }
        
        auto &scale = it.second;

        v_fftshift(scale->timeDomain.data(), fftSize);

        m_scaleData.at(fftSize)->fft.forward(scale->timeDomain.data(),
                                             scale->real.data(),
                                             scale->imag.data());

        // For the classify scale we always want the full range, as
        // all the magnitudes (though not phases) are potentially
        // relevant to classification and formant analysis. But this
        // case here only happens if we don't haveValidReadahead - the
        // normal case is above and just copies from the previous
        // readahead.
        if (fftSize == classify) {
            //!!! and because not all the phases are relevant, there
            //!!! is room for an optimisation here, though this is
            //!!! used only when ratio changes
            v_cartesian_to_polar(scale->mag.data(),
                                 scale->phase.data(),
                                 scale->real.data(),
                                 scale->imag.data(),
                                 fftSize/2 + 1);
            v_scale(scale->mag.data(),
                    1.0 / double(fftSize),
                    scale->mag.size());
            continue;
        }
        
        //!!! should this be a map?
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

    if (m_parameters.options & RubberBandStretcher::OptionFormantPreserved) {
        analyseFormant(c);
        adjustFormant(c);
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
/*
    if (c == 0) {
        double pb = cd->nextSegmentation.percussiveBelow;
        double pa = cd->nextSegmentation.percussiveAbove;
        double ra = cd->nextSegmentation.residualAbove;
        int pbb = binForFrequency(pb, classify, m_parameters.sampleRate);
        int pab = binForFrequency(pa, classify, m_parameters.sampleRate);
        int rab = binForFrequency(ra, classify, m_parameters.sampleRate);
        std::cout << "pb = " << pb << ", pbb = " << pbb << std::endl;
        std::cout << "pa = " << pa << ", pab = " << pab << std::endl;
        std::cout << "ra = " << ra << ", rab = " << rab << std::endl;
        std::cout << "s:";
        for (int i = 0; i <= classify/2; ++i) {
            if (i > 0) std::cout << ",";
            if (i < pbb || (i >= pab && i <= rab)) {
                std::cout << "1";
            } else {
                std::cout << "0";
            }
        }
        std::cout << std::endl;
    }
*/
    bool specialCaseUnity = true;
        
    m_guide.updateGuidance(getEffectiveRatio(),
                           prevOuthop,
                           classifyScale->mag.data(),
                           classifyScale->prevMag.data(),
                           cd->readahead.mag.data(),
                           cd->segmentation,
                           cd->prevSegmentation,
                           cd->nextSegmentation,
                           specialCaseUnity,
                           cd->guidance);
/*
    if (c == 0) {
        if (cd->guidance.kick.present) {
            std::cout << "k:2" << std::endl;
        } else if (cd->guidance.preKick.present) {
            std::cout << "k:1" << std::endl;
        } else {
            std::cout << "k:0" << std::endl;
        }
    }
*/
}

void
R3StretcherImpl::analyseFormant(int c)
{
    auto &cd = m_channelData.at(c);
    auto &f = *cd->formant;

    int fftSize = f.fftSize;
    int binCount = fftSize/2 + 1;
    
    auto &scale = cd->scales.at(fftSize);
    auto &scaleData = m_scaleData.at(fftSize);

    scaleData->fft.inverseCepstral(scale->mag.data(), f.cepstra.data());
    
    int cutoff = int(floor(m_parameters.sampleRate / 650.0));
    if (cutoff < 1) cutoff = 1;

    f.cepstra[0] /= 2.0;
    f.cepstra[cutoff-1] /= 2.0;
    for (int i = cutoff; i < fftSize; ++i) {
        f.cepstra[i] = 0.0;
    }
    v_scale(f.cepstra.data(), 1.0 / double(fftSize), cutoff);

    scaleData->fft.forward(f.cepstra.data(), f.envelope.data(), f.spare.data());

    v_exp(f.envelope.data(), binCount);
    v_square(f.envelope.data(), binCount);

    for (int i = 0; i < binCount; ++i) {
        if (f.envelope[i] > 1.0e10) f.envelope[i] = 1.0e10;
    }
}

void
R3StretcherImpl::adjustFormant(int c)
{
    auto &cd = m_channelData.at(c);
        
    for (auto &it : cd->scales) {
        
        int fftSize = it.first;
        auto &scale = it.second;

        int highBin = int(floor(fftSize * 10000.0 / m_parameters.sampleRate));
        double targetFactor = double(cd->formant->fftSize) / double(fftSize);
        double formantScale = m_formantScale;
        if (formantScale == 0.0) formantScale = 1.0 / m_pitchScale;
        double sourceFactor = targetFactor / formantScale;
        double maxRatio = 60.0;
        double minRatio = 1.0 / maxRatio;

        for (const auto &b : m_guideConfiguration.fftBandLimits) {
            if (b.fftSize != fftSize) continue;
            for (int i = b.b0min; i < b.b1max && i < highBin; ++i) {
                double source = cd->formant->envelopeAt(i * sourceFactor);
                double target = cd->formant->envelopeAt(i * targetFactor);
                if (target > 0.0) {
                    double ratio = source / target;
                    if (ratio < minRatio) ratio = minRatio;
                    if (ratio > maxRatio) ratio = maxRatio;
                    scale->mag[i] *= ratio;
                }
            }
        }
    }
}

void
R3StretcherImpl::adjustPreKick(int c)
{
    auto &cd = m_channelData.at(c);
    auto fftSize = cd->guidance.fftBands[0].fftSize;
    if (cd->guidance.preKick.present) {
        auto &scale = cd->scales.at(fftSize);
        int from = binForFrequency(cd->guidance.preKick.f0,
                                   fftSize, m_parameters.sampleRate);
        int to = binForFrequency(cd->guidance.preKick.f1,
                                 fftSize, m_parameters.sampleRate);
        for (int i = from; i <= to; ++i) {
            double diff = scale->mag[i] - scale->prevMag[i];
            if (diff > 0.0) {
                scale->pendingKick[i] = diff;
                scale->mag[i] -= diff;
            }
        }
    } else if (cd->guidance.kick.present) {
        auto &scale = cd->scales.at(fftSize);
        int from = binForFrequency(cd->guidance.preKick.f0,
                                   fftSize, m_parameters.sampleRate);
        int to = binForFrequency(cd->guidance.preKick.f1,
                                 fftSize, m_parameters.sampleRate);
        for (int i = from; i <= to; ++i) {
            scale->mag[i] += scale->pendingKick[i];
            scale->pendingKick[i] = 0.0;
        }
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

        double winscale = double(outhop) / scaleData->windowScaleFactor;

        // The frequency filter is applied naively in the frequency
        // domain. Aliasing is reduced by the shorter resynthesis
        // window
        
        int lowBin = binForFrequency(band.f0, fftSize, m_parameters.sampleRate);
        int highBin = binForFrequency(band.f1, fftSize, m_parameters.sampleRate);
        if (highBin % 2 == 0 && highBin > 0) --highBin;

        for (int i = 0; i < lowBin; ++i) {
            scale->mag[i] = 0.0;
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

