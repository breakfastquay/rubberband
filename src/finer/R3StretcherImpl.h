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

#ifndef RUBBERBAND_R3_STRETCHERIMPL_H
#define RUBBERBAND_R3_STRETCHERIMPL_H

#include "BinSegmenter.h"
#include "Guide.h"
#include "Peak.h"
#include "PhaseAdvance.h"

#include "../common/StretchCalculator.h"
#include "../common/Resampler.h"
#include "../common/FFT.h"
#include "../common/FixedVector.h"
#include "../common/Allocators.h"
#include "../common/Window.h"

#include <map>
#include <memory>
#include <functional>

namespace RubberBand
{

class R3StretcherImpl
{
public:
    struct Parameters {
        double sampleRate;
        int channels;
        std::function<void(const std::string &)> logger;
        Parameters(double _sampleRate, int _channels,
                   std::function<void(const std::string &)> _log = &logCerr) :
            sampleRate(_sampleRate), channels(_channels), logger(_log) { }
    };
    
    R3StretcherImpl(Parameters parameters,
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

        int ringBufferSize = m_guideConfiguration.longestFftSize * 2;

        for (int c = 0; c < m_parameters.channels; ++c) {
            m_channelData.push_back(std::make_shared<ChannelData>
                                    (segmenterParameters,
                                     classifierParameters,
                                     ringBufferSize));
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
        
        calculateHop();
    }
    
    ~R3StretcherImpl() { }

    void reset();
    
    void setTimeRatio(double ratio);
    void setPitchScale(double scale);

    double getTimeRatio() const;
    double getPitchScale() const;

    size_t getSamplesRequired() const;
    void process(const float *const *input, size_t samples, bool final);
    int available() const;
    size_t retrieve(float *const *output, size_t samples) const;

    size_t getLatency() const;
    size_t getChannelCount() const;
    
protected:
    struct ClassificationReadaheadData {
        FixedVector<double> timeDomain;
        FixedVector<double> mag;
        FixedVector<double> phase;
        ClassificationReadaheadData(int _fftSize) :
            timeDomain(_fftSize, 0.f),
            mag(_fftSize/2 + 1, 0.f),
            phase(_fftSize/2 + 1, 0.f)
        { }

    private:
        ClassificationReadaheadData(const ClassificationReadaheadData &) =delete;
        ClassificationReadaheadData &operator=(const ClassificationReadaheadData &) =delete;
    };
    
    struct ChannelScaleData {
        int fftSize;
        int bufSize; // size of every freq-domain array here: fftSize/2 + 1
        //!!! review later which of these we are actually using!
        FixedVector<double> timeDomain;
        FixedVector<double> real;
        FixedVector<double> imag;
        FixedVector<double> mag;
        FixedVector<double> phase;
        FixedVector<double> advancedPhase;
        FixedVector<int> troughs; //!!! not used in every scale
        FixedVector<double> prevMag; //!!! not used in every scale
        FixedVector<double> prevAdvancedPhase;
        FixedVector<double> accumulator;

        ChannelScaleData(int _fftSize, int _longestFftSize) :
            fftSize(_fftSize),
            bufSize(fftSize/2 + 1),
            timeDomain(fftSize, 0.f),
            real(bufSize, 0.f),
            imag(bufSize, 0.f),
            mag(bufSize, 0.f),
            phase(bufSize, 0.f),
            advancedPhase(bufSize, 0.f),
            troughs(bufSize, 0),
            prevMag(bufSize, 0.f),
            prevAdvancedPhase(bufSize, 0.f),
            accumulator(_longestFftSize, 0.f)
        { }

    private:
        ChannelScaleData(const ChannelScaleData &) =delete;
        ChannelScaleData &operator=(const ChannelScaleData &) =delete;
    };

    struct ChannelData {
        std::map<int, std::shared_ptr<ChannelScaleData>> scales;
        ClassificationReadaheadData readahead;
        std::unique_ptr<BinSegmenter> segmenter;
        BinSegmenter::Segmentation segmentation;
        BinSegmenter::Segmentation prevSegmentation;
        BinSegmenter::Segmentation nextSegmentation;
        Guide::Guidance guidance;
        FixedVector<double> mixdown;
        std::unique_ptr<RingBuffer<float>> inbuf;
        std::unique_ptr<RingBuffer<float>> outbuf;
        ChannelData(BinSegmenter::Parameters segmenterParameters,
                    BinClassifier::Parameters classifierParameters,
                    int ringBufferSize) :
            scales(),
            readahead(segmenterParameters.fftSize),
            segmenter(new BinSegmenter(segmenterParameters,
                                       classifierParameters)),
            segmentation(), prevSegmentation(), nextSegmentation(),
            mixdown(ringBufferSize, 0.f), //!!! could be shorter (bound is the max fft size I think)
            inbuf(new RingBuffer<float>(ringBufferSize)),
            outbuf(new RingBuffer<float>(ringBufferSize)) { }
    };

    struct ChannelAssembly {
        // Vectors of bare pointers, used to package container data
        // from different channels into arguments for PhaseAdvance
        FixedVector<double *> mag;
        FixedVector<double *> phase;
        FixedVector<Guide::Guidance *> guidance;
        FixedVector<double *> outPhase;
        ChannelAssembly(int channels) :
            mag(channels, nullptr), phase(channels, nullptr),
            guidance(channels, nullptr), outPhase(channels, nullptr) { }
    };

    struct ScaleData {
        FFT fft;
        Window<double> analysisWindow;
        Window<double> synthesisWindow;
        GuidedPhaseAdvance guided;
        ScaleData(GuidedPhaseAdvance::Parameters guidedParameters) :
            fft(guidedParameters.fftSize),
            analysisWindow(HannWindow, guidedParameters.fftSize),
            synthesisWindow(HannWindow, guidedParameters.fftSize/2),
            guided(guidedParameters) { }
    };
    
    Parameters m_parameters;

    double m_timeRatio;
    double m_pitchScale;
    
    std::vector<std::shared_ptr<ChannelData>> m_channelData;
    std::map<int, std::shared_ptr<ScaleData>> m_scaleData;
    Guide m_guide;
    Guide::Configuration m_guideConfiguration;
    ChannelAssembly m_channelAssembly;
    Peak<double, std::less<double>> m_troughPicker;
    std::unique_ptr<StretchCalculator> m_calculator;
    std::unique_ptr<Resampler> m_resampler;
    int m_inhop;
    int m_prevOuthop;
    bool m_draining;

    void consume();
    void calculateHop();

    double getEffectiveRatio() const {
        return m_timeRatio * m_pitchScale;
    }
    
    static void logCerr(const std::string &message) {
        std::cerr << "RubberBandStretcher: " << message << std::endl;
    }
};

}

#endif
