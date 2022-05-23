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
    
    R3StretcherImpl(Parameters parameters) :
        m_parameters(parameters),
        m_guide(Guide::Parameters(m_parameters.sampleRate)),
        m_guideConfiguration(m_guide.getConfiguration()),
        m_channelAssembly(m_parameters.channels),
        m_troughPicker(m_guideConfiguration.classificationFftSize / 2 + 1)
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
                    std::make_shared<ChannelScaleData>(fftSize);
            }
        }
        
        for (auto band: m_guideConfiguration.fftBandLimits) {
            int fftSize = band.fftSize;
            GuidedPhaseAdvance::Parameters guidedParameters
                (fftSize, m_parameters.sampleRate, m_parameters.channels,
                 m_parameters.logger);
            m_scaleData[fftSize] = std::make_shared<ScaleData>(guidedParameters);
        }
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
    struct ChannelScaleData {
        int fftSize;
        int bufSize; // size of every freq-domain array here: fftSize/2 + 1
        //!!! review later which of these we are actually using!
        FixedVector<float> timeDomainFrame;
        FixedVector<float> mag;
        FixedVector<float> phase;
        FixedVector<double> outPhase; //!!! "advanced"?
        FixedVector<int> nextTroughs; //!!! not used in every scale
        FixedVector<float> prevMag; //!!! not used in every scale
        FixedVector<double> prevOutPhase;
        FixedVector<float> accumulator;

        ChannelScaleData(int _fftSize) :
            fftSize(_fftSize),
            bufSize(fftSize/2 + 1),
            timeDomainFrame(fftSize, 0.f),
            mag(bufSize, 0.f),
            phase(bufSize, 0.f),
            outPhase(bufSize, 0.f),
            nextTroughs(bufSize, 0),
            prevMag(bufSize, 0.f),
            prevOutPhase(bufSize, 0.f),
            accumulator(fftSize, 0.f)
        { }

    private:
        ChannelScaleData(const ChannelScaleData &) =delete;
        ChannelScaleData &operator=(const ChannelScaleData &) =delete;
    };

    struct ChannelData {
        std::map<int, std::shared_ptr<ChannelScaleData>> scales;
        std::unique_ptr<BinSegmenter> segmenter;
        BinSegmenter::Segmentation segmentation;
        BinSegmenter::Segmentation prevSegmentation;
        BinSegmenter::Segmentation nextSegmentation;
        Guide::Guidance guidance;
        FixedVector<float> mixdown;
        std::unique_ptr<RingBuffer<float>> inbuf;
        std::unique_ptr<RingBuffer<float>> outbuf;
        ChannelData(BinSegmenter::Parameters segmenterParameters,
                    BinClassifier::Parameters classifierParameters,
                    int ringBufferSize) :
            scales(),
            segmenter(new BinSegmenter(segmenterParameters,
                                       classifierParameters)),
            segmentation(), prevSegmentation(), nextSegmentation(),
            mixdown(ringBufferSize, 0.f), //!!! could be much shorter (bound is the max outhop)
            inbuf(new RingBuffer<float>(ringBufferSize)),
            outbuf(new RingBuffer<float>(ringBufferSize)) { }
    };

    struct ChannelAssembly {
        // Vectors of bare pointers, used to package container data
        // from different channels into arguments for PhaseAdvance
        FixedVector<float *> mag;
        FixedVector<float *> phase;
        FixedVector<Guide::Guidance *> guidance;
        FixedVector<double *> outPhase;
        ChannelAssembly(int channels) :
            mag(channels, nullptr), phase(channels, nullptr),
            guidance(channels, nullptr), outPhase(channels, nullptr) { }
    };

    struct ScaleData {
        FFT fft;
        Window<float> analysisWindow;
        Window<float> synthesisWindow;
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
    Peak<float, std::less<float>> m_troughPicker;

    void consume();
    
    static void logCerr(const std::string &message) {
        std::cerr << "RubberBandStretcher: " << message << std::endl;
    }
};

}

#endif
