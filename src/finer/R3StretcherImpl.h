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
        m_guideConfiguration(m_guide.getConfiguration())
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
                m_ffts[fftSize] = std::make_shared<FFT>(fftSize);
                m_channelData[c]->scales[fftSize] =
                    std::make_shared<ChannelScaleData>(fftSize);
            }
        }
    }
    
    ~R3StretcherImpl() { }

    void reset();
    
    void setTimeRatio(double ratio);
    void setPitchScale(double scale);

    double getTimeRatio() const;
    double getPitchScale() const;

protected:
    struct ChannelScaleData {
        int fftSize;
        int bufSize; // size of every freq-domain array here: fftSize/2 + 1
        FixedVector<float> mag;
        FixedVector<float> phase;
        FixedVector<int> nearestPeaks;
        FixedVector<int> nearestTroughs;
        FixedVector<float> prevOutMag;
        FixedVector<double> prevOutPhase;
        FixedVector<int> prevNearestPeaks;
        FixedVector<float> timeDomainFrame;
        Window<float> analysisWindow;
        Window<float> synthesisWindow;

        ChannelScaleData(int _fftSize) :
            fftSize(_fftSize), bufSize(fftSize/2 + 1),
            mag(bufSize, 0.f), phase(bufSize, 0.f),
            nearestPeaks(bufSize, 0), nearestTroughs(bufSize, 0),
            prevOutMag(bufSize, 0.f), prevOutPhase(bufSize, 0.f),
            prevNearestPeaks(bufSize, 0), timeDomainFrame(fftSize, 0.f),
            analysisWindow(HannWindow, fftSize),
            synthesisWindow(HannWindow, fftSize/2)
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
        RingBuffer<float> inbuf;
        RingBuffer<float> outbuf;
        ChannelData(BinSegmenter::Parameters segmenterParameters,
                    BinClassifier::Parameters classifierParameters,
                    int ringBufferSize) :
            scales(),
            segmenter(new BinSegmenter(segmenterParameters,
                                       classifierParameters)),
            segmentation(), prevSegmentation(), nextSegmentation(),
            inbuf(ringBufferSize), outbuf(ringBufferSize) { }
    };

    Parameters m_parameters;

    double m_timeRatio;
    double m_pitchScale;

    std::vector<std::shared_ptr<ChannelData>> m_channelData;
    std::map<int, std::shared_ptr<FFT>> m_ffts;
    Guide m_guide;
    Guide::Configuration m_guideConfiguration;

    static void logCerr(const std::string &message) {
        std::cerr << "RubberBandStretcher: " << message << std::endl;
    }
};

}

#endif
