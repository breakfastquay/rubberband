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

#include <map>
#include <memory>

#include "BinSegmenter.h"
#include "Guide.h"
#include "Peak.h"
#include "PhaseAdvance.h"

#include "../common/FFT.h"
#include "../common/FixedVector.h"
#include "../common/Allocators.h"
#include "../common/Window.h"

namespace RubberBand
{

class R3StretcherImpl
{
public:
    R3StretcherImpl(double sampleRate, int channels) :
        m_sampleRate(sampleRate), m_channels(channels),
        m_guide(Guide::Parameters(sampleRate)),
        m_guideConfiguration(m_guide.getConfiguration())
    { }
    ~R3StretcherImpl();

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
    };

    double m_sampleRate;
    int m_channels;

    double m_timeRatio;
    double m_pitchScale;

    std::vector<ChannelData> m_channelData;
    std::map<int, std::shared_ptr<FFT>> m_ffts;
    Guide m_guide;
    Guide::Configuration m_guideConfiguration;
};

}

#endif
