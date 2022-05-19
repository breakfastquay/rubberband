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

#include "../common/FFT.h"
#include "../common/Allocators.h"

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
    double m_sampleRate;
    int m_channels;

    double m_timeRatio;
    double m_pitchScale;

    struct ChannelScaleData {
        int fftSize;
        int bufSize; // size of every array here: fftSize/2 + 1
        float *mag;
        float *phase;
        int *nearestPeaks;
        int *nearestTroughs;
        float *prevOutMag;
        float *prevOutPhase;
        int *prevNearestPeaks;

        ChannelScaleData(int _fftSize) :
            fftSize(_fftSize), bufSize(_fftSize/2 + 1),
            mag(allocate_and_zero<float>(size_t(bufSize))),
            phase(allocate_and_zero<float>(size_t(bufSize))),
            nearestPeaks(allocate_and_zero<int>(size_t(bufSize))),
            nearestTroughs(allocate_and_zero<int>(size_t(bufSize))),
            prevOutMag(allocate_and_zero<float>(size_t(bufSize))),
            prevOutPhase(allocate_and_zero<float>(size_t(bufSize))),
            prevNearestPeaks(allocate_and_zero<int>(size_t(bufSize))) { }

        ~ChannelScaleData() {
            deallocate(mag);
            deallocate(phase);
            deallocate(nearestPeaks);
            deallocate(nearestTroughs);
            deallocate(prevOutMag);
            deallocate(prevOutPhase);
            deallocate(prevNearestPeaks);
        }

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
        
    };
    
    std::map<int, std::shared_ptr<FFT>> m_ffts;
    Guide m_guide;
    Guide::Configuration m_guideConfiguration;
};

}

#endif
