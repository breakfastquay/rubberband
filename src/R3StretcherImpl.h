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

#include "dsp/BinSegmenter.h"
#include "dsp/FFT.h"
#include "system/Allocators.h"

#include "dsp/Peak.h"

namespace RubberBand
{

class R3StretcherImpl
{
public:
    R3StretcherImpl(int sampleRate, int channels);
    ~R3StretcherImpl();

    void reset();
    
    void setTimeRatio(double ratio);
    void setPitchScale(double scale);

    double getTimeRatio() const;
    double getPitchScale() const;

protected:
    int m_sampleRate;
    int m_channels;

    double m_timeRatio;
    double m_pitchScale;

    struct FftBand {
        int fftSize;
        float f0;
        float f1;
        FftBand(int _s, float _f0, float _f1) :
            fftSize(_s), f0(_f0), f1(_f1) { }
    };

    struct PhaseLockBand {
        int p;
        float beta;
        float f0;
        float f1;
        PhaseLockBand(int _p, float _beta, float _f0, float _f1) :
            p(_p), beta(_beta), f0(_f0), f1(_f1) { }
    };

    struct Range {
        bool present;
        float f0;
        float f1;
        Range() : present(false), f0(0.f), f1(0.f) { }
    };
        
    struct Guidance {
        FftBand fftBands[3];
        PhaseLockBand phaseLockBands[5];
        Range kick;
        Range lowPercussive;
        Range phaseReset;
        Range highPercussive;
        Range channelLock;
    };

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
        Guidance guidance;
    };
    
    std::map<int, std::shared_ptr<FFT>> m_ffts;

};

}

#endif
