/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2009 Chris Cannam.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "PercussiveAudioCurve.h"

#include "system/Allocators.h"
#include "system/VectorOps.h"

#include <cmath>

namespace RubberBand
{

PercussiveAudioCurve::PercussiveAudioCurve(size_t sampleRate, size_t windowSize) :
    AudioCurveCalculator(sampleRate, windowSize)
{
    m_prevMag = allocate_and_zero<double>(m_windowSize/2 + 1);
}

PercussiveAudioCurve::~PercussiveAudioCurve()
{
    deallocate(m_prevMag);
}

void
PercussiveAudioCurve::reset()
{
    v_zero(m_prevMag, m_windowSize/2 + 1);
}

void
PercussiveAudioCurve::setWindowSize(size_t newSize)
{
    m_prevMag = reallocate(m_prevMag, m_windowSize, newSize);
    m_windowSize = newSize;
    reset();
}

float
PercussiveAudioCurve::processFloat(const float *R__ mag, size_t increment)
{
    static float threshold = powf(10.f, 0.15f); // 3dB rise in square of magnitude
    static float zeroThresh = powf(10.f, -8);

    size_t count = 0;
    size_t nonZeroCount = 0;

    const int sz = m_windowSize / 2;

    for (int n = 1; n <= sz; ++n) {
        bool above = ((mag[n] / m_prevMag[n]) >= threshold);
        if (above) ++count;
        if (mag[n] > zeroThresh) ++nonZeroCount;
    }

    v_convert(m_prevMag, mag, sz + 1);

    if (nonZeroCount == 0) return 0;
    else return float(count) / float(nonZeroCount);
}

double
PercussiveAudioCurve::processDouble(const double *R__ mag, size_t increment)
{
    static double threshold = powf(10., 0.15); // 3dB rise in square of magnitude
    static double zeroThresh = powf(10., -8);

    size_t count = 0;
    size_t nonZeroCount = 0;

    const int sz = m_windowSize / 2;

    for (int n = 1; n <= sz; ++n) {
        bool above = ((mag[n] / m_prevMag[n]) >= threshold);
        if (above) ++count;
        if (mag[n] > zeroThresh) ++nonZeroCount;
    }

    v_copy(m_prevMag, mag, sz + 1);

    if (nonZeroCount == 0) return 0;
    else return double(count) / double(nonZeroCount);
}


}

