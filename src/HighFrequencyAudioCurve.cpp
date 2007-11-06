/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band
    An audio time-stretching and pitch-shifting library.
    Copyright 2007 Chris Cannam.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "HighFrequencyAudioCurve.h"

namespace RubberBand
{

HighFrequencyAudioCurve::HighFrequencyAudioCurve(size_t sampleRate, size_t blockSize) :
    AudioCurve(sampleRate, blockSize)
{
    m_prevMag = new double[m_blockSize/2 + 1];

    for (size_t i = 0; i <= m_blockSize/2; ++i) {
        m_prevMag[i] = 0.f;
    }
}

HighFrequencyAudioCurve::~HighFrequencyAudioCurve()
{
    delete[] m_prevMag;
}

void
HighFrequencyAudioCurve::reset()
{
    for (size_t i = 0; i <= m_blockSize/2; ++i) {
        m_prevMag[i] = 0;
    }
}

void
HighFrequencyAudioCurve::setBlockSize(size_t newSize)
{
    m_blockSize = newSize;
}

float
HighFrequencyAudioCurve::process(float *mag, size_t increment)
{
    float result = 0.0;

    for (size_t n = 0; n <= m_blockSize / 2; ++n) {
        result += mag[n];
    }

    return result;
}

}

