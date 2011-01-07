/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2011 Chris Cannam.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "AudioCurveCalculator.h"

#include <iostream>

namespace RubberBand
{

static const int MaxPerceivedFreq = 16000;

AudioCurveCalculator::AudioCurveCalculator(Parameters parameters) :
    m_sampleRate(parameters.sampleRate),
    m_fftSize(parameters.fftSize)
{
    recalculateLastPerceivedBin();
}

AudioCurveCalculator::~AudioCurveCalculator()
{
}

void
AudioCurveCalculator::setSampleRate(int newRate)
{
    m_sampleRate = newRate;
    recalculateLastPerceivedBin();
}

void
AudioCurveCalculator::setFftSize(int newSize)
{
    m_fftSize = newSize;
    recalculateLastPerceivedBin();
}

void
AudioCurveCalculator::recalculateLastPerceivedBin()
{
    if (m_sampleRate == 0) {
        m_lastPerceivedBin = 0;
        return;
    }
    m_lastPerceivedBin = ((MaxPerceivedFreq * m_fftSize) / m_sampleRate);
    if (m_lastPerceivedBin > m_fftSize/2) {
        m_lastPerceivedBin = m_fftSize/2;
    }
}


}
