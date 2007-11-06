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

#include "ConstantAudioCurve.h"

namespace RubberBand
{

ConstantAudioCurve::ConstantAudioCurve(size_t sampleRate, size_t blockSize) :
    AudioCurve(sampleRate, blockSize)
{
}

ConstantAudioCurve::~ConstantAudioCurve()
{
}

void
ConstantAudioCurve::reset()
{
}

void
ConstantAudioCurve::setBlockSize(size_t newSize)
{
    m_blockSize = newSize;
}

float
ConstantAudioCurve::process(float *, size_t)
{
    return 1.f;
}

}

