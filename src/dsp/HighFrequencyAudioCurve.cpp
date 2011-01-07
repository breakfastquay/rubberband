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

#include "HighFrequencyAudioCurve.h"

namespace RubberBand
{


HighFrequencyAudioCurve::HighFrequencyAudioCurve(Parameters parameters) :
    AudioCurveCalculator(parameters)
{
}

HighFrequencyAudioCurve::~HighFrequencyAudioCurve()
{
}

void
HighFrequencyAudioCurve::reset()
{
}

float
HighFrequencyAudioCurve::processFloat(const float *R__ mag, int increment)
{
    float result = 0.0;

    const int sz = m_lastPerceivedBin;

    for (int n = 0; n <= sz; ++n) {
        result = result + mag[n] * n;
    }

    return result;
}

double
HighFrequencyAudioCurve::processDouble(const double *R__ mag, int increment)
{
    float result = 0.0;

    const int sz = m_lastPerceivedBin;

    for (int n = 0; n <= sz; ++n) {
        result = result + mag[n] * n;
    }

    return result;
}

}

