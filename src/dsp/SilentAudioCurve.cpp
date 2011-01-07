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

#include "SilentAudioCurve.h"

#include <cmath>

namespace RubberBand
{


SilentAudioCurve::SilentAudioCurve(Parameters parameters) :
    AudioCurveCalculator(parameters)
{
}

SilentAudioCurve::~SilentAudioCurve()
{
}

void
SilentAudioCurve::reset()
{
}

float
SilentAudioCurve::processFloat(const float *R__ mag, int)
{
    const int hs = m_lastPerceivedBin;
    static float threshold = powf(10.f, -6);

    for (int i = 0; i <= hs; ++i) {
        if (mag[i] > threshold) return 0.f;
    }
        
    return 1.f;
}

double
SilentAudioCurve::processDouble(const double *R__ mag, int)
{
    const int hs = m_lastPerceivedBin;
    static double threshold = pow(10.0, -6);

    for (int i = 0; i <= hs; ++i) {
        if (mag[i] > threshold) return 0.f;
    }
        
    return 1.f;
}

}

