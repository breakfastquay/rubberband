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

#include "ConstantAudioCurve.h"

namespace RubberBand
{


ConstantAudioCurve::ConstantAudioCurve(Parameters parameters) :
    AudioCurveCalculator(parameters)
{
}

ConstantAudioCurve::~ConstantAudioCurve()
{
}

void
ConstantAudioCurve::reset()
{
}

float
ConstantAudioCurve::processFloat(const float *R__, int)
{
    return 1.f;
}

double
ConstantAudioCurve::processDouble(const double *R__, int)
{
    return 1.0;
}

}

