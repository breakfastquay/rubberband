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

#ifndef _AUDIO_CURVE_CALCULATOR_H_
#define _AUDIO_CURVE_CALCULATOR_H_

#include <sys/types.h>


#include "system/sysutils.h"

namespace RubberBand 
{

class AudioCurveCalculator
{
public:
    AudioCurveCalculator(size_t sampleRate, size_t windowSize);
    virtual ~AudioCurveCalculator();

    size_t getSampleRate() const { return m_sampleRate; }
    size_t getWindowSize() const { return m_windowSize; }

    virtual void setWindowSize(size_t newSize) = 0;

    // You may not mix calls to the various process functions on a
    // given instance


    virtual float processFloat(const float *R__ mag, size_t increment) = 0;
    virtual double processDouble(const double *R__ mag, size_t increment) = 0;

    virtual void reset() = 0;

protected:
    size_t m_sampleRate;
    size_t m_windowSize;
};

}

#endif

