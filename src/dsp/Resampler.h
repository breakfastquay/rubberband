/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*- vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band Library
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2020 Particular Programs Ltd.

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

#ifndef _RUBBERBAND_RESAMPLER_H_
#define _RUBBERBAND_RESAMPLER_H_

#include "system/sysutils.h"

namespace RubberBand {

class ResamplerImpl;

class Resampler
{
public:
    enum Quality { Best, FastestTolerable, Fastest };
    enum Exception { ImplementationError };

    /**
     * Construct a resampler with the given quality level and channel
     * count.  maxBufferSize gives a bound on the maximum incount size
     * that may be passed to the resample function before the
     * resampler needs to reallocate its internal buffers.
     */
    Resampler(Quality quality, int channels, int maxBufferSize = 0,
              int debugLevel = 0);
    ~Resampler();

    /**
     * Resample the given multi-channel buffers, where incount is the
     * number of frames in the input buffers.  Returns the number of
     * frames written to the output buffers.
     */
    int resample(const float *const R__ *const R__ in,
                 float *const R__ *const R__ out,
                 int incount,
                 float ratio,
                 bool final = false);

    /**
     * Resample the given interleaved buffer, where incount is the
     * number of frames in the input buffer (i.e. it has incount *
     * getChannelCount() samples).  Returns the number of frames
     * written to the output buffer.
     */
    int resampleInterleaved(const float *const R__ in,
                            float *const R__ out,
                            int incount,
                            float ratio,
                            bool final = false);

    int getChannelCount() const;

    void reset();

protected:
    ResamplerImpl *d;
    int m_method;
};

}

#endif
