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

#ifndef _RUBBERBAND_VECTOR_OPS_COMPLEX_H_
#define _RUBBERBAND_VECTOR_OPS_COMPLEX_H_

#include "VectorOps.h"


namespace RubberBand {


template<typename T>
inline void c_phasor(T *real, T *imag, T phase)
{
    //!!! IPP contains ippsSinCos_xxx in ippvm.h -- these are
    //!!! fixed-accuracy, test and compare
#if defined __GNUC__
    if (sizeof(T) == sizeof(float)) {
        sincosf(phase, (float *)imag, (float *)real);
    } else {
        sincos(phase, (double *)imag, (double *)real);
    }
#else
    if (sizeof(T) == sizeof(float)) {
        *real = cosf(phase);
        *imag = sinf(phase);
    } else {
        *real = cos(phase);
        *imag = sin(phase);
    }
#endif
}

template<typename T>
inline void c_magphase(T *mag, T *phase, T real, T imag)
{
    *mag = sqrt(real * real + imag * imag);
    *phase = atan2(imag, real);
}

template<>
inline void c_magphase(float *mag, float *phase, float real, float imag)
{
    *mag = sqrtf(real * real + imag * imag);
    *phase = atan2f(imag, real);
}


template<typename T>
void v_polar_to_cartesian(T *const R__ real,
                          T *const R__ imag,
                          T *const R__ mag,
                          T *const R__ phase,
                          const int count)
{
    for (int i = 0; i < count; ++i) {
        c_phasor(real + i, imag + i, phase[i]);
    }
    v_multiply(real, mag, count);
    v_multiply(imag, mag, count);
}

template<typename T>
void v_polar_interleaved_to_cartesian_inplace(T *const R__ srcdst,
                                              const int count)
{
    T real, imag;
    for (int i = 0; i < count*2; i += 2) {
        c_phasor(&real, &imag, srcdst[i+1]);
        real *= srcdst[i];
        imag *= srcdst[i];
        srcdst[i] = real;
        srcdst[i+1] = imag;
    }
}

template<typename T>
void v_cartesian_to_polar(T *const R__ mag,
                          T *const R__ phase,
                          T *const R__ real,
                          T *const R__ imag,
                          const int count)
{
    for (int i = 0; i < count; ++i) {
        c_magphase(mag + i, phase + i, real[i], imag[i]);
    }
}

template<typename T>
void v_cartesian_interleaved_to_polar(T *const R__ mag,
                                      T *const R__ phase,
                                      const T *const R__ src,
                                      const int count)
{
    for (int i = 0; i < count; ++i) {
        c_magphase(mag + i, phase + i, src[i*2], src[i*2+1]);
    }
}


template<typename T>
void v_cartesian_to_polar_interleaved_inplace(T *const R__ srcdst,
                                              const int count)
{
    T mag, phase;
    for (int i = 0; i < count * 2; i += 2) {
        c_magphase(&mag, &phase, srcdst[i], srcdst[i+1]);
        srcdst[i] = mag;
        srcdst[i+1] = phase;
    }
}

}

#endif

