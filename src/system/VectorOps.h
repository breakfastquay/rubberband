/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2010 Chris Cannam.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#ifndef _RUBBERBAND_VECTOR_OPS_H_
#define _RUBBERBAND_VECTOR_OPS_H_



#include <cstring>
#include "sysutils.h"

namespace RubberBand {

// Note that all functions with a "target" vector have their arguments
// in the same order as memcpy and friends, i.e. target vector first.
// This is the reverse order from the IPP functions.

// The ideal here is to write the basic loops in such a way as to be
// auto-vectorizable by a sensible compiler (definitely gcc-4.3 on
// Linux, ideally also gcc-4.0 on OS/X).

template<typename T>
inline void v_zero(T *const R__ ptr,
                   const int count)
{
    const T value = T(0);
    for (int i = 0; i < count; ++i) {
        ptr[i] = value;
    }
}


template<typename T>
inline void v_zero_channels(T *const R__ *const R__ ptr,
                            const int channels,
                            const int count)
{
    for (int c = 0; c < channels; ++c) {
        v_zero(ptr[c], count);
    }
}

template<typename T>
inline void v_set(T *const R__ ptr,
                  const T value,
                  const int count)
{
    for (int i = 0; i < count; ++i) {
        ptr[i] = value;
    }
}

template<typename T>
inline void v_copy(T *const R__ dst,
                   const T *const R__ src,
                   const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = src[i];
    }
}


template<typename T>
inline void v_copy_channels(T *const R__ *const R__ dst,
                            const T *const R__ *const R__ src,
                            const int channels,
                            const int count)
{
    for (int c = 0; c < channels; ++c) {
        v_copy(dst[c], src[c], count);
    }
}

template<typename T>
inline void v_move(T *const R__ dst,
                   const T *const R__ src,
                   const int count)
{
    memmove(dst, src, count * sizeof(T));
}


template<typename T, typename U>
inline void v_convert(U *const R__ dst,
                      const T *const R__ src,
                      const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = U(src[i]);
    }
}

template<>
inline void v_convert(float *const R__ dst,
                      const float *const R__ src,
                      const int count)
{
    v_copy(dst, src, count);
}
template<>
inline void v_convert(double *const R__ dst,
                      const double *const R__ src,
                      const int count)
{
    v_copy(dst, src, count);
}


template<typename T, typename U>
inline void v_convert_channels(U *const R__ *const R__ dst,
                               const T *const R__ *const R__ src,
                               const int channels,
                               const int count)
{
    for (int c = 0; c < channels; ++c) {
        v_convert(dst[c], src[c], count);
    }
}

template<typename T>
inline void v_add(T *const R__ dst,
                  const T *const R__ src,
                  const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] += src[i];
    }
}


template<typename T>
inline void v_add_channels(T *const R__ *const R__ dst,
                           const T *const R__ *const R__ src,
                           const int channels, const int count)
{
    for (int c = 0; c < channels; ++c) {
        v_add(dst[c], src[c], count);
    }
}

template<typename T, typename G>
inline void v_add_with_gain(T *const R__ dst,
                            const T *const R__ src,
                            const int count,
                            const G gain)
{
    for (int i = 0; i < count; ++i) {
        dst[i] += src[i] * gain;
    }
}

template<typename T, typename G>
inline void v_add_channels_with_gain(T *const R__ *const R__ dst,
                                     const T *const R__ *const R__ src,
                                     const int channels,
                                     const int count,
                                     const G gain)
{
    for (int c = 0; c < channels; ++c) {
        v_add_with_gain(dst[c], src[c], count, gain);
    }
}

template<typename T>
inline void v_subtract(T *const R__ dst,
                       const T *const R__ src,
                       const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] -= src[i];
    }
}


template<typename T, typename G>
inline void v_scale(T *const R__ dst,
                    const G gain,
                    const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] *= gain;
    }
}


template<typename T>
inline void v_multiply(T *const R__ dst,
                       const T *const R__ src,
                       const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] *= src[i];
    }
}


template<typename T>
inline void v_multiply(T *const R__ dst,
                       const T *const R__ src1,
                       const T *const R__ src2,
                       const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = src1[i] * src2[i];
    }
}

template<typename T>
inline void v_divide(T *const R__ dst,
                     const T *const R__ src,
                     const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] /= src[i];
    }
}



template<typename T>
inline void v_multiply_and_add(T *const R__ dst,
                               const T *const R__ src1,
                               const T *const R__ src2,
                               const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] += src1[i] * src2[i];
    }
}


template<typename T>
inline void v_log(T *const R__ dst,
                  const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = log(dst[i]);
    }
}


template<typename T>
inline void v_exp(T *const R__ dst,
                  const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = exp(dst[i]);
    }
}


template<typename T>
inline void v_sqrt(T *const R__ dst,
                   const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = sqrt(dst[i]);
    }
}


template<typename T>
inline void v_square(T *const R__ dst,
                   const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = dst[i] * dst[i];
    }
}


template<typename T>
inline void v_abs(T *const R__ dst,
                  const int count)
{
    for (int i = 0; i < count; ++i) {
        dst[i] = fabs(dst[i]);
    }
}


template<typename T>
inline void v_interleave(T *const R__ dst,
                         const T *const R__ *const R__ src,
                         const int channels, 
                         const int count)
{
    int idx = 0;
    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < channels; ++j) {
            dst[idx++] = src[j][i];
        }
    }
}


template<typename T>
inline void v_deinterleave(T *const R__ *const R__ dst,
                           const T *const R__ src,
                           const int channels, 
                           const int count)
{
    int idx = 0;
    for (int i = 0; i < count; ++i) {
        for (int j = 0; j < channels; ++j) {
            dst[j][i] = src[idx++];
        }
    }
}


template<typename T>
inline void v_fftshift(T *const R__ ptr,
                       const int count)
{
    const int hs = count/2;
    for (int i = 0; i < hs; ++i) {
        T t = ptr[i];
        ptr[i] = ptr[i + hs];
        ptr[i + hs] = t;
    }
}

}

#endif
