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

#ifndef _RUBBERBAND_ALLOCATORS_H_
#define _RUBBERBAND_ALLOCATORS_H_

#include "VectorOps.h"

#include <new> // for std::bad_alloc
#include <stdlib.h>

#ifndef HAVE_POSIX_MEMALIGN
#ifndef _WIN32
#ifndef __APPLE__
#ifndef LACK_POSIX_MEMALIGN
#define HAVE_POSIX_MEMALIGN
#endif
#endif
#endif
#endif

#ifdef HAVE_POSIX_MEMALIGN
#include <sys/mman.h>
#endif

namespace RubberBand {

template <typename T>
T *allocate(size_t count)
{
    void *ptr = 0;
#ifdef HAVE_POSIX_MEMALIGN
    if (posix_memalign(&ptr, 16, count * sizeof(T))) {
        ptr = malloc(count * sizeof(T));
    }
#else
#ifdef _WIN32
    ptr = _aligned_malloc(count * sizeof(T), 16);
#else
    // Note that malloc always aligns to 16 byte boundaries on OS/X,
    // so we don't need posix_memalign there (which is fortunate,
    // since it doesn't exist)
    ptr = malloc(count * sizeof(T));
#endif
#endif
    if (!ptr) throw(std::bad_alloc());
    return (T *)ptr;
}

	
template <typename T>
T *allocate_and_zero(size_t count)
{
    T *ptr = allocate<T>(count);
    v_zero(ptr, count);
    return ptr;
}

template <typename T>
void deallocate(T *ptr)
{
#ifdef _WIN32
    if (ptr) _aligned_free((void *)ptr);
#else
    if (ptr) free((void *)ptr);
#endif
}

	
template <typename T>
T *reallocate(T *ptr, size_t oldcount, size_t count)
{
    T *newptr = 0;
    try {
        newptr = allocate<T>(count);
    } catch (std::bad_alloc) {
        if (ptr) deallocate<T>(ptr);
        throw;
    }
    if (oldcount && ptr) {
        v_copy(newptr, ptr, oldcount < count ? oldcount : count);
    }
    if (ptr) deallocate<T>(ptr);
    return newptr;
}

template <typename T>
T **allocate_channels(size_t channels, size_t count)
{
    T **ptr = allocate<T *>(channels);
    for (size_t c = 0; c < channels; ++c) {
        ptr[c] = allocate<T>(count);
    }
    return ptr;
}
	
template <typename T>
T **allocate_and_zero_channels(size_t channels, size_t count)
{
    T **ptr = allocate<T *>(channels);
    for (size_t c = 0; c < channels; ++c) {
        ptr[c] = allocate_and_zero<T>(count);
    }
    return ptr;
}

template <typename T>
void deallocate_channels(T **ptr, size_t channels)
{
    if (!ptr) return;
    for (size_t c = 0; c < channels; ++c) {
        deallocate<T>(ptr[c]);
    }
    deallocate<T *>(ptr);
}
	
template <typename T>
T **reallocate_channels(T **ptr,
                        size_t oldchannels, size_t oldcount,
                        size_t channels, size_t count)
{
    T **newptr = 0;
    try {
        newptr = allocate_channels<T>(channels, count);
    } catch (std::bad_alloc) {
        if (ptr) deallocate_channels<T>(ptr);
        throw;
    }
    if (oldcount && ptr) {
        v_copy_channels(newptr, ptr, channels, oldcount < count ? oldcount : count);
    } 
    if (ptr) deallocate_channels<T>(ptr);
    return newptr;
}

}

#endif

