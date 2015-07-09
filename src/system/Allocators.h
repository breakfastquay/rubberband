/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band Library
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2015 Particular Programs Ltd.

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

#ifdef LACK_BAD_ALLOC
namespace std { struct bad_alloc { }; }
#endif

namespace RubberBand {

template <typename T>
T *allocate(size_t count)
{
    void *ptr = 0;
    // 32-byte alignment is required for at least OpenMAX
    static const int alignment = 32;
#ifdef USE_OWN_ALIGNED_MALLOC
    // Alignment must be a power of two, bigger than the pointer
    // size. Stuff the actual malloc'd pointer in just before the
    // returned value.  This is the least desirable way to do this --
    // the other options below are all better
    size_t allocd = count * sizeof(T) + alignment;
    void *buf = malloc(allocd);
    if (buf) {
        char *adj = (char *)buf;
        while ((unsigned long long)adj & (alignment-1)) --adj;
        ptr = ((char *)adj) + alignment;
        ((void **)ptr)[-1] = buf;
    }
#else /* !USE_OWN_ALIGNED_MALLOC */
#ifdef HAVE_POSIX_MEMALIGN
    if (posix_memalign(&ptr, alignment, count * sizeof(T))) {
        ptr = malloc(count * sizeof(T));
    }
#else /* !HAVE_POSIX_MEMALIGN */
#ifdef __MSVC__
    ptr = _aligned_malloc(count * sizeof(T), alignment);
#else /* !__MSVC__ */
#ifndef MALLOC_IS_ALIGNED
#warning "No aligned malloc available or defined"
#endif
    // Note that malloc always aligns to 16 byte boundaries on OS/X
    ptr = malloc(count * sizeof(T));
#endif /* !__MSVC__ */
#endif /* !HAVE_POSIX_MEMALIGN */
#endif /* !USE_OWN_ALIGNED_MALLOC */
    if (!ptr) {
#ifndef NO_EXCEPTIONS
        throw(std::bad_alloc());
#else
        abort();
#endif
    }
    return (T *)ptr;
}

#ifdef HAVE_IPP

template <>
float *allocate(size_t count);

template <>
double *allocate(size_t count);

#endif
	
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
#ifdef USE_OWN_ALIGNED_MALLOC
    if (ptr) free(((void **)ptr)[-1]);
#else /* !USE_OWN_ALIGNED_MALLOC */
#ifdef __MSVC__
    if (ptr) _aligned_free((void *)ptr);
#else /* !__MSVC__ */
    if (ptr) free((void *)ptr);
#endif /* !__MSVC__ */
#endif /* !USE_OWN_ALIGNED_MALLOC */
}

#ifdef HAVE_IPP

template <>
void deallocate(float *);

template <>
void deallocate(double *);

#endif

/// Reallocate preserving contents but leaving additional memory uninitialised	
template <typename T>
T *reallocate(T *ptr, size_t oldcount, size_t count)
{
    T *newptr = allocate<T>(count);
    if (oldcount && ptr) {
        v_copy(newptr, ptr, oldcount < count ? oldcount : count);
    }
    if (ptr) deallocate<T>(ptr);
    return newptr;
}

/// Reallocate, zeroing all contents
template <typename T>
T *reallocate_and_zero(T *ptr, size_t oldcount, size_t count)
{
    ptr = reallocate(ptr, oldcount, count);
    v_zero(ptr, count);
    return ptr;
}
	
/// Reallocate preserving contents and zeroing any additional memory	
template <typename T>
T *reallocate_and_zero_extension(T *ptr, size_t oldcount, size_t count)
{
    ptr = reallocate(ptr, oldcount, count);
    if (count > oldcount) v_zero(ptr + oldcount, count - oldcount);
    return ptr;
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
    T **newptr = allocate_channels<T>(channels, count);
    if (oldcount && ptr) {
        v_copy_channels(newptr, ptr, channels, oldcount < count ? oldcount : count);
    } 
    if (ptr) deallocate_channels<T>(ptr, channels);
    return newptr;
}
	
template <typename T>
T **reallocate_and_zero_extend_channels(T **ptr,
                                        size_t oldchannels, size_t oldcount,
                                        size_t channels, size_t count)
{
    T **newptr = allocate_and_zero_channels<T>(channels, count);
    if (oldcount && ptr) {
        v_copy_channels(newptr, ptr, channels, oldcount < count ? oldcount : count);
    } 
    if (ptr) deallocate_channels<T>(ptr, channels);
    return newptr;
}

/// RAII class to call deallocate() on destruction
template <typename T>
class Deallocator
{
public:
    Deallocator(T *t) : m_t(t) { }
    ~Deallocator() { deallocate<T>(m_t); }
private:
    T *m_t;
};

}

#endif

