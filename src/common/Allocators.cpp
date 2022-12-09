/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band Library
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2022 Particular Programs Ltd.

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

#include "Allocators.h"

#ifdef HAVE_IPP
#include <ipps.h>
#endif

#include <iostream>
using std::cerr;
using std::endl;

namespace RubberBand {

#ifdef HAVE_IPP

template <>
float *allocate(size_t count)
{
    float *ptr = ippsMalloc_32f(count);
    if (!ptr) throw (std::bad_alloc());
    return ptr;
}

template <>
double *allocate(size_t count)
{
    double *ptr = ippsMalloc_64f(count);
    if (!ptr) throw (std::bad_alloc());
    return ptr;
}

template <>
void deallocate(float *ptr)
{
    if (ptr) ippsFree((void *)ptr);
}

template <>
void deallocate(double *ptr)
{
    if (ptr) ippsFree((void *)ptr);
}

#endif

}

