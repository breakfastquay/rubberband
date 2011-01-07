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

#ifndef _SAMPLE_FILTER_H_
#define _SAMPLE_FILTER_H_

#include <cassert>

namespace RubberBand
{

template <typename T>
class SampleFilter
{
public:
    SampleFilter(int size) : m_size(size) {
	assert(m_size > 0);
    }

    virtual ~SampleFilter() { }

    int getSize() const { return m_size; }

    virtual void push(T) = 0;
    virtual T get() const = 0;
    virtual void reset() = 0;

protected:
    const int m_size;

private:
    SampleFilter(const SampleFilter &);
    SampleFilter &operator=(const SampleFilter &);
};

}

#endif

