/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band Library
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2024 Particular Programs Ltd.

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

#ifndef RUBBERBAND_MOVING_MEDIAN_H
#define RUBBERBAND_MOVING_MEDIAN_H

#include "SampleFilter.h"
#include "FixedVector.h"
#include "Allocators.h"
#include "SingleThreadRingBuffer.h"

#include <algorithm>
#include <iostream>
#include <set>

namespace RubberBand
{

/**
 * A median or modal filter implemented using sorting, for continuous
 * values. Pushing a value updates the sorted record, after which you
 * can get() the median.  You can call drop() to drop the oldest value
 * without pushing a new one, for example to drain the filter at the
 * tail of the sequence.
 */
template <typename T>
class MovingMedian : public SampleFilter<T>
{
    static constexpr float fifty = 50.f;
    
public:
    MovingMedian(int filterLength, float percentile = fifty) :
        m_buffer(filterLength),
        m_fill(0),
        m_percentile(percentile)
    {
        m_sortspace.clear();
    }

    ~MovingMedian() { }

    MovingMedian(const MovingMedian &) =default;
    MovingMedian &operator=(const MovingMedian &) =default;

    int getSize() const {
        return m_buffer.getSize();
    }

    void setPercentile(float p) {
        m_percentile = p;
    }
    
    void push(T value) {
        if (std::isnan(value)) {
            std::cerr << "WARNING: MovingMedian: NaN encountered" << std::endl;
            value = T();
        }
        if (m_fill == getSize()) {
            T toDrop = m_buffer.readOne();
            dropAndPut(toDrop, value);
        } else {
            put(value);
        }
        m_buffer.writeOne(value);
    }

    void drop() {
        if (m_fill > 0) {
            T toDrop = m_buffer.readOne();
            drop(toDrop);
        }
    }
    
    /** Return the median of the values in the filter currently. If
     *  the median lies between two values, return the first of them.
     */
    T get() const {
        auto it = m_sortspace.begin();
        std::advance(it, std::distance(m_sortspace.begin(), m_sortspace.end()) / 2);
        return *it;
    }

    void reset() {
        m_buffer.reset();
        m_sortspace.clear(); 
        m_fill = 0;
    }

    // Convenience function that applies a given filter to an array
    // in-place. Array has length n. Modifies both the filter and the
    // array.
    //
    static void filter(MovingMedian<T> &f, T *v, int n) {
        f.reset();
        int flen = f.getSize();
        int i = -(flen / 2);
        int j = 0;
        while (i != n) {
            if (j < n) {
                f.push(v[j]);
            } else if (j >= flen) {
                f.drop();
            }
            if (i >= 0) {
                v[i] = f.get();
            }
            ++i;
            ++j;
        }
    }

    // As above but with a vector argument
    //
    static void filter(MovingMedian<T> &mm, std::vector<T> &v) {
        filter(mm, v.data(), v.size());
    }
    
private:
    SingleThreadRingBuffer<T> m_buffer;
    std::multiset<T> m_sortspace;
    int m_fill;
    float m_percentile;

    void dropAndPut(const T &toDrop, const T &toPut) {
        m_sortspace.erase(m_sortspace.find(toDrop));
        m_sortspace.insert(toPut);
    }

    void put(const T &toPut) {
        m_sortspace.insert(toPut);
        ++m_fill;
    }

    void drop(const T &toDrop) {
        m_sortspace.erase(m_sortspace.find(toDrop));
        --m_fill;
    }

};

template <typename T>
class MovingMedianStack
{
public:
    MovingMedianStack(int nfilters, int size) :
        m_stack(nfilters, { size })
    {
    }

    ~MovingMedianStack() {
    }

    int getSize() const {
        return m_stack[0].getSize();
    }

    void push(int filter, T value) {
        m_stack[filter].push(value);
    }

    T get(int filter) const {
        return m_stack[filter].get();
    }

    void reset() {
        for (auto &f: m_stack) {
            f.reset();
        }
    }
    
private:
    std::vector<MovingMedian<T>> m_stack;
};

}

#endif

