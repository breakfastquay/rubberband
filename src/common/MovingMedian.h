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

#ifndef RUBBERBAND_MOVING_MEDIAN_H
#define RUBBERBAND_MOVING_MEDIAN_H

#include "SampleFilter.h"
#include "FixedVector.h"
#include "Allocators.h"

#include <algorithm>
#include <iostream>

//#define DEBUG_MM 1

namespace RubberBand
{

template <typename T>
class MovingMedianStack
{
public:
    MovingMedianStack(int nfilters, int filterLength, float percentile = 50.f) :
        m_buffer(nfilters * filterLength * 2, {}),
        m_length(filterLength)
    {
        setPercentile(percentile);
    }

    ~MovingMedianStack() { 
    }

    int getNFilters() const {
        return m_buffer.size() / (m_length * 2);
    }
    
    int getSize() const {
        return m_length;
    }
    
    void setPercentile(float p) {
        m_index = int((m_length * p) / 100.f);
        if (m_index >= m_length) m_index = m_length-1;
        if (m_index < 0) m_index = 0;
    }

    void push(int filter, T value) {
        if (value != value) {
            std::cerr << "WARNING: MovingMedian: NaN encountered" << std::endl;
            value = T();
        }
        T *frame = frameFor(filter);
        T toDrop = frame[0];
	v_move(frame, frame+1, m_length-1);
	frame[m_length-1] = value;
        dropAndPut(filter, toDrop, value);
    }

    T get(int filter) const {
        const T *sorted = sortedFor(filter);
	return sorted[m_index];
    }

    void reset() {
	v_zero(m_buffer.data(), m_buffer.size());
    }
    
private:
    FixedVector<T> m_buffer;
    int m_length;
    int m_index;

    const T *frameFor(int filter) const {
        return m_buffer.data() + filter * m_length * 2;
    }
    T *frameFor(int filter) {
        return m_buffer.data() + filter * m_length * 2;
    }
    const T *sortedFor(int filter) const {
        return frameFor(filter) + m_length;
    }
    T *sortedFor(int filter) {
        return frameFor(filter) + m_length;
    }

    void dropAndPut(int filter, const T &toDrop, const T &toPut) {
	// precondition: sorted contains m_length values, one of which is toDrop
	// postcondition: sorted contains m_length values, one of which is toPut
        // (and one instance of toDrop has been removed)
        const int n = m_length;
        T *sorted = sortedFor(filter);
        int dropIx;
        if (toDrop <= *sorted) {
            // this is quite a common short-circuit in situations
            // where many values can be (the equivalent of) 0
            dropIx = 0;
        } else {
            dropIx = std::lower_bound(sorted, sorted + n, toDrop) - sorted;
        }

#ifdef DEBUG_MM
        std::cout << "\nbefore: [";
        for (int i = 0; i < n; ++i) {
            if (i > 0) std::cout << ",";
            std::cout << sorted[i];
        }
        std::cout << "]" << std::endl;

        std::cout << "toDrop = " << toDrop << ", dropIx = " << dropIx << std::endl;
        std::cout << "toPut = " << toPut << std::endl;
        if (sorted[dropIx] != toDrop) {
            throw std::runtime_error("element not found");
        }
#endif
        
        if (toPut > toDrop) {
            int i = dropIx;
            while (i+1 < n) {
                if (sorted[i+1] > toPut) {
                    break;
                }
                sorted[i] = sorted[i+1];
                ++i;
            }
            sorted[i] = toPut;
        } else if (toPut < toDrop) {
            int i = dropIx;
            while (true) {
                if (--i < 0 || sorted[i] < toPut) {
                    break;
                }
                sorted[i+1] = sorted[i];
            }
            sorted[i+1] = toPut;
        }

#ifdef DEBUG_MM
        std::cout << "after: [";
        for (int i = 0; i < n; ++i) {
            if (i > 0) std::cout << ",";
            std::cout << sorted[i];
        }
        std::cout << "]" << std::endl;

        if (!std::is_sorted(sorted, sorted + n)) {
            throw std::runtime_error("array is not sorted");
        }
#endif
    }

    MovingMedianStack(const MovingMedianStack &) =delete;
    MovingMedianStack &operator=(const MovingMedianStack &) =delete;
};

template <typename T>
class MovingMedian : public SampleFilter<T>
{
public:
    MovingMedian(int size, float percentile = 50.f) :
        m_mm(1, size, percentile)
    {
    }

    ~MovingMedian() {
    }

    int getSize() const {
        return m_mm.getSize();
    }

    void setPercentile(float p) {
        m_mm.setPercentile(p);
    }

    void push(T value) {
        m_mm.push(0, value);
    }

    T get() const {
        return m_mm.get(0);
    }

    void reset() {
        m_mm.reset();
    }

    // Convenience function that applies a given filter to an array
    // in-place. Array has length n. Modifies both the filter and the
    // array.
    //
    static void filter(MovingMedian<T> &mm, T *v, int n) {
        int fn = mm.getSize();
        int lag = fn / 2;
        mm.reset();
        int i = 0;
        for (; i < lag; ++i) {
            if (i < n) mm.push(v[i]);
        }
        for (; i < n; ++i) {
            mm.push(v[i]);
            v[i-lag] = mm.get();
        }
        for (; i < lag; ++i) {
            // just for the unusual case where lag > n
            mm.push(T());
            (void)mm.get();
        }
        for (; i < n + lag; ++i) {
            mm.push(T());
            v[i-lag] = mm.get();
        }
    }

    // As above but with a vector argument
    //
    static void filter(MovingMedian<T> &mm, std::vector<T> &v) {
        filter(mm, v.data(), v.size());
    }
    
private:
    MovingMedianStack<T> m_mm;

    MovingMedian(const MovingMedian &) =delete;
    MovingMedian &operator=(const MovingMedian &) =delete;
};

}

#endif

