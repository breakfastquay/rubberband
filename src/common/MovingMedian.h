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
#include "SingleThreadRingBuffer.h"

#include <algorithm>
#include <iostream>

//#define DEBUG_MM 1

namespace RubberBand
{

template <typename T>
class MovingMedian : public SampleFilter<T>
{
public:
    MovingMedian(int filterLength, float percentile = 50.f) :
        m_buffer(filterLength),
        m_sortspace(filterLength, {})
    {
        setPercentile(percentile);
    }

    ~MovingMedian() { 
    }

    MovingMedian(const MovingMedian &) =default;
    MovingMedian &operator=(const MovingMedian &) =default;

    int getSize() const {
        return m_buffer.getSize();
    }
    
    void setPercentile(float p) {
        int length = getSize();
        m_index = int((length * p) / 100.f);
        if (m_index >= length) m_index = length-1;
        if (m_index < 0) m_index = 0;
    }

    void push(T value) {
        if (value != value) {
            std::cerr << "WARNING: MovingMedian: NaN encountered" << std::endl;
            value = T();
        }
        if (m_buffer.getWriteSpace() == 0) {
            T toDrop = m_buffer.readOne();
            dropAndPut(toDrop, value);
            m_buffer.writeOne(value);
        } else {
            put(value);
            m_buffer.writeOne(value);
        }
    }

    T get() const {
	return m_sortspace[m_index];
    }

    void reset() {
        m_buffer.reset();
	v_zero(m_sortspace.data(), m_sortspace.size());
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
    SingleThreadRingBuffer<T> m_buffer;
    std::vector<T> m_sortspace;
    int m_index;

    void dropAndPut(const T &toDrop, const T &toPut) {
	// precondition: sorted contains m_length values, one of which is toDrop
	// postcondition: sorted contains m_length values, one of which is toPut
        // (and one instance of toDrop has been removed)

        // This implementation was timed for rather short filters (no
        // longer than maybe 16 items). Two binary searches plus a
        // memmove should be faster for longer ones.
        
        const int n = getSize();
        T *sorted = m_sortspace.data();
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
            throw std::logic_error("element not found");
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
            throw std::logic_error("array is not sorted");
        }
#endif
    }

    void put(const T &toPut) {
	// precondition: sorted contains fewer than m_length values,
	// packed at the start
	// postcondition: sorted contains up to m_length values,
        // packed at the start, one of which is toPut
        const int n = m_buffer.getReadSpace(); // items in sorted

#ifdef DEBUG_MM
        if (n >= m_length) {
            throw std::logic_error("length mismatch");
        }
#endif
            
        T *sorted = m_sortspace.data();
        int putIx = std::lower_bound(sorted, sorted + n, toPut) - sorted;

#ifdef DEBUG_MM
        std::cout << "\nbefore: [";
        for (int i = 0; i < n; ++i) {
            if (i > 0) std::cout << ",";
            std::cout << sorted[i];
        }
        std::cout << "]" << std::endl;

        std::cout << "toPut = " << toPut << ", putIx = " << putIx << std::endl;
#endif

        if (putIx < n) {
            v_move(sorted + putIx + 1, sorted + putIx, n - putIx);
        }
        sorted[putIx] = toPut;

#ifdef DEBUG_MM
        std::cout << "after: [";
        for (int i = 0; i < n + 1; ++i) {
            if (i > 0) std::cout << ",";
            std::cout << sorted[i];
        }
        std::cout << "]" << std::endl;

        if (!std::is_sorted(sorted, sorted + n)) {
            throw std::logic_error("array is not sorted");
        }
#endif
    }
};

template <typename T>
class MovingMedianStack
{
public:
    MovingMedianStack(int nfilters, int size, float percentile = 50.f) :
        m_stack(nfilters, { size, percentile })
    {
    }

    ~MovingMedianStack() {
    }

    int getSize() const {
        return m_stack[0].getSize();
    }

    void setPercentile(float p) {
        for (auto &f: m_stack) {
            f.setPercentile(p);
        }
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

