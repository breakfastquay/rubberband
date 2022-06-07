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
#include "Allocators.h"

#include <algorithm>
#include <iostream>

//#define DEBUG_MM 1

namespace RubberBand
{

template <typename T>
class MovingMedian : public SampleFilter<T>
{
    typedef SampleFilter<T> P;

public:
    MovingMedian(int size, float percentile = 50.f) :
        SampleFilter<T>(size),
	m_frame(allocate_and_zero<T>(size * 2)),
	m_sorted(m_frame + size)
    {
        setPercentile(percentile);
    }

    ~MovingMedian() { 
	deallocate(m_frame);
    }

    void setPercentile(float p) {
        m_index = int((P::m_size * p) / 100.f);
        if (m_index >= P::m_size) m_index = P::m_size-1;
        if (m_index < 0) m_index = 0;
    }

    void push(T value) {
        if (value != value) {
            std::cerr << "WARNING: MovingMedian: NaN encountered" << std::endl;
            value = T();
        }
        T toDrop = m_frame[0];
	v_move(m_frame, m_frame+1, P::m_size-1);
	m_frame[P::m_size-1] = value;
        dropAndPut(toDrop, value);
    }

    T get() const {
	return m_sorted[m_index];
    }

    void reset() {
	v_zero(m_frame, P::m_size);
	v_zero(m_sorted, P::m_size);
    }

    // Convenience function that applies a given filter to an array
    // in-place. Array must have length equal to getSize(). Modifies
    // both the filter and the array.
    //
    static void filter(MovingMedian<T> &mm, T *v) {
        int n = mm.getSize();
        int lag = n / 2;
        mm.reset();
        for (int i = 0; i < lag; ++i) {
            mm.push(v[i]);
        }
        for (int i = lag; i < n; ++i) {
            mm.push(v[i]);
            v[i-lag] = mm.get();
        }
        for (int i = n; i < n + lag; ++i) {
            mm.push(T());
            v[i-lag] = mm.get();
        }
    }
    
private:
    T *const m_frame;
    T *const m_sorted;
    int m_index;

    void dropAndPut(const T &toDrop, const T &toPut) {
	// precondition: m_sorted contains m_size values, one of which is toDrop
	// postcondition: m_sorted contains m_size values, one of which is toPut
        // (and one instance of toDrop has been removed)
        int n = P::m_size;
        int dropIx;
        if (toDrop <= m_sorted[0]) {
            // this is quite a common short-circuit in situations
            // where many values can be (the equivalent of) 0
            dropIx = 0;
        } else {
            dropIx = std::lower_bound(m_sorted, m_sorted + n, toDrop) - m_sorted;
        }

#ifdef DEBUG_MM
        std::cout << "\nbefore: [";
        for (int i = 0; i < P::m_size; ++i) {
            if (i > 0) std::cout << ",";
            std::cout << m_sorted[i];
        }
        std::cout << "]" << std::endl;

        std::cout << "toDrop = " << toDrop << ", dropIx = " << dropIx << std::endl;
        std::cout << "toPut = " << toPut /* << ", putIx = " << putIx */ << std::endl;
        if (m_sorted[dropIx] != toDrop) {
            throw std::runtime_error("element not found");
        }
#endif
        
        if (toPut > toDrop) {
            int i = dropIx;
            while (i+1 < n) {
                if (m_sorted[i+1] > toPut) {
                    break;
                }
                m_sorted[i] = m_sorted[i+1];
                ++i;
            }
            m_sorted[i] = toPut;
        } else if (toPut < toDrop) {
            int i = dropIx;
            while (i >= 0) {
                --i;
                if (i < 0 || m_sorted[i] < toPut) {
                    break;
                }
                m_sorted[i+1] = m_sorted[i];
            }
            m_sorted[i+1] = toPut;
        }

#ifdef DEBUG_MM
        std::cout << "after: [";
        for (int i = 0; i < P::m_size; ++i) {
            if (i > 0) std::cout << ",";
            std::cout << m_sorted[i];
        }
        std::cout << "]" << std::endl;

        if (!std::is_sorted(m_sorted, m_sorted + n)) {
            throw std::runtime_error("array is not sorted");
        }
#endif
    }

    MovingMedian(const MovingMedian &) =delete;
    MovingMedian &operator=(const MovingMedian &) =delete;
};

}

#endif

