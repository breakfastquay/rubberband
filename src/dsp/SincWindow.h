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

#ifndef _RUBBERBAND_SINC_WINDOW_H_
#define _RUBBERBAND_SINC_WINDOW_H_

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <map>

#include "system/sysutils.h"
#include "system/VectorOps.h"
#include "system/Allocators.h"

namespace RubberBand {

template <typename T>
class SincWindow
{
public:
    /**
     * Construct a sinc windower which produces a window of the given
     * size containing the values of sinc(x) with x=0 at the centre,
     * such that the distance from -pi to pi (the point at which the
     * sinc function first crosses zero, for negative and positive
     * arguments respectively) is p samples.
     */
    SincWindow(int size, int p) : m_size(size), m_p(p) { encache(); }
    SincWindow(const SincWindow &w) : m_size(w.m_size), m_p(w.m_p) { encache(); }
    SincWindow &operator=(const SincWindow &w) {
	if (&w == this) return *this;
	m_size = w.m_size;
	m_p = w.m_p;
	encache();
	return *this;
    }
    virtual ~SincWindow() { delete[] m_cache; }
    
    inline void cut(T *const R__ block) const {
        v_multiply(block, m_cache, m_size);
    }

    inline void cut(const T *const R__ src, T *const R__ dst) const {
        v_multiply(dst, src, m_cache, m_size);
    }

    inline void add(T *const R__ dst, T scale) const {
        v_add_with_gain(dst, m_cache, m_size, scale);
    }

    inline T getArea() const { return m_area; }
    inline T getValue(int i) const { return m_cache[i]; }

    inline int getSize() const { return m_size; }
    inline int getP() const { return m_p; }

protected:
    int m_size;
    int m_p;
    T *R__ m_cache;
    T m_area;
    
    void encache();
};

template <typename T>
void SincWindow<T>::encache()
{
    const int n = m_size;
    T *mult = allocate<T>(n);
    v_set(mult, T(1.0), n);
    int i;

    for (i = 0; i < n; ++i) {
	T extent = T(n)/2.;
	T arg = (T(i) - extent) * (2. * M_PI) / m_p;
	if (arg != 0.) {
	    mult[i] *= sin(arg) / arg;
	}
    }
	
    m_cache = mult;

    m_area = 0;
    for (i = 0; i < n; ++i) {
	std::cout << i << ":" << m_cache[i] << " ";
        m_area += m_cache[i];
    }
    std::cout << std::endl;
    m_area /= n;
}

}

#endif
