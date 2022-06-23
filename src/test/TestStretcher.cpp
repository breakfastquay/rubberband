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

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "../../rubberband/RubberBandStretcher.h"

#include <iostream>

#include <cmath>

using namespace RubberBand;
using namespace std;
namespace tt = boost::test_tools;

BOOST_AUTO_TEST_SUITE(TestStretcher)

BOOST_AUTO_TEST_CASE(sinusoid_unchanged_single_offline)
{
    int n = 10000;
    float freq = 440.f;
    int rate = 44100;
    RubberBandStretcher stretcher(rate, 1);
    vector<float> in(n), out(n);
    for (int i = 0; i < n; ++i) {
        in[i] = sinf(float(i) * freq * M_PI * 2.f / float(rate));
    }
    float *inp = in.data(), *outp = out.data();

    stretcher.setMaxProcessSize(n);
    stretcher.setExpectedInputDuration(n);
    BOOST_TEST(stretcher.available() == 0);

    stretcher.study(&inp, n, true);
    BOOST_TEST(stretcher.available() == 0);

    stretcher.process(&inp, n, true);
    BOOST_TEST(stretcher.available() == n);

    size_t got = stretcher.retrieve(&outp, n);
    BOOST_TEST(got == n);
    BOOST_TEST(stretcher.available() == -1);

    // We now have n samples of a simple sinusoid with stretch factor
    // 1.0; obviously we expect the output to be essentially the same
    // thing. It will have lower precision for a while at the start
    // and end because of windowing factors, so we check those with a
    // threshold of 0.1; in the middle we expect better precision.

    // This syntax for comparing containers with a certain tolerance
    // using BOOST_TEST is just bonkers. I can't find the << syntax to
    // combine manipulators documented anywhere other than in a
    // release note, but it does work. Well, sort of - it works this
    // way around but not as per_element << tolerance. And
    // tolerance(0.1) doesn't do what you'd expect if the things
    // you're comparing are floats (it sets the tolerance for doubles,
    // leaving float comparison unchanged). Clever... too clever.
    
    BOOST_TEST(in == out,
               tt::tolerance(0.1f) << tt::per_element());
    
    BOOST_TEST(vector<float>(in.begin() + 1000, in.begin() + n - 1000) ==
               vector<float>(out.begin() + 1000, out.begin() + n - 1000),
               tt::tolerance(0.0001f) << tt::per_element());
}

BOOST_AUTO_TEST_SUITE_END()
