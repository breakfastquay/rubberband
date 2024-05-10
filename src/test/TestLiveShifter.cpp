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

#ifndef BOOST_TEST_DYN_LINK
#define BOOST_TEST_DYN_LINK
#endif
#include <boost/test/unit_test.hpp>

#include "../../rubberband/RubberBandLiveShifter.h"

#include <iostream>

#include <cmath>

using namespace RubberBand;

using std::vector;
using std::cerr;
using std::endl;

namespace tt = boost::test_tools;

BOOST_AUTO_TEST_SUITE(TestLiveShifter)

static void dump(const vector<float> &in,
                 const vector<float> &out,
                 const vector<float> &expected,
                 int delay)
{
    std::cerr << "dump: delay reported as " << delay << std::endl;
    
    // The prefix is to allow grep on the test output
        
    std::cout << "IN,sample,V" << std::endl;
    for (int i = 0; i < int(in.size()); ++i) {
        std::cout << "IN," << i << "," << in[i] << std::endl;
    }
        
    std::cout << "OUT,sample,V" << std::endl;
    for (int i = 0; i < int(out.size()); ++i) {
        std::cout << "OUT," << i << "," << out[i] << std::endl;
    }

    std::cout << "SHIFTED,sample,V" << std::endl;
    for (int i = 0; i + delay < int(out.size()); ++i) {
        std::cout << "SHIFTED," << i << "," << out[i + delay] << std::endl;
    }
    
    std::cout << "EXPECTED,sample,V" << std::endl;
    for (int i = 0; i < int(expected.size()); ++i) {
        std::cout << "EXPECTED," << i << "," << expected[i] << std::endl;
    }
    
    std::cout << "DIFF,sample,V" << std::endl;
    for (int i = 0; i + delay < int(expected.size()); ++i) {
        std::cout << "DIFF," << i << "," << out[i + delay] - expected[i] << std::endl;
    }
}

static void check_sinusoid_unchanged(int n, int rate, float freq,
                                     RubberBandLiveShifter::Options options,
                                     bool printDebug)
{
    if (printDebug) {
        RubberBandLiveShifter::setDefaultDebugLevel(2);
    }
    
    RubberBandLiveShifter shifter(rate, 1, options);
    
    int blocksize = shifter.getBlockSize();
    BOOST_TEST(blocksize == 512);

    n = (n / blocksize + 1) * blocksize;
    
    vector<float> in(n), out(n);
    for (int i = 0; i < n; ++i) {
        in[i] = 0.5f * sinf(float(i) * freq * M_PI * 2.f / float(rate));
    }

    for (int i = 0; i < n; i += blocksize) {
        float *inp = in.data() + i;
        float *outp = out.data() + i;
        shifter.shift(&inp, &outp);
    }

    int delay = shifter.getStartDelay();
    
    // We now have n samples of a simple sinusoid with stretch factor
    // 1.0; obviously we expect the output to be essentially the same
    // thing. It will have lower precision for a while at the start,
    // so we check that with a threshold of 0.1; after that we expect
    // better precision.

    int slackpart = 2048;
    float slackeps = 1.0e-1f;
    float eps = 1.0e-3f;

#ifdef USE_BQRESAMPLER
    eps = 1.0e-2f;
#endif
    
    for (int i = 0; i < slackpart; ++i) {
        float fin = in[i];
        float fout = out[delay + i];
        float err = fabsf(fin - fout);
        if (err > slackeps) {
            std::cerr << "Error at index " << i << " exceeds slack eps "
                      << slackeps << ": output " << fout << " - input "
                      << fin << " = " << fout - fin << std::endl;
            BOOST_TEST(err < eps);
            break;
        }
    }
    
    for (int i = slackpart; i < n - delay; ++i) {
        float fin = in[i];
        float fout = out[delay + i];
        float err = fabsf(fin - fout);
        if (err > eps) {
            std::cerr << "Error at index " << i << " exceeds tight eps "
                      << eps << ": output " << fout << " - input "
                      << fin << " = " << fout - fin << std::endl;
            BOOST_TEST(err < eps);
            break;
        }
    }

    if (printDebug) {
        RubberBandLiveShifter::setDefaultDebugLevel(0);
        dump(in, out, in, delay);
    }
}

static void check_sinusoid_shifted(int n, int rate, float freq, float shift,
                                   RubberBandLiveShifter::Options options,
                                   bool printDebug)
{
    if (printDebug) {
        RubberBandLiveShifter::setDefaultDebugLevel(2);
    }
    
    RubberBandLiveShifter shifter(rate, 1, options);

    shifter.setPitchScale(shift);
    
    int blocksize = shifter.getBlockSize();
    BOOST_TEST(blocksize == 512);

    n = (n / blocksize + 1) * blocksize;
    
    vector<float> in(n), out(n), expected(n);
    for (int i = 0; i < n; ++i) {
        in[i] = 0.5f * sinf(float(i) * freq * M_PI * 2.f / float(rate));
        expected[i] = 0.5f * sinf(float(i) * freq * shift * M_PI * 2.f / float(rate));
    }

    for (int i = 0; i < n; i += blocksize) {
        float *inp = in.data() + i;
        float *outp = out.data() + i;
        shifter.shift(&inp, &outp);
    }

    int delay = shifter.getStartDelay();

    std::cerr << "delay reported as " << delay << std::endl;    
    
    // We now have n samples of a simple sinusoid with stretch factor
    // 1.0; obviously we expect the output to be essentially the same
    // thing. It will have lower precision for a while at the start,
    // so we check that with a threshold of 0.1; after that we expect
    // better precision.

    int slackpart = 2048;
    float slackeps = 1.0e-1f;
    float eps = 1.0e-3f;

#ifdef USE_BQRESAMPLER
    eps = 1.0e-2f;
#endif
    
    for (int i = 0; i < slackpart; ++i) {
        float fin = expected[i];
        float fout = out[delay + i];
        float err = fabsf(fin - fout);
        if (err > slackeps) {
            std::cerr << "Error at index " << i << " exceeds slack eps "
                      << slackeps << ": output " << fout << " - expected "
                      << fin << " = " << fout - fin << std::endl;
            BOOST_TEST(err < eps);
            break;
        }
    }
    
    for (int i = slackpart; i < n - delay; ++i) {
        float fin = expected[i];
        float fout = out[delay + i];
        float err = fabsf(fin - fout);
        if (err > eps) {
            std::cerr << "Error at index " << i << " exceeds tight eps "
                      << eps << ": output " << fout << " - expected "
                      << fin << " = " << fout - fin << std::endl;
            BOOST_TEST(err < eps);
            break;
        }
    }

    if (printDebug) {
        RubberBandLiveShifter::setDefaultDebugLevel(0);
        dump(in, out, expected, delay);
    }
}

BOOST_AUTO_TEST_CASE(sinusoid_unchanged)
{
    int n = 20000;
    check_sinusoid_unchanged(n, 44100, 440.f, 0, false);
    check_sinusoid_unchanged(n, 48000, 260.f, 0, false);
}

BOOST_AUTO_TEST_CASE(sinusoid_down_octave)
{
    int n = 20000;
    check_sinusoid_shifted(n, 44100, 440.f, 0.5f, 0, true);
//    check_sinusoid_shifted(n, 48000, 260.f, 0.5f, 0, false);
}

BOOST_AUTO_TEST_SUITE_END()
