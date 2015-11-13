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

#include "rubberband/RubberBandStretcher.h"

#include <iostream>
#include <cstdlib>
#include <cmath>

using namespace std;
using namespace RubberBand;

int main(int argc, char **argv)
{
    if (argc != 3) {
	//!!! + other options: window, hq, sample rate
	cerr << "usage: " << argv[0] << " timeratio pitchshift" << endl;
	return 2;
    }

    double ratio = atof(argv[1]); 
    double pitchshift = atof(argv[2]);

    double frequencyshift = 1.0;
    if (pitchshift != 0.0) {
        frequencyshift *= pow(2.0, pitchshift / 12);
    }
    
    RubberBandStretcher ts(44100, 1,
			   RubberBandStretcher::OptionProcessRealTime,
                           ratio, frequencyshift);
    
    int reported = ts.getLatency();
    int required = ts.getSamplesRequired();
    int latency = -1;
    int delay = -1;
    float prev = 0.f;

    cout << "reported latency: " << reported << endl;
    cout << "initial required: " << required << endl;

    for (int i = 0; ; ++i) {
	float f = 0.f;
        float *ff = &f;
        if (i == 0) f = 1.f;
	if (latency < 0 && ts.available() > 0) {
	    latency = i;
            cout << "measured latency: " << latency << endl;
	}
	ts.process(&ff, 1, false);
	if (ts.available() > 0) {
	    ts.retrieve(&ff, 1);
	    if (fabsf(f) < prev) {
		//!!! is this right, or do we have an off-by-one error?
		delay = i - latency;
                cout << "measured delay: " << delay << endl;
		break;
	    }
	    prev = fabsf(f);
	}
    }

    return 0;
}
