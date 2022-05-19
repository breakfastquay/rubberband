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

#ifndef RUBBERBAND_BIN_SEGMENTER_H
#define RUBBERBAND_BIN_SEGMENTER_H

#include "BinClassifier.h"

#include <vector>

namespace RubberBand {

class BinSegmenter
{
public:
    struct Segmentation {
        double percussiveBelow;
        double percussiveAbove;
        double residualAbove;
        Segmentation(double _pb, double _pa, double _ra) :
            percussiveBelow(_pb), percussiveAbove(_pa), residualAbove(_ra) { }
    };

    struct Parameters {
        int fftSize;
        double sampleRate;
        Parameters(int _fftSize, double _sampleRate) :
            fftSize(_fftSize), sampleRate(_sampleRate) { }
    };
    
    BinSegmenter(Parameters parameters,
                 BinClassifier::Parameters classifierParameters) :
        m_parameters(parameters),
        m_classifierParameters(classifierParameters),
        m_classifier(classifierParameters),
        m_classification(classifierParameters.binCount,
                         BinClassifier::Classification::Silent),
        m_numeric(classifierParameters.binCount, 0),
        m_classFilter(classifierParameters.binCount / 64)
    {
    }

    Segmentation segment(const float *const mag) {
        int n = m_classifierParameters.binCount;
        m_classifier.classify(mag, m_classification.data());
        for (int i = 0; i < n; ++i) {
            switch (m_classification[i]) {
            case BinClassifier::Classification::Harmonic:
                m_numeric[i] = 0; break;
            case BinClassifier::Classification::Percussive:
                m_numeric[i] = 1; break;
            default:
                m_numeric[i] = 2; break;
            }
        }
        MovingMedian<int>::filter(m_classFilter, m_numeric.data());
        double f0 = 0.0;
        for (int i = 1; i < n; ++i) {
            if (m_numeric[i] != 1) {
                f0 = frequencyForBin(i);
                break;
            }
        }
        double nyquist = m_parameters.sampleRate / 2.0;
        int top = binForFrequency(16000.0);
        if (top >= n) top = n-1;
        double f1 = nyquist;
        double f2 = nyquist;
        bool inPercussive = false;
        for (int i = top; i > 0; --i) {
            if (m_numeric[i] == 1) { // percussive
                if (!inPercussive) {
                    inPercussive = true;
                    f2 = frequencyForBin(i);
                    continue;
                }
            } else if (m_numeric[i] == 0) { // harmonic
                if (inPercussive) {
                    f1 = frequencyForBin(i);
                }
                break; // always when harmonic reached
            }
        }
        return Segmentation(f0, f1, f2);
    }

protected:
    Parameters m_parameters;
    BinClassifier::Parameters m_classifierParameters;
    BinClassifier m_classifier;
    std::vector<BinClassifier::Classification> m_classification;
    std::vector<int> m_numeric;
    MovingMedian<int> m_classFilter;

    int binForFrequency(double f) {
        return int(round(f * double(m_parameters.fftSize) /
                         m_parameters.sampleRate));
    }
    double frequencyForBin(int b) {
        return (double(b) * m_parameters.sampleRate)
            / double(m_parameters.fftSize);
    }
    
    BinSegmenter(const BinSegmenter &) =delete;
    BinSegmenter &operator=(const BinSegmenter &) =delete;
};

}

#endif
