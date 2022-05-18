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

#ifndef RUBBERBAND_BIN_CLASSIFIER_H
#define RUBBERBAND_BIN_CLASSIFIER_H

#include "../system/Allocators.h"
#include "../dsp/MovingMedian.h"
#include "../base/RingBuffer.h"

#include <vector>
#include <memory>

namespace RubberBand {

class BinClassifier {

    enum class Classification {
        Harmonic = 0,
        Percussive = 1,
        Residual = 2,
        Silent = 3
    };

    struct Parameters {
        int binCount;
        int horizontalFilterLength;
        int horizontalFilterLag;
        int verticalFilterLength;
        double harmonicThreshold;
        double percussiveThreshold;
        float silenceThreshold;
        Parameters(int _binCount, int _horizontalFilterLength,
                   int _horizontalFilterLag, int _verticalFilterLength,
                   double _harmonicThreshold, double _percussiveThreshold,
                   float _silenceThreshold) :
            binCount(_binCount),
            horizontalFilterLength(_horizontalFilterLength),
            horizontalFilterLag(_horizontalFilterLag),
            verticalFilterLength(_verticalFilterLength),
            harmonicThreshold(_harmonicThreshold),
            percussiveThreshold(_percussiveThreshold),
            silenceThreshold(_silenceThreshold) { }
    };
    
    BinClassifier(Parameters parameters) :
        m_parameters(parameters),
        m_vfQueue(parameters.horizontalFilterLag)
    {
        int n = m_parameters.binCount;

        for (int i = 0; i < n; ++i) {
            m_hFilters.push_back(std::make_shared<MovingMedian<float>>
                                 (m_parameters.horizontalFilterLength));
        }

        m_vFilter = std::make_unique<MovingMedian<float>>
            (m_parameters.verticalFilterLength);

        m_hf = allocate_and_zero<float>(n);
        m_vf = allocate_and_zero<float>(n);
        
        for (int i = 0; i < m_parameters.horizontalFilterLag; ++i) {
            float *entry = allocate_and_zero<float>(n);
            m_vfQueue.write(&entry, 1);
        }
    }

    ~BinClassifier()
    {
        while (m_vfQueue.getReadSpace() > 0) {
            float *entry = m_vfQueue.readOne();
            deallocate(entry);
        }

        deallocate(m_hf);
        deallocate(m_vf);
    }

    void classify(const float *const mag, Classification *classification) {
        const int n = m_parameters.binCount;

        for (int i = 0; i < n; ++i) {
            m_hFilters[i]->push(mag[i]);
            m_hf[i] = m_hFilters[i]->get();
        }
        
        m_vFilter->reset();
        int vFilterLag = m_parameters.verticalFilterLength / 2;
        
        for (int i = 0; i < vFilterLag; ++i) {
            m_vFilter->push(mag[i]);
        }
        for (int i = vFilterLag; i < n; ++i) {
            m_vFilter->push(mag[i]);
            m_vf[i-vFilterLag] = m_vFilter->get();
        }
        for (int i = n; i < n + vFilterLag; ++i) {
            m_vFilter->push(0.f);
            m_vf[i-vFilterLag] = m_vFilter->get();
        }

        if (m_parameters.horizontalFilterLag > 0) {
            float *lagged = m_vfQueue.readOne();
            m_vfQueue.write(&m_vf, 1);
            m_vf = lagged;
        }

        double eps = 1.0e-7;
        
        for (int i = 0; i < n; ++i) {
            Classification c;
            if (mag[i] < m_parameters.silenceThreshold) {
                c = Classification::Silent;
            } else if (double(m_hf[i]) / (double(m_vf[i]) + eps) >
                       m_parameters.harmonicThreshold) {
                c = Classification::Harmonic;
            } else if (double(m_vf[i]) / (double(m_hf[i]) + eps) >
                       m_parameters.percussiveThreshold) {
                c = Classification::Percussive;
            } else {
                c = Classification::Residual;
            }
            classification[i] = c;
        }
    }

protected:
    Parameters m_parameters;
    std::vector<std::shared_ptr<MovingMedian<float>>> m_hFilters;
    std::unique_ptr<MovingMedian<float>> m_vFilter;
    float *m_hf;
    float *m_vf;
    RingBuffer<float *> m_vfQueue;
};

}

#endif
