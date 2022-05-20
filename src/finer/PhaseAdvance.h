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

#ifndef RUBBERBAND_PHASE_ADVANCE_H
#define RUBBERBAND_PHASE_ADVANCE_H

#include "Guide.h"

namespace RubberBand
{

class GuidedPhaseAdvance
{
public:
    struct Parameters {
        int fftSize;
        double sampleRate;
        int channels;
        Parameters(int _fftSize, double _sampleRate, int _channels) :
            fftSize(_fftSize), sampleRate(_sampleRate), channels(_channels) { }
    };
    
    GuidedPhaseAdvance(Parameters parameters) :
        m_parameters(parameters),
        m_blockSize(parameters.fftSize / 2 + 1),
        m_peakPicker(m_blockSize) {
        size_t ch = m_parameters.channels;
        m_currentPeaks = allocate_and_zero_channels<int>(ch, m_blockSize);
        m_prevPeaks = allocate_and_zero_channels<int>(ch, m_blockSize);
        m_greatestChannel = allocate_and_zero<int>(m_blockSize);
        m_prevInPhase = allocate_and_zero_channels<float>(ch, m_blockSize);
        m_prevOutPhase = allocate_and_zero_channels<double>(ch, m_blockSize);
        m_unlocked = allocate_and_zero_channels<double>(ch, m_blockSize);
    }

    ~GuidedPhaseAdvance() {
        size_t ch = m_parameters.channels;
        deallocate_channels(m_currentPeaks, ch);
        deallocate_channels(m_prevPeaks, ch);
        deallocate(m_greatestChannel);
        deallocate_channels(m_prevInPhase, ch);
        deallocate_channels(m_prevOutPhase, ch);
        deallocate_channels(m_unlocked, ch);
    }
    
    void advance(double *const *outPhase,
                 const float *const *mag,
                 const float *const *phase,
                 const Guide::Configuration &configuration,
                 const Guide::Guidance *const *guidance,
                 int inhop,
                 int outhop) {

        int myFftBand = 0;
        int i = 0;
        for (const auto &fband : guidance[0]->fftBands) {
            if (fband.fftSize == m_parameters.fftSize) {
                myFftBand = i;
                break;
            }
            ++i;
        }

        int bs = m_parameters.fftSize / 2 + 1;
        int channels = m_parameters.channels;
        double ratio = double(outhop) / double(inhop);

        int lowest = binForFrequency
            (configuration.fftBandLimits[myFftBand].f0min);
        int highest = binForFrequency
            (configuration.fftBandLimits[myFftBand].f1max);
        
        for (int c = 0; c < channels; ++c) {
            for (int i = lowest; i <= highest; ++i) {
                m_currentPeaks[c][i] = i;
            }
            for (const auto &band : guidance[c]->phaseLockBands) {
                int startBin = binForFrequency(band.f0);
                int endBin = binForFrequency(band.f1);
                if (startBin > highest || endBin < lowest) continue;
                int count = endBin - startBin;
                m_peakPicker.findNearestAndNextPeaks(mag[c], startBin, count,
                                                     band.p, m_currentPeaks[c]);
            }
        }

        if (channels > 1) {
            for (int i = lowest; i <= highest; ++i) {
                int gc = 0;
                float gmag = mag[0][i];
                for (int c = 1; c < channels; ++c) {
                    if (mag[c][i] > gmag) {
                        gmag = mag[c][i];
                        gc = c;
                    }
                }
                m_greatestChannel[i] = gc;
            }
        } else {
            v_zero(m_greatestChannel, bs);
        }

        double omegaFactor = 2.0 * M_PI * double(inhop) /
            double(m_parameters.fftSize);
        for (int c = 0; c < channels; ++c) {
            for (int i = lowest; i <= highest; ++i) {
                double omega = omegaFactor * double(i);
                double expected = m_prevInPhase[c][i] + omega;
                double error = princarg(phase[c][i] - expected);
                double advance = ratio * (omega + error);
                m_unlocked[c][i] = m_prevOutPhase[c][i] + advance;
            }
        }

        for (int c = 0; c < channels; ++c) {
            const Guide::Guidance *g = guidance[c];
            int phaseLockBand = 0;
            for (int i = lowest; i <= highest; ++i) {
                double f = frequencyForBin(i);
                while (f > g->phaseLockBands[phaseLockBand].f1) {
                    ++phaseLockBand;
                }
                double ph = 0.0;
                if (inRange(f, g->phaseReset) || inRange(f, g->kick)) {
                    ph = phase[c][i];
                } else if (inRange (f, g->highPercussive)) {
                    ph = m_unlocked[c][i];
                } else {
                    int peak = m_currentPeaks[c][i];
                    int prevPeak = m_prevPeaks[c][peak];
                    int peakCh = c;
                    if (inRange (f, g->channelLock)) {
                        int other = m_greatestChannel[i];
                        if (other != c) {
                            int otherPeak = m_currentPeaks[other][i];
                            int otherPrevPeak = m_prevPeaks[other][otherPeak];
                            if (otherPrevPeak == prevPeak) {
                                peakCh = other;
                            }
                        }
                    }
                    double peakAdvance =
                        m_unlocked[peakCh][peak] - m_prevOutPhase[peakCh][peak];
                    double peakNew =
                        m_prevOutPhase[peakCh][prevPeak] + peakAdvance;
                    double diff =
                        double(phase[c][i]) - double(phase[peakCh][peak]);
                    double beta =
                        double(g->phaseLockBands[phaseLockBand].beta);
                    ph = peakNew + beta * diff;
                }
                outPhase[c][i] = ph;
            }
        }
                
        for (int c = 0; c < channels; ++c) {
            for (int i = lowest; i <= highest; ++i) {
                m_prevInPhase[c][i] = phase[c][i];
            }
        }
        for (int c = 0; c < channels; ++c) {
            for (int i = lowest; i <= highest; ++i) {
                m_prevOutPhase[c][i] = outPhase[c][i];
            }
        }

        int **tmp = m_prevPeaks;
        m_prevPeaks = m_currentPeaks;
        m_currentPeaks = m_prevPeaks;
    }

protected:
    Parameters m_parameters;
    int m_blockSize;
    Peak<float> m_peakPicker;
    int **m_currentPeaks;
    int **m_prevPeaks;
    int *m_greatestChannel;
    float **m_prevInPhase;
    double **m_prevOutPhase;
    double **m_unlocked;

    int binForFrequency(double f) const {
        return int(round(f * double(m_parameters.fftSize) /
                         m_parameters.sampleRate));
    }
    double frequencyForBin(int b) const {
        return (double(b) * m_parameters.sampleRate)
            / double(m_parameters.fftSize);
    }
    bool inRange(double f, const Guide::Range &r) {
        return r.present && f >= r.f0 && f < r.f1;
    }
};

}

#endif
