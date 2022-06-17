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

#include "../common/mathmisc.h"

#include <sstream>
#include <functional>

namespace RubberBand
{

class GuidedPhaseAdvance
{
public:
    struct Parameters {
        int fftSize;
        double sampleRate;
        int channels;
        std::function<void(const std::string &)> logger;
        Parameters(int _fftSize, double _sampleRate, int _channels,
                   std::function<void(const std::string &)> _log) :
            fftSize(_fftSize), sampleRate(_sampleRate),
            channels(_channels), logger(_log) { }
    };
    
    GuidedPhaseAdvance(Parameters parameters) :
        m_parameters(parameters),
        m_binCount(parameters.fftSize / 2 + 1),
        m_peakPicker(m_binCount),
        m_reported(false) {
        int ch = m_parameters.channels;
        m_currentPeaks = allocate_and_zero_channels<int>(ch, m_binCount);
        m_prevPeaks = allocate_and_zero_channels<int>(ch, m_binCount);
        m_greatestChannel = allocate_and_zero<int>(m_binCount);
        m_prevInPhase = allocate_and_zero_channels<double>(ch, m_binCount);
        m_prevOutPhase = allocate_and_zero_channels<double>(ch, m_binCount);
        m_unlocked = allocate_and_zero_channels<double>(ch, m_binCount);

        for (int c = 0; c < ch; ++c) {
            for (int i = 0; i < m_binCount; ++i) {
                m_prevPeaks[c][i] = i;
            }
        }
    }

    ~GuidedPhaseAdvance() {
        int ch = m_parameters.channels;
        deallocate_channels(m_currentPeaks, ch);
        deallocate_channels(m_prevPeaks, ch);
        deallocate(m_greatestChannel);
        deallocate_channels(m_prevInPhase, ch);
        deallocate_channels(m_prevOutPhase, ch);
        deallocate_channels(m_unlocked, ch);
    }

    void reset() {
        int ch = m_parameters.channels;
        v_zero_channels(m_prevPeaks, ch, m_binCount);
        v_zero_channels(m_prevInPhase, ch, m_binCount);
        v_zero_channels(m_prevOutPhase, ch, m_binCount);
    }
    
    void advance(double *const *outPhase,
                 const double *const *mag,
                 const double *const *phase,
                 const double *const *prevMag,
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

        int lowest = configuration.fftBandLimits[myFftBand].b0min;
        int highest = configuration.fftBandLimits[myFftBand].b1max;
        
        if (!m_reported) {
            std::ostringstream ostr;
            ostr << "PhaseAdvance: fftSize = " << m_parameters.fftSize
                 << ": bins = " << bs << ", channels = " << channels
                 << ", inhop = "<< inhop << ", outhop = " << outhop
                 << ", ratio = " << ratio << std::endl;
            ostr << "PhaseAdvance: lowest possible bin = " << lowest
                 << " (" << configuration.fftBandLimits[myFftBand].f0min
                 << "Hz), highest = " << highest
                 << " (" << configuration.fftBandLimits[myFftBand].f1max
                 << "Hz)" << std::endl;
            m_parameters.logger(ostr.str());
            m_reported = true;
        }
        
        for (int c = 0; c < channels; ++c) {
            for (int i = lowest; i <= highest; ++i) {
                m_currentPeaks[c][i] = i;
            }
            for (const auto &band : guidance[c]->phaseLockBands) {
                int startBin = binForFrequency
                    (band.f0, m_parameters.fftSize, m_parameters.sampleRate);
                int endBin = binForFrequency
                    (band.f1, m_parameters.fftSize, m_parameters.sampleRate);
                if (startBin > highest || endBin < lowest) continue;
                int count = endBin - startBin + 1;
                m_peakPicker.findNearestAndNextPeaks(mag[c],
                                                     startBin, count,
                                                     band.p, m_currentPeaks[c],
                                                     nullptr);
            }
            m_peakPicker.findNearestAndNextPeaks(prevMag[c],
                                                 lowest, highest - lowest + 1,
                                                 1, m_prevPeaks[c],
                                                 nullptr);
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
                double f = frequencyForBin
                    (i, m_parameters.fftSize, m_parameters.sampleRate);
                while (f > g->phaseLockBands[phaseLockBand].f1) {
                    ++phaseLockBand;
                }
                double ph = 0.0;
                if (inRange(f, g->phaseReset) || inRange(f, g->kick)) {
                    ph = phase[c][i];
                } else if (inhop == outhop) {
                    ph = m_unlocked[c][i];
                } else if (inRange (f, g->highUnlocked)) {
                    ph = m_unlocked[c][i];
                } else {
                    int peak = m_currentPeaks[c][i];
                    int prevPeak = m_prevPeaks[c][peak];
                    int peakCh = c;
                    if (inRange(f, g->channelLock)) {
                        int other = m_greatestChannel[i];
                        if (other != c &&
                            inRange(f, guidance[other]->channelLock)) {
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
                outPhase[c][i] = princarg(ph);
            }
        }
                
        for (int c = 0; c < channels; ++c) {
            for (int i = lowest; i <= highest; ++i) {
                m_prevInPhase[c][i] = phase[c][i];
            }
            for (int i = lowest; i <= highest; ++i) {
                m_prevOutPhase[c][i] = outPhase[c][i];
            }
        }
    }

protected:
    Parameters m_parameters;
    int m_binCount;
    Peak<double> m_peakPicker;
    int **m_currentPeaks;
    int **m_prevPeaks;
    int *m_greatestChannel;
    double **m_prevInPhase;
    double **m_prevOutPhase;
    double **m_unlocked;
    bool m_reported;

    bool inRange(double f, const Guide::Range &r) {
        return r.present && f >= r.f0 && f < r.f1;
    }

    GuidedPhaseAdvance(const GuidedPhaseAdvance &) =delete;
    GuidedPhaseAdvance &operator=(const GuidedPhaseAdvance &) =delete;
};

}

#endif
