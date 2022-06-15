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

#ifndef RUBBERBAND_GUIDE_H
#define RUBBERBAND_GUIDE_H

#include <functional>
#include <sstream>

namespace RubberBand
{

class Guide
{
public:
    struct FftBand {
        int fftSize;
        double f0;
        double f1;
        FftBand(int _s, double _f0, double _f1) :
            fftSize(_s), f0(_f0), f1(_f1) { }
        FftBand() :
            fftSize(0), f0(0.f), f1(0.f) { }
    };

    struct PhaseLockBand {
        int p;
        double beta;
        double f0;
        double f1;
        PhaseLockBand(int _p, double _beta, double _f0, double _f1) :
            p(_p), beta(_beta), f0(_f0), f1(_f1) { }
        PhaseLockBand() :
            p(0), beta(1.0), f0(0.f), f1(0.f) { }
    };

    struct Range {
        bool present;
        double f0;
        double f1;
        Range(bool _present, double _f0, double _f1) :
            present(_present), f0(_f0), f1(_f1) { }
        Range() :
            present(false), f0(0.f), f1(0.f) { }
    };
        
    struct Guidance {
        FftBand fftBands[3];
        PhaseLockBand phaseLockBands[4];
        Range kick;
        Range preKick;
        Range highUnlocked;
        Range phaseReset;
        Range channelLock;
    };

    struct BandLimits {
        int fftSize;
        double f0min;
        double f1max;
        int b0min;
        int b1max;
        BandLimits(int _fftSize, double _rate, double _f0min, double _f1max) :
            fftSize(_fftSize), f0min(_f0min), f1max(_f1max),
            b0min(int(floor(f0min * fftSize / _rate))),
            b1max(int(ceil(f1max * fftSize / _rate))) { }
        BandLimits() :
            fftSize(0), f0min(0.f), f1max(0.f), b0min(0), b1max(0) { }
    };

    struct Configuration {
        int longestFftSize;
        int shortestFftSize;
        int classificationFftSize;
        BandLimits fftBandLimits[3];
        Configuration(int _longestFftSize, int _shortestFftSize,
                      int _classificationFftSize) :
            longestFftSize(_longestFftSize),
            shortestFftSize(_shortestFftSize),
            classificationFftSize(_classificationFftSize) { }
    };
    
    struct Parameters {
        double sampleRate;
        std::function<void(const std::string &)> logger;
        Parameters(double _sampleRate,
                   std::function<void(const std::string &)> _log) :
            sampleRate(_sampleRate),
            logger(_log) { }
    };

    Guide(Parameters parameters) :
        m_parameters(parameters),
        m_configuration(roundUp(int(ceil(parameters.sampleRate / 16.0))),
                        roundUp(int(ceil(parameters.sampleRate / 64.0))),
                        roundUp(int(ceil(parameters.sampleRate / 32.0)))),
        m_minLower(500.0), m_minHigher(4000.0),
        m_defaultLower(700.0), m_defaultHigher(4800.0),
        m_maxLower(1100.0), m_maxHigher(7000.0)
    {
        double rate = m_parameters.sampleRate;

        int bandFftSize = roundUp(int(ceil(rate/16.0)));
        m_configuration.fftBandLimits[0] =
            BandLimits(bandFftSize, rate, 0.0, m_maxLower);

        // This is the classification and fallback FFT: we need the
        // full range for it
        bandFftSize = roundUp(int(ceil(rate/32.0)));
        m_configuration.fftBandLimits[1] =
            BandLimits(bandFftSize, rate, 0.0, rate / 2.0);
        
        bandFftSize = roundUp(int(ceil(rate/64.0)));
        m_configuration.fftBandLimits[2] =
            BandLimits(bandFftSize, rate, m_minHigher, rate/2.0);
    }

    const Configuration &getConfiguration() const {
        return m_configuration;
    }
    
    void updateGuidance(double ratio,
                        int outhop,
                        const double *const magnitudes,
                        const double *const prevMagnitudes,
                        const double *const nextMagnitudes,
                        const BinSegmenter::Segmentation &segmentation,
                        const BinSegmenter::Segmentation &prevSegmentation,
                        const BinSegmenter::Segmentation &nextSegmentation,
                        bool specialCaseUnity,
                        Guidance &guidance) const {

        guidance.kick.present = false;
        guidance.preKick.present = false;
        guidance.highUnlocked.present = false;
        guidance.phaseReset.present = false;

        double nyquist = m_parameters.sampleRate / 2.0;

        guidance.fftBands[0].fftSize = roundUp(int(ceil(nyquist/8.0)));
        guidance.fftBands[1].fftSize = roundUp(int(ceil(nyquist/16.0)));
        guidance.fftBands[2].fftSize = roundUp(int(ceil(nyquist/32.0)));

        if (specialCaseUnity && (fabs(ratio - 1.0) < 1.0e-6)) {
            guidance.fftBands[0].f0 = 0.0;
            guidance.fftBands[0].f1 = 0.0;
            guidance.fftBands[1].f0 = 0.0;
            guidance.fftBands[1].f1 = m_minHigher;
            guidance.fftBands[2].f0 = m_minHigher;
            guidance.fftBands[2].f1 = nyquist;
            for (int i = 0; i < int(sizeof(guidance.phaseLockBands) /
                                    sizeof(guidance.phaseLockBands[0])); ++i) {
                guidance.phaseLockBands[i].p = 0;
                guidance.phaseLockBands[i].beta = 1.0;
                guidance.phaseLockBands[i].f0 = nyquist;
                guidance.phaseLockBands[i].f1 = nyquist;
            }
            guidance.phaseLockBands[0].f0 = 0.0;
            guidance.phaseLockBands[0].f1 = nyquist;
            guidance.channelLock.present = false;
            return;
        }

        guidance.channelLock.present = true;
        guidance.channelLock.f0 = 0.0;
        guidance.channelLock.f1 = 600.0;

        bool kick =
            (segmentation.percussiveBelow > 40.0) &&
            (prevSegmentation.percussiveBelow < 40.0) &&
            checkPotentialKick(magnitudes, prevMagnitudes);

        bool futureKick = !kick &&
            (nextSegmentation.percussiveBelow > 40.0) &&
            (segmentation.percussiveBelow < 40.0) &&
            checkPotentialKick(nextMagnitudes, magnitudes);
/*
        std::cout << "d:"
                  << prevSegmentation.percussiveBelow << ","
                  << segmentation.percussiveBelow << ","
                  << nextSegmentation.percussiveBelow << ","
                  << checkPotentialKick(magnitudes, prevMagnitudes) << ","
                  << checkPotentialKick(nextMagnitudes, magnitudes) << ","
                  << (kick ? "K" : "N") << ","
                  << (futureKick ? "F" : "N") << std::endl;
*/        
        if (kick) {
            guidance.kick.present = true;
            guidance.kick.f0 = 0.0;
            guidance.kick.f1 = segmentation.percussiveBelow;
        } else if (futureKick) {
            guidance.preKick.present = true;
            guidance.preKick.f0 = 0.0;
            guidance.preKick.f1 = nextSegmentation.percussiveBelow;
        }
        
        if (segmentation.residualAbove > segmentation.percussiveAbove) {
            guidance.highUnlocked.present = true;
            guidance.highUnlocked.f0 = segmentation.percussiveAbove;
            guidance.highUnlocked.f1 = segmentation.residualAbove;
        }

        double bigGap = 4000.0;
        if (segmentation.residualAbove >
            segmentation.percussiveAbove + bigGap &&
            prevSegmentation.residualAbove <
            prevSegmentation.percussiveAbove + bigGap) {
            guidance.phaseReset.present = true;
            guidance.phaseReset.f0 = std::min(segmentation.percussiveAbove,
                                              nextSegmentation.percussiveAbove);
            guidance.phaseReset.f1 = std::max(segmentation.residualAbove,
                                              nextSegmentation.residualAbove);
            if (guidance.phaseReset.f0 < 200.0) {
                guidance.phaseReset.f0 = 0.0;
            }
        }

        double prevLower = guidance.fftBands[0].f1;
        double lower = descendToValley(prevLower, magnitudes);
        if (lower > m_maxLower || lower < m_minLower) {
            lower = m_defaultLower;
        }
        
        double prevHigher = guidance.fftBands[1].f1;
        double higher = descendToValley(prevHigher, magnitudes);
        if (higher > m_maxHigher || higher < m_minHigher) {
            higher = m_defaultHigher;
        }

        guidance.fftBands[0].f0 = 0.0;
        guidance.fftBands[0].f1 = lower;

//        std::cout << "x:" << lower << std::endl;
        
        guidance.fftBands[1].f0 = lower;
        guidance.fftBands[1].f1 = higher;
        
        guidance.fftBands[2].f0 = higher;
        guidance.fftBands[2].f1 = nyquist;

        if (outhop > 256) {
            guidance.fftBands[1].f1 = nyquist;
            guidance.fftBands[2].f0 = nyquist;
        }
        
        double mid = std::max(lower, 1600.0);

        guidance.phaseLockBands[0].p = 1;
        guidance.phaseLockBands[0].beta = betaFor(300.0, ratio);
        guidance.phaseLockBands[0].f0 = 0.0;
        guidance.phaseLockBands[0].f1 = lower;
        
        guidance.phaseLockBands[1].p = 2;
        guidance.phaseLockBands[1].beta = betaFor(1600.0, ratio);
        guidance.phaseLockBands[1].f0 = lower;
        guidance.phaseLockBands[1].f1 = mid;
        
        guidance.phaseLockBands[2].p = 3;
        guidance.phaseLockBands[2].beta = betaFor(5000.0, ratio);
        guidance.phaseLockBands[2].f0 = mid;
        guidance.phaseLockBands[2].f1 = higher;
        
        guidance.phaseLockBands[3].p = 4;
        guidance.phaseLockBands[3].beta = betaFor(10000.0, ratio);
        guidance.phaseLockBands[3].f0 = higher;
        guidance.phaseLockBands[3].f1 = nyquist;
        
        if (outhop > 256) {
            guidance.phaseLockBands[3].p = 3;
        }

        if (ratio > 2.0) {
            
            // For very long stretches, diffuse is better than
            // metallic - gradually unlock the higher frequencies and
            // reduce the channel lock
            
            double channelLimit = guidance.channelLock.f1;
            channelLimit = channelLimit - (ratio - 2.0) * 150.0;
            if (channelLimit < 100.0) channelLimit = 100.0;
            guidance.channelLock.f1 = channelLimit;
            
            double unlockedAbove = 12000.0 - (ratio - 2.0) * 400.0;
            if (unlockedAbove < channelLimit) unlockedAbove = channelLimit;
            if (guidance.highUnlocked.present) {
                guidance.highUnlocked.f0 = std::min(guidance.highUnlocked.f0,
                                                    unlockedAbove);
            } else {
                guidance.highUnlocked.f0 = unlockedAbove;
            }
            guidance.highUnlocked.f1 = nyquist;
            guidance.highUnlocked.present = true;
        }

        /*
        std::ostringstream str;
        str << "Guidance: FFT bands: ["
            << guidance.fftBands[0].fftSize << " from "
            << guidance.fftBands[0].f0 << " to " << guidance.fftBands[0].f1
            << ", "
            << guidance.fftBands[1].fftSize << " from "
            << guidance.fftBands[1].f0 << " to " << guidance.fftBands[1].f1
            << ", "
            << guidance.fftBands[2].fftSize << " from "
            << guidance.fftBands[2].f0 << " to " << guidance.fftBands[2].f1
            << "]; phase reset range: ["
            << guidance.phaseReset.present << " from "
            << guidance.phaseReset.f0 << " to " << guidance.phaseReset.f1
            << "]" << std::endl;
        m_parameters.logger(str.str());
        */
    }

protected:
    Parameters m_parameters;
    Configuration m_configuration;

    double m_minLower;
    double m_minHigher;
    double m_defaultLower;
    double m_defaultHigher;
    double m_maxLower;
    double m_maxHigher;
    
    // near-dupe with R2 RubberBandStretcher::Impl
    int roundUp(int value) const {
        if (value < 1) return 1;
        if (!(value & (value - 1))) return value;
        int bits = 0;
        while (value) { ++bits; value >>= 1; }
        value = 1 << bits;
        return value;
    }
    
    bool checkPotentialKick(const double *const magnitudes,
                            const double *const prevMagnitudes) const {
        int b = binForFrequency(200.0, m_configuration.classificationFftSize,
                                m_parameters.sampleRate);
        double here = 0.0, there = 0.0;
        for (int i = 1; i <= b; ++i) {
            here += magnitudes[i];
        }
        for (int i = 1; i <= b; ++i) {
            there += prevMagnitudes[i];
        }
        return (here > 10.e-3 && here > there * 1.4);
    }

    double descendToValley(double f, const double *const magnitudes) const {
        int b = binForFrequency(f, m_configuration.classificationFftSize,
                                m_parameters.sampleRate);
        for (int i = 0; i < 3; ++i) {
            if (magnitudes[b+1] < magnitudes[b]) {
                ++b;
            } else if (magnitudes[b-1] < magnitudes[b]) {
                --b;
            } else {
                break;
            }
        }
        double sf = frequencyForBin(b, m_configuration.classificationFftSize,
                                    m_parameters.sampleRate);
        return sf;
    }

    double betaFor(double f, double ratio) const {
        double b = (2.0 + ratio) / 3.0;
        double limit = 10000.0;
        if (f > limit) {
            return b;
        } else {
            return 1.0 + f * (b - 1.0) / limit;
        }
    }
};

}

#endif
