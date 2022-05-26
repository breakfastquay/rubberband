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
        PhaseLockBand phaseLockBands[5];
        Range kick;
        Range lowPercussive;
        Range highPercussive;
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
        m_minLower(350.0), m_minHigher(2400.0),
        m_defaultLower(700.0), m_defaultHigher(4800.0),
        m_maxLower(1100.0), m_maxHigher(7000.0)
    {
        double rate = m_parameters.sampleRate;

        int bandFftSize = roundUp(int(ceil(rate/16.0)));
        m_configuration.fftBandLimits[0] =
            BandLimits(bandFftSize, rate, 0.0, m_maxLower);

        bandFftSize = roundUp(int(ceil(rate/32.0)));
        m_configuration.fftBandLimits[1] =
            BandLimits(bandFftSize, rate, m_minLower, m_maxHigher);
        
        bandFftSize = roundUp(int(ceil(rate/64.0)));
        m_configuration.fftBandLimits[2] =
            BandLimits(bandFftSize, rate, m_minHigher, rate/2.0);
    }

    const Configuration &getConfiguration() const {
        return m_configuration;
    }
    
    void calculate(double ratio,
                   const double *const magnitudes,
                   const int *const troughs,
                   const double *const prevMagnitudes,
                   const BinSegmenter::Segmentation &segmentation,
                   const BinSegmenter::Segmentation &prevSegmentation,
                   const BinSegmenter::Segmentation &nextSegmentation,
                   Guidance &guidance) const {

        bool potentialKick = checkPotentialKick(magnitudes, prevMagnitudes);

        guidance.kick.present = false;
        guidance.lowPercussive.present = false;
        guidance.highPercussive.present = false;
        guidance.phaseReset.present = false;

        guidance.channelLock.present = true;
        guidance.channelLock.f0 = 0.0;
        guidance.channelLock.f1 = 600.0;

        if (segmentation.percussiveBelow > 40.0) {
            guidance.lowPercussive.present = true;
            guidance.lowPercussive.f0 = 0.0;
            guidance.lowPercussive.f1 = segmentation.percussiveBelow;
        }

        if (potentialKick && prevSegmentation.percussiveBelow < 40.0) {
            guidance.kick = guidance.lowPercussive;
        }
        
        if (segmentation.residualAbove > segmentation.percussiveAbove) {
            guidance.highPercussive.present = true;
            guidance.highPercussive.f0 = segmentation.percussiveAbove;
            guidance.highPercussive.f1 = segmentation.residualAbove;
        }

        double bigGap = 4000.0;
        if (ratio > 1.0 &&
            segmentation.residualAbove >
            segmentation.percussiveAbove + bigGap &&
            prevSegmentation.residualAbove <
            prevSegmentation.percussiveAbove + bigGap) {
            guidance.phaseReset.present = true;
            guidance.phaseReset.f0 = std::min(segmentation.percussiveAbove,
                                              nextSegmentation.percussiveAbove);
            if (guidance.phaseReset.f0 < 200.0) {
                guidance.phaseReset.f0 = 0.0;
            }
            guidance.phaseReset.f1 = std::max(segmentation.residualAbove,
                                              nextSegmentation.residualAbove);
        }

        double higher = snapToTrough(m_defaultHigher, troughs);
        if (higher > m_maxHigher) higher = m_maxHigher;

        double lower = snapToTrough(m_defaultLower, troughs);
        if (lower > m_maxLower) lower = m_maxLower;

        double nyquist = m_parameters.sampleRate / 2.0;

        guidance.fftBands[0].fftSize = roundUp(int(ceil(nyquist/8.0)));
        guidance.fftBands[0].f0 = 0.0;
        guidance.fftBands[0].f1 = lower;
        
        guidance.fftBands[1].fftSize = roundUp(int(ceil(nyquist/16.0)));
        guidance.fftBands[1].f0 = lower;
        guidance.fftBands[1].f1 = higher;
        
        guidance.fftBands[2].fftSize = roundUp(int(ceil(nyquist/32.0)));
        guidance.fftBands[2].f0 = higher;
        guidance.fftBands[2].f1 = nyquist;
        
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
    
    int binForFrequency(double f) const {
        return int(round(f * double(m_configuration.classificationFftSize) /
                         m_parameters.sampleRate));
    }
    double frequencyForBin(int b) const {
        return (double(b) * m_parameters.sampleRate)
            / double(m_configuration.classificationFftSize);
    }

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
        int b = binForFrequency(200.0);
        double here = 0.0, there = 0.0;
        for (int i = 1; i <= b; ++i) {
            here += magnitudes[i];
        }
        for (int i = 1; i <= b; ++i) {
            there += prevMagnitudes[i];
        }
        return (here > 10.e-3 && here > there * 1.4);
    }

    double snapToTrough(double f, const int *const troughs) const {
        return frequencyForBin(troughs[binForFrequency(f)]);
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
