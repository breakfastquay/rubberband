/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band
    An audio time-stretching and pitch-shifting library.
    Copyright 2007 Chris Cannam.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#ifndef _RUBBERBANDSTRETCHER_H_
#define _RUBBERBANDSTRETCHER_H_

#include "TimeStretcher.h"

#include <vector>

namespace RubberBand
{

class RubberBandStretcher : public TimeStretcher
{
public:

    typedef int Options;
    
    static const int OptionProcessOffline   = 0x00000000;
    static const int OptionProcessRealTime  = 0x00000001;
    
    static const int OptionStretchElastic   = 0x00000000;
    static const int OptionStretchPrecise   = 0x00000010;
    
    static const int OptionTransientsCrisp  = 0x00000000;
    static const int OptionTransientsSmooth = 0x00000100;

    static const int OptionPhasePeakLocked  = 0x00000000;
    static const int OptionPhaseIndependent = 0x00001000;
    
    static const int OptionThreadingAuto    = 0x00000000;
    static const int OptionThreadingNone    = 0x00010000;

    static const int OptionWindowStandard   = 0x00000000;
    static const int OptionWindowShort      = 0x00100000;
    static const int OptionWindowLong       = 0x00200000;

    static const int DefaultOptions         = 0x00000000;
    static const int PercussiveOptions      = OptionWindowShort | \
                                              OptionPhaseIndependent;

    RubberBandStretcher(size_t sampleRate,
                        size_t channels,
                        Options options = DefaultOptions,
                        double initialTimeRatio = 1.0,
                        double initialPitchScale = 1.0);
    virtual ~RubberBandStretcher();
    
    virtual void reset();
    virtual void setTimeRatio(double ratio);
    virtual void setPitchScale(double scale); //!!!??? pitch ratio?

    virtual double getTimeRatio() const;
    virtual double getPitchScale() const;

    virtual size_t getLatency() const;

    virtual void setTransientsOption(Options options);
    virtual void setPhaseOption(Options options);

    virtual void setExpectedInputDuration(size_t samples);
    virtual void setMaxProcessBlockSize(size_t samples);
    virtual size_t getSamplesRequired() const;

    // if samples == 0, input may be null
    virtual void study(const float *const *input, size_t samples, bool final);
    virtual void process(const float *const *input, size_t samples, bool final);

    virtual int available() const; // returns -1 if all data processed and nothing further will be available
    virtual size_t retrieve(float *const *output, size_t samples) const;

    virtual float getFrequencyCutoff(int n) const;
    virtual void setFrequencyCutoff(int n, float f);
    
    //!!! ideally, this stuff wouldn't be here...

    virtual size_t getInputIncrement() const;
    virtual std::vector<int> getOutputIncrements() const; //!!! document particular meaning in RT mode
    virtual std::vector<float> getLockCurve() const; //!!! document particular meaning in RT mode
    virtual std::vector<int> getExactTimePoints() const; //!!! meaningless in RT mode

    virtual size_t getChannelCount() const;
    
    virtual void calculateStretch();

    virtual void setDebugLevel(int level);

protected:
    class Impl;
    Impl *m_d;
};

}

#endif
