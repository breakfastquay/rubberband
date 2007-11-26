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

    /**
     * Processing options for the timestretcher.  The preferred
     * options should normally be set in the constructor, as a bitwise
     * OR of the option flags.  The default value (DefaultOptions) is
     * intended to give good results in most situations.
     *
     * 1. Flags prefixed OptionProcess determine how the timestretcher
     * will be invoked.  These options may not be changed after
     * construction.
     * 
     *   OptionProcessOffline - Run the stretcher in offline mode.  In
     *   this mode the input data needs to be provided twice, once to
     *   study(), which calculates a stretch profile for the audio,
     *   and once to process(), which stretches it.
     *
     *   OptionProcessRealTime - Run the stretcher in real-time mode.
     *   In this mode only process() should be called, and the
     *   stretcher adjusts dynamically in response to the input audio.
     * 
     * The Process setting is likely to depend on your architecture:
     * non-real-time operation on seekable files: Offline; real-time
     * or streaming operation: RealTime.
     *
     * 2. Flags prefixed OptionStretch control the profile used for
     * variable timestretching.  Rubber Band always adjusts the
     * stretch profile to minimise stretching of busy broadband
     * transient sounds, but the degree to which it does so is
     * adjustable.  These options may not be changed after
     * construction.
     *
     *   OptionStretchElastic - Only meaningful in offline mode, and
     *   the default in that mode.  The audio will be stretched at a
     *   variable rate, aimed at preserving the quality of transient
     *   sounds as much as possible.  The timings of low activity
     *   regions between transients may be less exact than when the
     *   precise flag is set.
     * 
     *   OptionStretchPrecise - Although still using a variable
     *   stretch rate, the audio will be stretched so as to maintain
     *   as close as possible to a linear stretch ratio throughout.
     *   Timing may be better than when using OptionStretchElastic, at
     *   slight cost to the sound quality of transients.  This setting
     *   is always used when running in real-time mode.
     *
     * 3. Flags prefixed OptionTransients control the component
     * frequency phase-reset mechanism that may be used at transient
     * points to provide clarity and realism to percussion and other
     * significant transient sounds.  These options may be changed
     * after construction when running in real-time mode, but not when
     * running in offline mode.
     * 
     *   OptionTransientsCrisp - Reset component phases at the peak of
     *   each transient (the start of a significant note or percussive
     *   event).  This, the default setting, usually results in a
     *   clear-sounding output; but it is not always consistent, and
     *   may cause interruptions in stable sounds present at the same
     *   time as transient events.
     *
     *   OptionTransientsMixed - Reset component phases at the peak of
     *   each transient, outside a frequency range typical of musical
     *   fundamental frequencies.  The results may be more regular for
     *   mixed stable and percussive notes than OptionTransientsCrisp,
     *   but with a "phasier" sound.  The balance may sound very good
     *   for certain types of music and fairly bad for others.
     *
     *   OptionTransientsSmooth - Do not reset component phases at any
     *   point.  The results will be smoother and more regular but may
     *   be less clear than with either of the other transients flags.
     *
     * 4. Flags prefixed OptionPhase control the adjustment of
     * component frequency phases from one analysis window to the next
     * during non-transient segments.  These options may be changed at
     * any time.
     *
     *   OptionPhaseAdaptive - Lock the adjustments of phase for
     *   frequencies close to peak frequencies to those of the peak,
     *   but reduce the degree of locking as the stretch ratio gets
     *   longer.  This, the default setting, should give a good
     *   balance between clarity and smoothness in most situations.
     *
     *   OptionPhasePeakLocked - Lock the adjustments of phase for
     *   frequencies close to peak frequencies to those of the peak.
     *   This should give a clear result in situations with relatively
     *   low stretch ratios, but a relatively metallic sound at longer
     *   stretches.
     *
     *   OptionPhaseIndependent - Do not lock phase adjustments to
     *   peak frequencies.  This usually results in a softer, phasier
     *   sound.
     *
     * 5. Options prefixed OptionThreading control the threading model
     * of the stretcher.  These options may not be changed after
     * construction.
     *
     *   OptionThreadingAuto - Permit the stretcher to determine its
     *   own threading model.  Usually this means using one processing
     *   thread per audio channel in offline mode if the stretcher is
     *   able to determine that more than one CPU is available, and
     *   one thread only in realtime mode.
     *
     *   OptionThreadingNever - Never use more than one thread.
     *  
     *   OptionThreadingAlways - Use multiple threads in any situation
     *   where OptionThreadingAuto would do so, except omit the check
     *   for multiple CPUs and instead assume it to be true.
     *
     * 6. Options prefixed OptionWindow control the window size for
     * FFT processing.  The window size actually used will depend on
     * many factors, but it can be influenced.  These options may not
     * be changed after construction.
     *
     *   OptionWindowStandard - Use the default window size.  The
     *   actual size will vary depending on other parameters.  This
     *   option is expected to produce better results than the other
     *   window options in most situations.
     *
     *   OptionWindowShort - Use a shorter window.  This may result in
     *   crisper sound for audio that depends strongly on its timing
     *   qualities.
     *
     *   OptionWindowLong - Use a longer window.  This is likely to
     *   result in a smoother sound at the expense of clarity and
     *   timing.
     */
    typedef int Options;
    
    static const int OptionProcessOffline   = 0x00000000;
    static const int OptionProcessRealTime  = 0x00000001;
    
    static const int OptionStretchElastic   = 0x00000000;
    static const int OptionStretchPrecise   = 0x00000010;
    
    static const int OptionTransientsCrisp  = 0x00000000;
    static const int OptionTransientsMixed  = 0x00000100;
    static const int OptionTransientsSmooth = 0x00000200;

    static const int OptionPhaseAdaptive    = 0x00000000;
    static const int OptionPhasePeakLocked  = 0x00001000;
    static const int OptionPhaseIndependent = 0x00002000;
    
    static const int OptionThreadingAuto    = 0x00000000;
    static const int OptionThreadingNever   = 0x00010000;
    static const int OptionThreadingAlways  = 0x00020000;

    static const int OptionWindowStandard   = 0x00000000;
    static const int OptionWindowShort      = 0x00100000;
    static const int OptionWindowLong       = 0x00200000;

    static const int DefaultOptions         = 0x00000000;
    static const int PercussiveOptions      = OptionWindowShort | \
                                              OptionPhaseIndependent;

    /**
     * Construct a time-and-pitch-scaling object to run at the given
     * sample rate, with the given number of channels.  Processing
     * options and the time and pitch scaling ratios may be provided.
     * The time and pitch ratios may be changed after construction,
     * but most of the options may not.  See the option documentation
     * above for more details.
     */
    RubberBandStretcher(size_t sampleRate,
                        size_t channels,
                        Options options = DefaultOptions,
                        double initialTimeRatio = 1.0,
                        double initialPitchScale = 1.0);
    virtual ~RubberBandStretcher();
    
    virtual void reset();
    virtual void setTimeRatio(double ratio);
    virtual void setPitchScale(double scale);

    virtual double getTimeRatio() const;
    virtual double getPitchScale() const;

    virtual size_t getLatency() const;

    virtual void setTransientsOption(Options options);
    virtual void setPhaseOption(Options options);

    virtual void setExpectedInputDuration(size_t samples);
    virtual void setMaxProcessSize(size_t samples);
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
    virtual std::vector<float> getPhaseResetCurve() const; //!!! document particular meaning in RT mode
    virtual std::vector<int> getExactTimePoints() const; //!!! meaningless in RT mode

    virtual size_t getChannelCount() const;
    
    virtual void calculateStretch();

    virtual void setDebugLevel(int level);

    static void setDefaultDebugLevel(int level);

protected:
    class Impl;
    Impl *m_d;
};

}

#endif
