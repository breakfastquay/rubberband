/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2008 Chris Cannam.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "StretcherImpl.h"

#include "rubberband.h"

struct RubberBandState_
{
    RubberBand::StretcherImpl *m_impl;
};

RubberBandState rubberband_new(unsigned int sampleRate,
                               unsigned int channels,
                               RubberBandOptions options, //!!! sort out RubberBand namespacing
                               double initialTimeRatio,
                               double initialPitchScale)
{
    RubberBandState_ *state = new RubberBandState_();
    state->m_impl = new RubberBand::StretcherImpl
        (sampleRate, channels, options,
         initialTimeRatio, initialPitchScale);
    return state;
}

void rubberband_delete(RubberBandState state)
{
    delete state->m_impl;
    delete state;
}

void rubberband_reset(RubberBandState state)
{
    state->m_impl->reset();
}

void rubberband_set_time_ratio(RubberBandState state, double ratio)
{
    state->m_impl->setTimeRatio(ratio);
}

void rubberband_set_pitch_scale(RubberBandState state, double scale)
{
    state->m_impl->setPitchScale(scale);
}

double rubberband_get_time_ratio(const RubberBandState state) 
{
    return state->m_impl->getTimeRatio();
}

double rubberband_get_pitch_scale(const RubberBandState state)
{
    return state->m_impl->getPitchScale();
}

unsigned int rubberband_get_latency(const RubberBandState state) 
{
    return state->m_impl->getLatency();
}

void rubberband_set_transients_option(RubberBandState state, RubberBandOptions options)
{
    state->m_impl->setTransientsOption(options);
}

void rubberband_set_phase_option(RubberBandState state, RubberBandOptions options)
{
    state->m_impl->setPhaseOption(options);
}

void rubberband_set_formant_option(RubberBandState state, RubberBandOptions options)
{
    state->m_impl->setFormantOption(options);
}

void rubberband_set_pitch_option(RubberBandState state, RubberBandOptions options)
{
    state->m_impl->setPitchOption(options);
}

void rubberband_set_expected_input_duration(RubberBandState state, unsigned int samples)
{
    state->m_impl->setExpectedInputDuration(samples);
}

unsigned int rubberband_get_samples_required(const RubberBandState state)
{
    return state->m_impl->getSamplesRequired();
}

void rubberband_set_max_process_size(RubberBandState state, unsigned int samples)
{
    state->m_impl->setMaxProcessSize(samples);
}

void rubberband_study(RubberBandState state, const float *const *input, unsigned int samples, int final)
{
    state->m_impl->study(input, samples, final != 0);
}

void rubberband_process(RubberBandState state, const float *const *input, unsigned int samples, int final)
{
    state->m_impl->process(input, samples, final != 0);
}

int rubberband_available(const RubberBandState state)
{
    return state->m_impl->available();
}

unsigned int rubberband_retrieve(const RubberBandState state, float *const *output, unsigned int samples)
{
    return state->m_impl->retrieve(output, samples);
}

unsigned int rubberband_get_channel_count(const RubberBandState state)
{
    return state->m_impl->getChannelCount();
}

void rubberband_calculate_stretch(RubberBandState state)
{
    state->m_impl->calculateStretch();
}

void rubberband_set_debug_level(RubberBandState state, int level)
{
    state->m_impl->setDebugLevel(level);
}

void rubberband_set_default_debug_level(int level)
{
    RubberBand::StretcherImpl::setDefaultDebugLevel(level);
}
