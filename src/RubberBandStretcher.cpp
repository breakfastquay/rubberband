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

#include "faster/StretcherImpl.h"
#include "finer/R3StretcherImpl.h"


namespace RubberBand {

//#define FASTER 1

RubberBandStretcher::RubberBandStretcher(size_t sampleRate,
                                         size_t channels,
                                         Options options,
                                         double initialTimeRatio,
                                         double initialPitchScale) :
    m_d
    (!(options & OptionEngineFiner) ?
     new Impl(sampleRate, channels, options,
              initialTimeRatio, initialPitchScale)
     : nullptr),
    m_r3d
    ((options & OptionEngineFiner) ?
     new R3StretcherImpl(R3StretcherImpl::Parameters
                         (sampleRate, channels, options),
                         initialTimeRatio, initialPitchScale)
     : nullptr)
{
}

RubberBandStretcher::~RubberBandStretcher()
{
    delete m_d;
    delete m_r3d;
}

void
RubberBandStretcher::reset()
{
    if (m_d) m_d->reset();
    else m_r3d->reset();
}

void
RubberBandStretcher::setTimeRatio(double ratio)
{
    if (m_d) m_d->setTimeRatio(ratio);
    else m_r3d->setTimeRatio(ratio);
}

void
RubberBandStretcher::setPitchScale(double scale)
{
    if (m_d) m_d->setPitchScale(scale);
    else m_r3d->setPitchScale(scale);
}

double
RubberBandStretcher::getTimeRatio() const
{
    if (m_d) return m_d->getTimeRatio();
    else return m_r3d->getTimeRatio();
}

double
RubberBandStretcher::getPitchScale() const
{
    if (m_d) return m_d->getPitchScale();
    else return m_r3d->getPitchScale();
}

size_t
RubberBandStretcher::getLatency() const
{
    if (m_d) return m_d->getLatency();
    else return m_r3d->getLatency();
}

void
RubberBandStretcher::setTransientsOption(Options options) 
{
    if (m_d) m_d->setTransientsOption(options);
}

void
RubberBandStretcher::setDetectorOption(Options options) 
{
    if (m_d) m_d->setDetectorOption(options);
}

void
RubberBandStretcher::setPhaseOption(Options options) 
{
    if (m_d) m_d->setPhaseOption(options);
}

void
RubberBandStretcher::setFormantOption(Options options)
{
    if (m_d) m_d->setFormantOption(options);
}

void
RubberBandStretcher::setPitchOption(Options options)
{
    if (m_d) m_d->setPitchOption(options);
}

void
RubberBandStretcher::setExpectedInputDuration(size_t samples) 
{
    if (m_d) m_d->setExpectedInputDuration(samples);
}

void
RubberBandStretcher::setMaxProcessSize(size_t samples)
{
    if (m_d) m_d->setMaxProcessSize(samples); //!!!  definitely need for r3d
}

void
RubberBandStretcher::setKeyFrameMap(const std::map<size_t, size_t> &mapping)
{
    if (m_d) m_d->setKeyFrameMap(mapping);
    //!!!
}

size_t
RubberBandStretcher::getSamplesRequired() const
{
    if (m_d) return m_d->getSamplesRequired();
    else return m_r3d->getSamplesRequired();
}

void
RubberBandStretcher::study(const float *const *input, size_t samples,
                           bool final)
{
    if (m_d) m_d->study(input, samples, final);
    //!!!
}

void
RubberBandStretcher::process(const float *const *input, size_t samples,
                             bool final)
{
    if (m_d) m_d->process(input, samples, final);
    else m_r3d->process(input, samples, final);
}

int
RubberBandStretcher::available() const
{
    if (m_d) return m_d->available();
    else return m_r3d->available();
}

size_t
RubberBandStretcher::retrieve(float *const *output, size_t samples) const
{
    if (m_d) return m_d->retrieve(output, samples);
    else return m_r3d->retrieve(output, samples);
}

float
RubberBandStretcher::getFrequencyCutoff(int n) const
{
    if (m_d) return m_d->getFrequencyCutoff(n);
    else return {};
}

void
RubberBandStretcher::setFrequencyCutoff(int n, float f) 
{
    if (m_d) m_d->setFrequencyCutoff(n, f);
}

size_t
RubberBandStretcher::getInputIncrement() const
{
    if (m_d) return m_d->getInputIncrement();
    else return {};
}

std::vector<int>
RubberBandStretcher::getOutputIncrements() const
{
    if (m_d) return m_d->getOutputIncrements();
    else return {};
}

std::vector<float>
RubberBandStretcher::getPhaseResetCurve() const
{
    if (m_d) return m_d->getPhaseResetCurve();
    else return {};
}

std::vector<int>
RubberBandStretcher::getExactTimePoints() const
{
    if (m_d) return m_d->getExactTimePoints();
    else return {};
}

size_t
RubberBandStretcher::getChannelCount() const
{
    if (m_d) return m_d->getChannelCount();
    else return m_r3d->getChannelCount();
}

void
RubberBandStretcher::calculateStretch()
{
    if (m_d) m_d->calculateStretch();
}

void
RubberBandStretcher::setDebugLevel(int level)
{
    if (m_d) m_d->setDebugLevel(level);
}

void
RubberBandStretcher::setDefaultDebugLevel(int level)
{
    Impl::setDefaultDebugLevel(level);
}

}

