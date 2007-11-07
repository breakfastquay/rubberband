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

#include "RubberBandVampPlugin.h"

#include "StretchCalculator.h"

#include <cmath>

using std::string;
using std::vector;
using std::cerr;
using std::endl;

class RubberBandVampPlugin::Impl
{
public:
    size_t m_stepSize;
    size_t m_blockSize;
    size_t m_sampleRate;

    float m_timeRatio;
    float m_pitchRatio;

    bool m_realtime;
    bool m_elasticTiming;
    bool m_crispTransients;
    bool m_threadingAllowed;

    RubberBand::RubberBandStretcher *m_stretcher;

    int m_incrementsOutput;
    int m_aggregateIncrementsOutput;
    int m_divergenceOutput;
    int m_lockDfOutput;
    int m_smoothedLockDfOutput;
    int m_lockPointsOutput;

    size_t m_counter;
    size_t m_accumulatedIncrement;

    FeatureSet processOffline(const float *const *inputBuffers,
                              Vamp::RealTime timestamp);

    FeatureSet getRemainingFeaturesOffline();

    FeatureSet processRealTime(const float *const *inputBuffers,
                               Vamp::RealTime timestamp);

    FeatureSet getRemainingFeaturesRealTime();

    FeatureSet createFeatures(size_t inputIncrement,
                              std::vector<int> &outputIncrements,
                              std::vector<float> &lockDf,
                              std::vector<int> &exactPoints,
                              std::vector<float> &smoothedDF,
                              size_t baseCount,
                              bool includeFinal);
};


RubberBandVampPlugin::RubberBandVampPlugin(float inputSampleRate) :
    Plugin(inputSampleRate)
{
    m_d = new Impl();
    m_d->m_stepSize = 0;
    m_d->m_timeRatio = 1.0;
    m_d->m_pitchRatio = 1.0;
    m_d->m_realtime = false;
    m_d->m_elasticTiming = true;
    m_d->m_crispTransients = true;
    m_d->m_threadingAllowed = true;
    m_d->m_stretcher = 0;
    m_d->m_sampleRate = lrintf(m_inputSampleRate);
}

RubberBandVampPlugin::~RubberBandVampPlugin()
{
    delete m_d->m_stretcher;
    delete m_d;
}

string
RubberBandVampPlugin::getIdentifier() const
{
    return "rubberband";
}

string
RubberBandVampPlugin::getName() const
{
    return "Rubber Band Timestretch Analysis";
}

string
RubberBandVampPlugin::getDescription() const
{
    return "Carry out analysis phases of time stretcher process";
}

string
RubberBandVampPlugin::getMaker() const
{
    return "Rubber Band"; ///!!!
}

int
RubberBandVampPlugin::getPluginVersion() const
{
    return 1;
}

string
RubberBandVampPlugin::getCopyright() const
{
    return "";//!!!
}

RubberBandVampPlugin::OutputList
RubberBandVampPlugin::getOutputDescriptors() const
{
    OutputList list;

    size_t rate = 0;
    if (m_d->m_stretcher) {
        rate = lrintf(m_inputSampleRate / m_d->m_stretcher->getInputIncrement());
    }

    OutputDescriptor d;
    d.identifier = "increments";
    d.name = "Output Increments";
    d.description = ""; //!!!
    d.unit = "samples";
    d.hasFixedBinCount = true;
    d.binCount = 1;
    d.hasKnownExtents = false; //!!!
    d.isQuantized = true;
    d.quantizeStep = 1.0;
    d.sampleType = OutputDescriptor::VariableSampleRate;
    d.sampleRate = rate;
    m_d->m_incrementsOutput = list.size();
    list.push_back(d);

    d.identifier = "aggregate_increments";
    d.name = "Accumulated Output Increments";
    d.description = ""; //!!!
    d.sampleRate = 0;
    m_d->m_aggregateIncrementsOutput = list.size();
    list.push_back(d);

    d.identifier = "divergence";
    d.name = "Divergence from Linear";
    d.description = ""; //!!!
    d.isQuantized = false;
    d.sampleRate = 0;
    m_d->m_divergenceOutput = list.size();
    list.push_back(d);

    d.identifier = "lockdf";
    d.name = "Lock Point Detection Function";
    d.description = ""; //!!!
    d.unit = "";
    d.sampleRate = rate;
    m_d->m_lockDfOutput = list.size();
    list.push_back(d);

    d.identifier = "smoothedlockdf";
    d.name = "Smoothed Lock Point Detection Function";
    d.description = ""; //!!!
    d.unit = "";
    m_d->m_smoothedLockDfOutput = list.size();
    list.push_back(d);

    d.identifier = "lockpoints";
    d.name = "Phase Lock Points";
    d.description = ""; //!!!
    d.unit = "";
    d.hasFixedBinCount = true;
    d.binCount = 0;
    d.hasKnownExtents = false;
    d.isQuantized = false;
    d.sampleRate = 0;
    m_d->m_lockPointsOutput = list.size();
    list.push_back(d);

    return list;
}

RubberBandVampPlugin::ParameterList
RubberBandVampPlugin::getParameterDescriptors() const
{
    ParameterList list;

    ParameterDescriptor d;
    d.identifier = "timeratio";
    d.name = "Timestretch Ratio";
    d.description = ""; //!!!
    d.unit = "";
    d.minValue = 0.0000001;
    d.maxValue = 1000;
    d.defaultValue = 1.0;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "pitchratio";
    d.name = "Pitch Scaling Ratio";
    d.description = ""; //!!!
    d.unit = "";
    d.minValue = 0.0000001;
    d.maxValue = 1000;
    d.defaultValue = 1.0;
    d.isQuantized = false;
    list.push_back(d);

    d.identifier = "mode";
    d.name = "Processing Mode";
    d.description = ""; //!!!
    d.unit = "";
    d.minValue = 0;
    d.maxValue = 1;
    d.defaultValue = 0;
    d.isQuantized = true;
    d.quantizeStep = 1;
    d.valueNames.clear();
    d.valueNames.push_back("Offline");
    d.valueNames.push_back("Real Time");
    list.push_back(d);

    d.identifier = "stretchtype";
    d.name = "Stretch Flexibility";
    d.description = ""; //!!!
    d.unit = "";
    d.minValue = 0;
    d.maxValue = 1;
    d.defaultValue = 0;
    d.isQuantized = true;
    d.quantizeStep = 1;
    d.valueNames.clear();
    d.valueNames.push_back("Elastic");
    d.valueNames.push_back("Precise");
    list.push_back(d);

    d.identifier = "transientmode";
    d.name = "Transient Handling";
    d.description = ""; //!!!
    d.unit = "";
    d.minValue = 0;
    d.maxValue = 1;
    d.defaultValue = 0;
    d.isQuantized = true;
    d.quantizeStep = 1;
    d.valueNames.clear();
    d.valueNames.push_back("Crisp");
    d.valueNames.push_back("Soft");
    list.push_back(d);

    d.identifier = "threadingmode";
    d.name = "Threading";
    d.description = ""; //!!!
    d.unit = "";
    d.minValue = 0;
    d.maxValue = 1;
    d.defaultValue = 0;
    d.isQuantized = true;
    d.quantizeStep = 1;
    d.valueNames.clear();
    d.valueNames.push_back("Enabled");
    d.valueNames.push_back("Disabled");
    list.push_back(d);

    return list;
}

float
RubberBandVampPlugin::getParameter(std::string id) const
{
    if (id == "timeratio") return m_d->m_timeRatio;
    if (id == "pitchratio") return m_d->m_pitchRatio;
    if (id == "mode") return m_d->m_realtime ? 1 : 0;
    if (id == "stretchtype") return m_d->m_elasticTiming ? 0 : 1;
    if (id == "transientmode") return m_d->m_crispTransients ? 0 : 1;
    if (id == "threadingmode") return m_d->m_threadingAllowed ? 0 : 1;
    return 0.f;
}

void
RubberBandVampPlugin::setParameter(std::string id, float value)
{
    if (id == "timeratio") {
        m_d->m_timeRatio = value;
    } else if (id == "pitchratio") {
        m_d->m_pitchRatio = value;
    } else {
        bool set = (value > 0.5);
        if (id == "mode") m_d->m_realtime = set;
        else if (id == "stretchtype") m_d->m_elasticTiming = !set;
        else if (id == "transientmode") m_d->m_crispTransients = !set;
        else if (id == "threadingmode") m_d->m_threadingAllowed = !set;
    }
}

bool
RubberBandVampPlugin::initialise(size_t channels, size_t stepSize, size_t blockSize)
{
    if (channels < getMinChannelCount() ||
	channels > getMaxChannelCount()) return false;

    m_d->m_stepSize = std::min(stepSize, blockSize);
    m_d->m_blockSize = stepSize;

    RubberBand::RubberBandStretcher::Options options = 0;

    if (m_d->m_realtime)
         options |= RubberBand::RubberBandStretcher::OptionProcessRealTime;
    else options |= RubberBand::RubberBandStretcher::OptionProcessOffline;

    if (m_d->m_elasticTiming)
         options |= RubberBand::RubberBandStretcher::OptionStretchElastic;
    else options |= RubberBand::RubberBandStretcher::OptionStretchPrecise;
 
    if (m_d->m_crispTransients) 
         options |= RubberBand::RubberBandStretcher::OptionTransientsCrisp;
    else options |= RubberBand::RubberBandStretcher::OptionTransientsSmooth;

    if (m_d->m_threadingAllowed)
         options |= RubberBand::RubberBandStretcher::OptionThreadingAuto;
    else options |= RubberBand::RubberBandStretcher::OptionThreadingNone;

    delete m_d->m_stretcher;
    m_d->m_stretcher = new RubberBand::RubberBandStretcher
        (m_d->m_sampleRate, channels, options);
    m_d->m_stretcher->setDebugLevel(2);
    m_d->m_stretcher->setTimeRatio(m_d->m_timeRatio);
    m_d->m_stretcher->setPitchScale(m_d->m_pitchRatio);

    m_d->m_counter = 0;
    m_d->m_accumulatedIncrement = 0;

    return true;
}

void
RubberBandVampPlugin::reset()
{
//    delete m_stretcher;  //!!! or just if (m_stretcher) m_stretcher->reset();
//    m_stretcher = new RubberBand::RubberBandStretcher(lrintf(m_inputSampleRate), channels);
    if (m_d->m_stretcher) m_d->m_stretcher->reset();
}

RubberBandVampPlugin::FeatureSet
RubberBandVampPlugin::process(const float *const *inputBuffers,
                              Vamp::RealTime timestamp)
{
    if (m_d->m_realtime) {
        return m_d->processRealTime(inputBuffers, timestamp);
    } else {
        return m_d->processOffline(inputBuffers, timestamp);
    }        
}

RubberBandVampPlugin::FeatureSet
RubberBandVampPlugin::getRemainingFeatures()
{
    if (m_d->m_realtime) {
        return m_d->getRemainingFeaturesRealTime();
    } else {
        return m_d->getRemainingFeaturesOffline();
    }
}

RubberBandVampPlugin::FeatureSet
RubberBandVampPlugin::Impl::processOffline(const float *const *inputBuffers,
                                           Vamp::RealTime timestamp)
{
    if (!m_stretcher) {
	cerr << "ERROR: RubberBandVampPlugin::processOffline: "
	     << "RubberBandVampPlugin has not been initialised"
	     << endl;
	return FeatureSet();
    }

    m_stretcher->study(inputBuffers, m_blockSize, false);
    return FeatureSet();
}

RubberBandVampPlugin::FeatureSet
RubberBandVampPlugin::Impl::getRemainingFeaturesOffline()
{
    m_stretcher->study(0, 0, true);

    m_stretcher->calculateStretch();

    int rate = m_sampleRate;

    RubberBand::StretchCalculator sc(rate,
                                     m_stretcher->getInputIncrement(),
                                     true);

    size_t inputIncrement = m_stretcher->getInputIncrement();
    std::vector<int> outputIncrements = m_stretcher->getOutputIncrements();
    std::vector<float> lockDf = m_stretcher->getLockCurve();
    std::vector<int> peaks = m_stretcher->getExactTimePoints();
    std::vector<float> smoothedDf = sc.smoothDF(lockDf);

    FeatureSet features = createFeatures
        (inputIncrement, outputIncrements, lockDf, peaks, smoothedDf,
         0, true);

    return features;
}

RubberBandVampPlugin::FeatureSet
RubberBandVampPlugin::Impl::processRealTime(const float *const *inputBuffers,
                                            Vamp::RealTime timestamp)
{
    // This function is not in any way a real-time function (i.e. it
    // has no requirement to be RT safe); it simply operates the
    // stretcher in RT mode.

    if (!m_stretcher) {
	cerr << "ERROR: RubberBandVampPlugin::processRealTime: "
	     << "RubberBandVampPlugin has not been initialised"
	     << endl;
	return FeatureSet();
    }

    m_stretcher->process(inputBuffers, m_blockSize, false);
    
    size_t inputIncrement = m_stretcher->getInputIncrement();
    std::vector<int> outputIncrements = m_stretcher->getOutputIncrements();
    std::vector<float> lockDf = m_stretcher->getLockCurve();
    std::vector<float> smoothedDf; // not meaningful in RT mode
    std::vector<int> dummyPoints;
    FeatureSet features = createFeatures
        (inputIncrement, outputIncrements, lockDf, dummyPoints, smoothedDf, 
         m_counter, false);
    m_counter += outputIncrements.size();

    return features;
}

RubberBandVampPlugin::FeatureSet
RubberBandVampPlugin::Impl::getRemainingFeaturesRealTime()
{
    return FeatureSet();
}

RubberBandVampPlugin::FeatureSet
RubberBandVampPlugin::Impl::createFeatures(size_t inputIncrement,
                                           std::vector<int> &outputIncrements,
                                           std::vector<float> &lockDf,
                                           std::vector<int> &exactPoints,
                                           std::vector<float> &smoothedDf,
                                           size_t baseCount,
                                           bool includeFinal)
{
    size_t actual = m_accumulatedIncrement;

    double overallRatio = m_timeRatio * m_pitchRatio;

    char label[200];

    FeatureSet features;

    int rate = m_sampleRate;

    size_t epi = 0;

    for (size_t i = 0; i < outputIncrements.size(); ++i) {

        size_t frame = (baseCount + i) * inputIncrement;

        int oi = outputIncrements[i];
        bool hardLock = false;
        bool softLock = false;

        if (oi < 0) {
            oi = -oi;
            hardLock = true;
        }

        if (epi < exactPoints.size() && int(i) == exactPoints[epi]) {
            softLock = true;
            ++epi;
        }

        double linear = (frame * overallRatio);

        Vamp::RealTime t = Vamp::RealTime::frame2RealTime(frame, rate);

        Feature feature;
        feature.hasTimestamp = true;
        feature.timestamp = t;
        feature.values.push_back(oi);
        feature.label = Vamp::RealTime::frame2RealTime(oi, rate).toText();
        features[m_incrementsOutput].push_back(feature);

        feature.values.clear();
        feature.values.push_back(actual);
        feature.label = Vamp::RealTime::frame2RealTime(actual, rate).toText();
        features[m_aggregateIncrementsOutput].push_back(feature);

        feature.values.clear();
        feature.values.push_back(actual - linear);

        sprintf(label, "expected %ld, actual %ld, difference %ld (%s ms)",
                long(linear), long(actual), long(actual - linear),
                // frame2RealTime expects an integer frame number,
                // hence our multiplication factor
                (Vamp::RealTime::frame2RealTime
                 (lrintf((actual - linear) * 1000), rate) / 1000)
                .toText().c_str());
        feature.label = label;

        features[m_divergenceOutput].push_back(feature);
        actual += oi;
        
        char buf[30];

        if (i < lockDf.size()) {
            feature.values.clear();
            feature.values.push_back(lockDf[i]);
            sprintf(buf, "%d", baseCount + i);
            feature.label = buf;
            features[m_lockDfOutput].push_back(feature);
        }

        if (i < smoothedDf.size()) {
            feature.values.clear();
            feature.values.push_back(smoothedDf[i]);
            features[m_smoothedLockDfOutput].push_back(feature);
        }

        if (hardLock) {
            feature.values.clear();
            feature.label = "Phase Reset";
            features[m_lockPointsOutput].push_back(feature);
        } else if (softLock) {
            feature.values.clear();
            feature.label = "Time Sync";
            features[m_lockPointsOutput].push_back(feature);
        }            
    }

    if (includeFinal) {
        Vamp::RealTime t = Vamp::RealTime::frame2RealTime
            (inputIncrement * (baseCount + outputIncrements.size()), rate);
        Feature feature;
        feature.hasTimestamp = true;
        feature.timestamp = t;
        feature.label = Vamp::RealTime::frame2RealTime(actual, rate).toText();
        feature.values.clear();
        feature.values.push_back(actual);
        features[m_aggregateIncrementsOutput].push_back(feature);

        float linear = ((baseCount + outputIncrements.size())
                        * inputIncrement * overallRatio);
        feature.values.clear();
        feature.values.push_back(actual - linear);
        feature.label =  // see earlier comment
            (Vamp::RealTime::frame2RealTime //!!! update this as earlier label
             (lrintf((actual - linear) * 1000), rate) / 1000)
            .toText();
        features[m_divergenceOutput].push_back(feature);
    }

    m_accumulatedIncrement = actual;

    return features;
}

