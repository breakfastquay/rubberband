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

#ifndef RUBBERBAND_R3_STRETCHERIMPL_H
#define RUBBERBAND_R3_STRETCHERIMPL_H

#include "BinSegmenter.h"
#include "Guide.h"
#include "Peak.h"
#include "PhaseAdvance.h"

#include "../common/StretchCalculator.h"
#include "../common/Resampler.h"
#include "../common/FFT.h"
#include "../common/FixedVector.h"
#include "../common/Allocators.h"
#include "../common/Window.h"

#include "../../rubberband/RubberBandStretcher.h"

#include <map>
#include <memory>
#include <functional>

namespace RubberBand
{

class R3StretcherImpl
{
public:
    struct Parameters {
        double sampleRate;
        int channels;
        RubberBandStretcher::Options options;
        std::function<void(const std::string &)> logger;
        Parameters(double _sampleRate, int _channels,
                   RubberBandStretcher::Options _options,
                   std::function<void(const std::string &)> _log = &logCout) :
            sampleRate(_sampleRate), channels(_channels), options(_options),
            logger(_log) { }
    };
    
    R3StretcherImpl(Parameters parameters,
                    double initialTimeRatio,
                    double initialPitchScale);
    ~R3StretcherImpl() { }

    void reset();
    
    void setTimeRatio(double ratio);
    void setPitchScale(double scale);
    void setFormantScale(double scale);

    double getTimeRatio() const;
    double getPitchScale() const;
    double getFormantScale() const;

    void setKeyFrameMap(const std::map<size_t, size_t> &);

    void setFormantOption(RubberBandStretcher::Options);
    
    void study(const float *const *input, size_t samples, bool final);
    size_t getSamplesRequired() const;
    void process(const float *const *input, size_t samples, bool final);
    int available() const;
    size_t retrieve(float *const *output, size_t samples) const;

    size_t getLatency() const;
    size_t getChannelCount() const;
    
protected:
    struct ClassificationReadaheadData {
        FixedVector<double> timeDomain;
        FixedVector<double> mag;
        FixedVector<double> phase;
        ClassificationReadaheadData(int _fftSize) :
            timeDomain(_fftSize, 0.f),
            mag(_fftSize/2 + 1, 0.f),
            phase(_fftSize/2 + 1, 0.f)
        { }

    private:
        ClassificationReadaheadData(const ClassificationReadaheadData &) =delete;
        ClassificationReadaheadData &operator=(const ClassificationReadaheadData &) =delete;
    };
    
    struct ChannelScaleData {
        int fftSize;
        int bufSize; // size of every freq-domain array here: fftSize/2 + 1
        FixedVector<double> timeDomain;
        FixedVector<double> real;
        FixedVector<double> imag;
        FixedVector<double> mag;
        FixedVector<double> phase;
        FixedVector<double> advancedPhase;
        FixedVector<double> prevMag;
        FixedVector<double> pendingKick;
        FixedVector<double> accumulator;

        ChannelScaleData(int _fftSize, int _longestFftSize) :
            fftSize(_fftSize),
            bufSize(fftSize/2 + 1),
            timeDomain(fftSize, 0.f),
            real(bufSize, 0.f),
            imag(bufSize, 0.f),
            mag(bufSize, 0.f),
            phase(bufSize, 0.f),
            advancedPhase(bufSize, 0.f),
            prevMag(bufSize, 0.f),
            pendingKick(bufSize, 0.f),
            accumulator(_longestFftSize, 0.f)
        { }

        void reset() {
            v_zero(prevMag.data(), prevMag.size());
            v_zero(accumulator.data(), accumulator.size());
        }

    private:
        ChannelScaleData(const ChannelScaleData &) =delete;
        ChannelScaleData &operator=(const ChannelScaleData &) =delete;
    };

    struct FormantData {
        int fftSize;
        FixedVector<double> cepstra;
        FixedVector<double> envelope;
        FixedVector<double> spare;

        FormantData(int _fftSize) :
            fftSize(_fftSize),
            cepstra(_fftSize, 0.0),
            envelope(_fftSize/2 + 1, 0.0),
            spare(_fftSize/2 + 1, 0.0) { }

        double envelopeAt(double bin) const {
            int b0 = int(floor(bin)), b1 = int(ceil(bin));
            if (b0 < 0 || b0 > fftSize/2) {
                return 0.0;
            } else if (b1 == b0 || b1 > fftSize/2) {
                return envelope.at(b0);
            } else {
                double diff = bin - double(b0);
                return envelope.at(b0) * (1.0 - diff) + envelope.at(b1) * diff;
            }
        }
    };

    struct ChannelData {
        std::map<int, std::shared_ptr<ChannelScaleData>> scales;
        ClassificationReadaheadData readahead;
        bool haveReadahead;
        std::unique_ptr<BinClassifier> classifier;
        FixedVector<BinClassifier::Classification> classification;
        FixedVector<BinClassifier::Classification> nextClassification;
        std::unique_ptr<BinSegmenter> segmenter;
        BinSegmenter::Segmentation segmentation;
        BinSegmenter::Segmentation prevSegmentation;
        BinSegmenter::Segmentation nextSegmentation;
        Guide::Guidance guidance;
        FixedVector<float> mixdown;
        FixedVector<float> resampled;
        std::unique_ptr<RingBuffer<float>> inbuf;
        std::unique_ptr<RingBuffer<float>> outbuf;
        std::unique_ptr<FormantData> formant;
        ChannelData(BinSegmenter::Parameters segmenterParameters,
                    BinClassifier::Parameters classifierParameters,
                    int longestFftSize,
                    int inRingBufferSize,
                    int outRingBufferSize) :
            scales(),
            readahead(segmenterParameters.fftSize),
            haveReadahead(false),
            classifier(new BinClassifier(classifierParameters)),
            classification(classifierParameters.binCount,
                           BinClassifier::Classification::Silent),
            nextClassification(classifierParameters.binCount,
                               BinClassifier::Classification::Silent),
            segmenter(new BinSegmenter(segmenterParameters)),
            segmentation(), prevSegmentation(), nextSegmentation(),
            mixdown(longestFftSize, 0.f), // though it could be shorter
            resampled(outRingBufferSize, 0.f),
            inbuf(new RingBuffer<float>(inRingBufferSize)),
            outbuf(new RingBuffer<float>(outRingBufferSize)),
            formant(new FormantData(segmenterParameters.fftSize)) { }
        void reset() {
            haveReadahead = false;
            classifier->reset();
            segmentation = BinSegmenter::Segmentation();
            prevSegmentation = BinSegmenter::Segmentation();
            nextSegmentation = BinSegmenter::Segmentation();
            inbuf->reset();
            outbuf->reset();
            for (auto &s : scales) {
                s.second->reset();
            }
        }
    };

    struct ChannelAssembly {
        // Vectors of bare pointers, used to package container data
        // from different channels into arguments for PhaseAdvance
        FixedVector<double *> mag;
        FixedVector<double *> phase;
        FixedVector<double *> prevMag;
        FixedVector<Guide::Guidance *> guidance;
        FixedVector<double *> outPhase;
        FixedVector<float *> mixdown;
        FixedVector<float *> resampled;
        ChannelAssembly(int channels) :
            mag(channels, nullptr), phase(channels, nullptr),
            prevMag(channels, nullptr), guidance(channels, nullptr),
            outPhase(channels, nullptr), mixdown(channels, nullptr),
            resampled(channels, nullptr) { }
    };

    struct ScaleData {
        int fftSize;
        FFT fft;
        Window<double> analysisWindow;
        Window<double> synthesisWindow;
        double windowScaleFactor;
        GuidedPhaseAdvance guided;
        ScaleData(GuidedPhaseAdvance::Parameters guidedParameters) :
            fftSize(guidedParameters.fftSize),
            fft(fftSize),
            analysisWindow(analysisWindowShape(fftSize),
                           analysisWindowLength(fftSize)),
            synthesisWindow(synthesisWindowShape(fftSize),
                            synthesisWindowLength(fftSize)),
            guided(guidedParameters),
            windowScaleFactor(0.0)
        {
            int asz = analysisWindow.getSize(), ssz = synthesisWindow.getSize();
            int off = (asz - ssz) / 2;
            for (int i = 0; i < ssz; ++i) {
                windowScaleFactor += analysisWindow.getValue(i + off) *
                    synthesisWindow.getValue(i);
            }
        }

        WindowType analysisWindowShape(int fftSize);
        int analysisWindowLength(int fftSize);
        WindowType synthesisWindowShape(int fftSize);
        int synthesisWindowLength(int fftSize);
    };
    
    Parameters m_parameters;

    std::atomic<double> m_timeRatio;
    std::atomic<double> m_pitchScale;
    std::atomic<double> m_formantScale;
    
    std::vector<std::shared_ptr<ChannelData>> m_channelData;
    std::map<int, std::shared_ptr<ScaleData>> m_scaleData;
    Guide m_guide;
    Guide::Configuration m_guideConfiguration;
    ChannelAssembly m_channelAssembly;
    std::unique_ptr<StretchCalculator> m_calculator;
    std::unique_ptr<Resampler> m_resampler;
    std::atomic<int> m_inhop;
    int m_prevInhop;
    int m_prevOuthop;
    int m_startSkip;

    size_t m_studyInputDuration;
    size_t m_totalTargetDuration;
    size_t m_processInputDuration;
    size_t m_totalOutputDuration;
    std::map<size_t, size_t> m_keyFrameMap;
    
    enum class ProcessMode {
        JustCreated,
        Studying,
        Processing,
        Finished
    };
    ProcessMode m_mode;

    void consume();
    void calculateHop();
    void updateRatioFromMap();
    void analyseChannel(int channel, int inhop, int prevInhop, int prevOuthop);
    void analyseFormant(int channel);
    void adjustFormant(int channel);
    void adjustPreKick(int channel);
    void synthesiseChannel(int channel, int outhop);

    double getEffectiveRatio() const {
        return m_timeRatio * m_pitchScale;
    }

    bool isRealTime() const {
        return m_parameters.options &
            RubberBandStretcher::OptionProcessRealTime;
    }
    
    static void logCout(const std::string &message) {
        std::cout << "RubberBandStretcher: " << message << std::endl;
    }
};

}

#endif
