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

#include "StretchCalculator.h"

#include <math.h>
#include <iostream>
#include <deque>
#include <set>
#include <cassert>

namespace RubberBand
{
	
StretchCalculator::StretchCalculator(size_t sampleRate,
                                     size_t inputIncrement,
                                     bool useHardPeaks) :
    m_sampleRate(sampleRate),
    m_increment(inputIncrement),
    m_prevDf(0),
    m_divergence(0),
    m_recovery(0),
    m_prevRatio(1.0),
    m_wasTransient(false),
    m_useHardPeaks(useHardPeaks)
{
    std::cerr << "StretchCalculator::StretchCalculator: useHardPeaks = " << useHardPeaks << std::endl;
}    

StretchCalculator::~StretchCalculator()
{
}

std::vector<int>
StretchCalculator::calculate(double ratio, size_t inputDuration,
                             const std::vector<float> &lockDf,
                             const std::vector<float> &stretchDf)
{
    // Method:

    //!!! This description is out of date.

    //!!! Rationalise naming -- generally wise to avoid the word
    //"frame" and instead use "block" / "sample" for processing frame /
    // audio frame.

    // 1. Pre-process the df array, and for each (say) one second's
    // worth of values, calculate the number of peaks that would
    // qualify for phase locking given the default threshold.  Then
    // reduce or increase the threshold by stages until that number is
    // in a sensible range (say 1-10 peaks per second -- the low end
    // is harder to estimate than the high end, so it may be better to
    // start with a high sensitivity and reduce it).

    // 2. Record the positions of peaks, and separately the positions
    // of those peaks that qualify for locking using the sliding
    // threshold window.  Don't permit two locked peaks within a very
    // short time frame (e.g. 30-50ms).

    // 3. Map each of the locked peaks (or any peaks that are over a
    // given intensity?), as well as the start and end points, to a
    // proportionate position in the newly stretched array so as to
    // ensure that their timing is strictly "correct".

    // 4. Calculate how much time is left in the stretch total, after
    // each of the locked frames has been allocated its static
    // allowance.  Also count the non-locked frames.

    // 5. For each region between two locked frames, calculate the
    // number of samples to allocate that region given the time
    // available for stretch and the number of non-locked frames.
    // Then divvy them up... how exactly?


    assert(lockDf.size() == stretchDf.size());
    
    m_lastPeaks = findPeaks(lockDf);
    std::vector<Peak> &peaks = m_lastPeaks;
    size_t totalCount = lockDf.size();

    std::vector<int> increments;

    size_t outputDuration = lrint(inputDuration * ratio);

    std::cerr << "debug level: " << m_debugLevel << std::endl;

    if (m_debugLevel > 0) {
        std::cerr << "StretchCalculator::calculate(): inputDuration " << inputDuration << ", ratio " << ratio << ", outputDuration " << outputDuration;
    }

    //!!! round down?
    outputDuration = lrint((lockDf.size() * m_increment) * ratio);

    if (m_debugLevel > 0) {
        std::cerr << " (rounded up to " << outputDuration << ")";
        std::cerr << ", df size " << lockDf.size() << std::endl;
    }

//    size_t stretchable = outputDuration - lockCount * m_increment;
    
    std::vector<size_t> fixedAudioFrames;
    for (size_t i = 0; i < peaks.size(); ++i) {
        fixedAudioFrames.push_back
            //!!! this should be rounding down, shouldn't it? not lrint?
            (lrint((double(peaks[i].frame) * outputDuration) / totalCount));
    }

//    size_t lockIndex = 0;

    if (m_debugLevel > 1) {
        std::cerr << "have " << peaks.size() << " fixed positions" << std::endl;
    }

    size_t totalInput = 0, totalOutput = 0;

    // so for each inter-lock region, we want to take the number of
    // output frames to be allocated and the detection function values
    // within the range, and produce a series of increments that sum
    // to the number of output frames, such that each increment is
    // displaced from the input increment by an amount inversely
    // proportional to the magnitude of the detection function at that
    // input step.  Ideally the detection function would have been
    // somewhat smoothed for this purpose but we'll start raw.

    //!!! Actually, we would possibly be better off using a fixed
    // smooth curve than the detection function itself.

    size_t regionTotalFrames = 0;

    for (size_t i = 0; i <= peaks.size(); ++i) {
        
        size_t regionStart, regionStartBlock, regionEnd, regionEndBlock;
        bool phaseLock = false;

        if (i == 0) {
            regionStartBlock = 0;
            regionStart = 0;
        } else {
            regionStartBlock = peaks[i-1].frame;
            regionStart = fixedAudioFrames[i-1];
            phaseLock = peaks[i-1].hard;
        }

        if (i == peaks.size()) {
            regionEndBlock = totalCount;
            regionEnd = outputDuration;
        } else {
            regionEndBlock = peaks[i].frame;
            regionEnd = fixedAudioFrames[i];
        }
        
        size_t regionDuration = regionEnd - regionStart;
        regionTotalFrames += regionDuration;

        std::vector<float> dfRegion;

        for (size_t j = regionStartBlock; j != regionEndBlock; ++j) {
            dfRegion.push_back(stretchDf[j]);
        }

        if (m_debugLevel > 1) {
            std::cerr << "distributeRegion from " << regionStartBlock << " to " << regionEndBlock << " (frames " << regionStart << " to " << regionEnd << ")" << std::endl;
        }

        dfRegion = smoothDF(dfRegion);
        
        std::vector<int> regionIncrements = distributeRegion
            (dfRegion, regionDuration, ratio, phaseLock);

        size_t totalForRegion = 0;

        for (size_t j = 0; j < regionIncrements.size(); ++j) {

            int incr = regionIncrements[j];

            if (j == 0 && phaseLock) increments.push_back(-incr);
            else increments.push_back(incr);

            if (incr > 0) totalForRegion += incr;
            else totalForRegion += -incr;

            totalInput += m_increment;
        }

        if (totalForRegion != regionDuration) {
            std::cerr << "*** WARNING: distributeRegion returned wrong duration " << totalForRegion << ", expected " << regionDuration << std::endl;
        }

        totalOutput += totalForRegion;
    }

    if (m_debugLevel > 0) {
        std::cerr << "total input increment = " << totalInput << " (= " << totalInput / m_increment << " blocks), output = " << totalOutput << ", ratio = " << double(totalOutput)/double(totalInput) << ", ideal output " << ceil(totalInput * ratio) << std::endl;
        std::cerr << "(region total = " << regionTotalFrames << ")" << std::endl;
    }
    return increments;
}

int
StretchCalculator::calculateSingle(double ratio,
                                   size_t inputDurationSoFar,
                                   float df)
{
    bool isTransient = false;

    //!!! We want to ensure, as close as possible, that the lock
    // points appear at _exactly_ the right frame numbers

    //!!! depends on block size.  larger block sizes need higher
    //thresholds.  since block size depends on ratio, I suppose we
    //could in theory calculate the threshold from the ratio directly.
    //For now we just frig it to work OK for a couple of common cases
    float transientThreshold = 0.35;
    if (ratio > 1) transientThreshold = 0.25;

    if (m_useHardPeaks && df > m_prevDf * 1.1 && df > transientThreshold) {
        isTransient = true;
    }

    if (m_debugLevel > 2) {
        std::cerr << "df = " << df << ", prevDf = " << m_prevDf
                  << ", thresh = " << transientThreshold << std::endl;
    }

    m_prevDf = df;

    if (isTransient && !m_wasTransient) {
        if (m_debugLevel > 1) {
            std::cerr << "StretchCalculator::calculateSingle: transient found at "
                      << inputDurationSoFar << std::endl;
        }
        m_divergence += m_increment - (m_increment * ratio);
        m_wasTransient = true;
        m_recovery = m_divergence / ((m_sampleRate / 10.0) / m_increment);
        return -m_increment;
    }

    if (m_prevRatio != ratio) {
        m_recovery = m_divergence / ((m_sampleRate / 10.0) / m_increment);
        m_prevRatio = ratio;
    }

    //!!! want transient amnesty as above (hard peak amnesty)
    m_wasTransient = false;

    int incr = lrint(m_increment * ratio - m_recovery);
    if (m_debugLevel > 2 || (m_debugLevel > 1 && m_divergence != 0)) {
        std::cerr << "divergence = " << m_divergence << ", recovery = " << m_recovery << ", incr = " << incr << ", ";
    }
    if (incr < (m_increment * ratio) / 2) {
        incr = (m_increment * ratio) / 2;
    } else if (incr > m_increment * ratio * 2) {
        incr = m_increment * ratio * 2;
    }

    double divdiff = (m_increment * ratio) - incr;

    if (m_debugLevel > 2 || (m_debugLevel > 1 && m_divergence != 0)) {
        std::cerr << "divdiff = " << divdiff << std::endl;
    }

    double prevDivergence = m_divergence;
    m_divergence -= divdiff;
    if ((prevDivergence < 0 && m_divergence > 0) ||
        (prevDivergence > 0 && m_divergence < 0)) {
        m_recovery = m_divergence / ((m_sampleRate / 10.0) / m_increment);
    }

    return incr;
}

void
StretchCalculator::reset()
{
    m_prevDf = 0;
    m_divergence = 0;
}

std::vector<StretchCalculator::Peak>
StretchCalculator::findPeaks(const std::vector<float> &rawDf)
{
    std::vector<float> df = smoothDF(rawDf);

    std::set<size_t> hardPeakCandidates;
    std::set<size_t> softPeakCandidates;

    if (m_useHardPeaks) {

    //!!! this should depend on duration based on output increment surely?
        size_t hardPeakAmnesty = lrint(ceil(double(m_sampleRate) /
                                            (20 * double(m_increment)))); // 0.05 sec ish
//    size_t hardPeakAmnesty = 5;

        size_t prevHardPeak = 0;
        std::cerr << "hardPeakAmnesty = " << hardPeakAmnesty << std::endl;
        for (size_t i = 1; i + 1 < df.size(); ++i) {

            //!!! this ratio configurable? dependent on block size and sr?

            if (df[i] < 0.1) continue;
            if (df[i] <= df[i-1] * 1.2) continue;

            if (df[i] > df[i-1] * 1.4 ||
                (df[i+1] > df[i] && df[i+1] > df[i-1] * 1.8) ||
                df[i] > 0.4) {
                if (!hardPeakCandidates.empty() &&
                    i < prevHardPeak + hardPeakAmnesty) {
                    continue;
                }
                size_t peakLocation = i;
                if (i + 1 < rawDf.size() &&
                    rawDf[i + 1] > rawDf[i] * 1.4) {
                    ++peakLocation;
                }
                if (m_debugLevel > 1) {
                    std::cerr << "hard peak at " << peakLocation << " (" << df[peakLocation] << " > " << df[peakLocation-1] << " * " << 1.4 << ")" << std::endl;
                }
                hardPeakCandidates.insert(peakLocation);
                prevHardPeak = peakLocation;
            }
        }
    }

    //!!! we don't yet do the right thing with soft peaks.  if
    //!useHardPeaks, we should be locking on soft peaks; if
    //useHardPeaks, we should be ignoring soft peaks if they occur
    //shortly after hard ones, otherwise either locking on them, or at
    //least making sure they fall at the correct sample time

//    int mediansize = lrint(ceil(double(m_sampleRate) /
//                                (4 * double(m_increment)))); // 0.25 sec ish
    size_t medianmaxsize = lrint(ceil(double(m_sampleRate) /
                                 double(m_increment))); // 1 sec ish
//    int mediansize = lrint(ceil(double(m_sampleRate) /
//                                (2 * double(m_increment)))); // 0.5 sec ish

    if (m_debugLevel > 1) {
        std::cerr << "mediansize = " << medianmaxsize << std::endl;
    }
    if (medianmaxsize < 7) {
        medianmaxsize = 7;
        if (m_debugLevel > 1) {
            std::cerr << "adjusted mediansize = " << medianmaxsize << std::endl;
        }
    }

    int minspacing = lrint(ceil(double(m_sampleRate) /
                                (20 * double(m_increment)))); // 0.05 sec ish
    
    std::deque<float> medianwin;
    std::vector<float> sorted;
    int softPeakAmnesty = 0;

    for (size_t i = 0; i < medianmaxsize/2; ++i) {
        medianwin.push_back(0);
    }
    for (size_t i = 0; i < medianmaxsize/2 && i < df.size(); ++i) {
        medianwin.push_back(df[i]);
    }

    size_t lastSoftPeak = 0;

    for (size_t i = 0; i < df.size(); ++i) {
        
        size_t mediansize = medianmaxsize;

        if (medianwin.size() < mediansize) {
            mediansize = medianwin.size();
        }

        size_t middle = medianmaxsize / 2;
        if (middle >= mediansize) middle = mediansize-1;

        size_t nextDf = i + mediansize - middle;

        if (mediansize < 2) {
            if (mediansize > medianmaxsize) { // absurd, but never mind that
//                std::cerr << "(<2) pop front ";
                medianwin.pop_front();
            }
            if (nextDf < df.size()) {
//                std::cerr << "(<2) push back " << df[nextDf] << " ";
                medianwin.push_back(df[nextDf]);
            } else {
                medianwin.push_back(0);
            }
//            std::cerr << "(<2) continue" << std::endl;
            continue;
        }

        if (m_debugLevel > 2) {
//            std::cerr << "have " << mediansize << " in median buffer" << std::endl;
        }

        sorted.clear();
        for (size_t j = 0; j < mediansize; ++j) {
            sorted.push_back(medianwin[j]);
        }
        std::sort(sorted.begin(), sorted.end());

        size_t n = 90; // percentile above which we pick peaks
        size_t index = (sorted.size() * n) / 100;
        if (index >= sorted.size()) index = sorted.size()-1;
        if (index == sorted.size()-1 && index > 0) --index;
        float thresh = sorted[index];

        if (m_debugLevel > 2) {
//            std::cerr << "medianwin[" << middle << "] = " << medianwin[middle] << ", thresh = " << thresh << std::endl;
            if (medianwin[middle] == 0.f) {
//                std::cerr << "contents: ";
                for (size_t j = 0; j < medianwin.size(); ++j) {
//                    std::cerr << medianwin[j] << " ";
                }
//                std::cerr << std::endl;
            }
        }

        if (medianwin[middle] > thresh &&
            medianwin[middle] > medianwin[middle-1] &&
            medianwin[middle] > medianwin[middle+1] &&
            softPeakAmnesty == 0) {

            size_t maxindex = middle;
            float maxval = medianwin[middle];

            for (size_t j = middle+1; j < mediansize; ++j) {
                if (medianwin[j] > maxval) {
                    maxval = medianwin[j];
                    maxindex = j;
                } else if (medianwin[j] < medianwin[middle]) {
                    break;
                }
            }

            //!!! we should distinguish between soft peaks (any found
            //using the above method) and hard peaks, which also show
            //a very rapid rise in detection function prior to the
            //peak (the first value after the rise is not necessarily
            //the peak itself, but it is probably where we should
            //locate the lock).  For hard peaks we need to lock in
            //time to preserve the shape of the transient (unless some
            //option is set to soft mode), for soft peaks we just want
            //to avoid poor timing positioning so we build up to the
            //lock at the exact peak moment.
            
//            size_t peak = i + maxindex - mediansize;
            size_t peak = i + maxindex - middle;

//            std::cerr << "i = " << i << ", maxindex = " << maxindex << ", middle = " << middle << ", so peak at " << peak << std::endl;

//            if (peak > 0) --peak; //!!! that's a fudge

            if (softPeakCandidates.empty() || lastSoftPeak != peak) {

                if (m_debugLevel > 1) {
                    std::cerr << "soft peak at " << peak << " (" << peak * m_increment << "): "
                              << medianwin[middle] << " > " << thresh << " and "
                              << medianwin[middle] << " > " << medianwin[middle-1] << " and "
                              << medianwin[middle] << " > " << medianwin[middle+1]
                              << std::endl;
                }

                if (peak >= df.size()) {
                    if (m_debugLevel > 2) {
                        std::cerr << "peak is beyond end"  << std::endl;
                    }
                } else {
                    softPeakCandidates.insert(peak);
                    lastSoftPeak = peak;
                }
            }

            softPeakAmnesty = minspacing + maxindex - middle;
            if (m_debugLevel > 2) {
                std::cerr << "amnesty = " << softPeakAmnesty << std::endl;
            }

        } else if (softPeakAmnesty > 0) --softPeakAmnesty;

//        std::cerr << "i = " << i << " ";
        if (mediansize >= medianmaxsize) {
//            std::cerr << "(>= " << medianmaxsize << ") pop front ";
            medianwin.pop_front();
        }
        if (nextDf < df.size()) {
//            std::cerr << "(" << nextDf << "<" << df.size() << ") push back " << df[nextDf] << " ";
            medianwin.push_back(df[nextDf]);
        } else {
            medianwin.push_back(0);
        }
//        std::cerr << "continue" << std::endl;
    }

    std::vector<Peak> peaks;

    //!!!
//    if (!softPeakCandidates.empty()) {
//        std::cerr << "clearing " << softPeakCandidates.size() << " soft peak candidates" << std::endl;
//    }
//    softPeakCandidates.clear();


    while (!hardPeakCandidates.empty() || !softPeakCandidates.empty()) {
        bool haveHardPeak = !hardPeakCandidates.empty();
        bool haveSoftPeak = !softPeakCandidates.empty();
        size_t hardPeak = (haveHardPeak ? *hardPeakCandidates.begin() : 0);
        size_t softPeak = (haveSoftPeak ? *softPeakCandidates.begin() : 0);
        Peak peak;
        peak.hard = false;
        peak.frame = softPeak;
        if (haveHardPeak &&
            (!haveSoftPeak || hardPeak <= softPeak)) {
            if (m_debugLevel > 2) {
                std::cerr << "Hard peak: " << hardPeak << std::endl;
            }
            peak.hard = true;
            peak.frame = hardPeak;
            hardPeakCandidates.erase(hardPeakCandidates.begin());
        } else {
            if (m_debugLevel > 2) {
                std::cerr << "Soft peak: " << softPeak << std::endl;
            }
        }            

        if (haveSoftPeak && peak.frame == softPeak) {
            softPeakCandidates.erase(softPeakCandidates.begin());
        }
            
        peaks.push_back(peak);
    }                

    return peaks;
}

std::vector<float>
StretchCalculator::smoothDF(const std::vector<float> &df)
{
    std::vector<float> smoothedDF;
    
    for (size_t i = 0; i < df.size(); ++i) {
        // three-value moving mean window for simple smoothing
        float total = 0.f, count = 0;
        if (i > 0) { total += df[i-1]; ++count; }
        total += df[i]; ++count;
        if (i+1 < df.size()) { total += df[i+1]; ++count; }
        float mean = total / count;
//        if (isnan(mean)) {
//            std::cerr << "ERROR: mean at " << i << " (of " << df.size() << ") is NaN: dfs are: "
//                      << df[i-1] << ", " << df[i] << ", " << df[i+1] << std::endl;
//        }
        smoothedDF.push_back(mean);
    }

    return smoothedDF;
}

std::vector<int>
StretchCalculator::distributeRegion(const std::vector<float> &dfIn,
                                    size_t duration, float ratio, bool lock)
{
    std::vector<float> df(dfIn);
    std::vector<int> increments;

    // The peak for the stretch detection function may appear after
    // the peak that we're using to calculate the start of the region.
    // We don't want that.  If we find a peak in the first half of
    // the region, we should set all the values up to that point to
    // the same value as the peak.

    //!!! this is not subtle enough, especially if the region is long
    //-- we want a bound that corresponds to acoustic perception of
    //the audible bounce

    for (size_t i = 1; i < df.size()/2; ++i) {
        if (df[i] < df[i-1]) {
            if (m_debugLevel > 1) {
                std::cerr << "stretch peak offset: " << i-1 << " (peak " << df[i-1] << ")" << std::endl;
            }
            for (size_t j = 0; j < i-1; ++j) {
                df[j] = df[i-1];
            }
            break;
        }
    }

    float maxDf = 0;

    for (size_t i = 0; i < df.size(); ++i) {
        if (i == 0 || df[i] > maxDf) maxDf = df[i];
    }

    // We want to try to ensure the last 100ms or so (if possible) are
    // tending back towards the maximum df, so that the stretchiness
    // reduces at the end of the stretched region.
    
    int reducedRegion = (0.1 * m_sampleRate) / m_increment;
    if (reducedRegion > df.size()/5) reducedRegion = df.size()/5;

    for (size_t i = 0; i < reducedRegion; ++i) {
        size_t index = df.size() - reducedRegion + i;
        df[index] = df[index] + ((maxDf - df[index]) * i) / reducedRegion;
    }

    long toAllot = long(duration) - long(m_increment * df.size());
//    bool negative = (toAllot < 0);
    
    if (m_debugLevel > 1) {
        std::cerr << "region of " << df.size() << " blocks, output duration " << duration << ", toAllot " << toAllot << std::endl;
    }

    size_t totalIncrement = 0;

    //!!! we need to place limits on the amount of displacement per
    //chunk.  if ratio < 0, no increment should be larger than
    //increment*ratio or smaller than increment*ratio/2; if ratio > 0,
    //none should be smaller than increment*ratio or larger than
    //increment*ratio*2.  We need to enforce this in the assignment of
    //displacements to allotments, not by trying to respond if
    //something turns out wrong

    //!!! ratio is only provided to this function for the purposes of
    //establishing this bound to the displacement
    
    // so if maxDisplacement / totalDisplacement > increment * ratio*2 - increment (for ratio > 1)
    // or maxDisplacement / totalDisplacement < increment * ratio/2 (for ratio < 1)
    // then we need to adjust... what?
    
    bool acceptableSquashRange = false;

    double totalDisplacement = 0;
    double maxDisplacement = 0; // min displacement will be 0 by definition

    maxDf = 0;
    float adj = 0;

    while (!acceptableSquashRange) {

        acceptableSquashRange = true;
        calculateDisplacements(df, maxDf, totalDisplacement, maxDisplacement,
                               adj);

        if (m_debugLevel > 1) {
            std::cerr << "totalDisplacement " << totalDisplacement << ", max " << maxDisplacement << " (maxDf " << maxDf << ", df count " << df.size() << ")" << std::endl;
        }

        if (totalDisplacement == 0) {
// Not usually a problem, in fact
//            std::cerr << "WARNING: totalDisplacement == 0 (duration " << duration << ", " << df.size() << " values in df)" << std::endl;
            if (!df.empty() && adj == 0) {
                acceptableSquashRange = false;
                adj = 1;
            }
            continue;
        }

        int extremeIncrement = m_increment + lrint((toAllot * maxDisplacement) / totalDisplacement);
        if (ratio < 1.0) {
            if (extremeIncrement > lrint(ceil(m_increment * ratio))) {
                std::cerr << "ERROR: extreme increment " << extremeIncrement << " > " << m_increment * ratio << " (I thought this couldn't happen?)" << std::endl;
            } else if (extremeIncrement < (m_increment * ratio) / 2) {
                if (m_debugLevel > 0) {
                    std::cerr << "WARNING: extreme increment " << extremeIncrement << " < " << (m_increment * ratio) / 2 << std::endl;
                }
                acceptableSquashRange = false;
            }
        } else {
            if (extremeIncrement > m_increment * ratio * 2) {
                if (m_debugLevel > 0) {
                    std::cerr << "WARNING: extreme increment " << extremeIncrement << " > " << m_increment * ratio * 2 << std::endl;
                }
                acceptableSquashRange = false;
            } else if (extremeIncrement < lrint(floor(m_increment * ratio))) {
                std::cerr << "ERROR: extreme increment " << extremeIncrement << " < " << m_increment * ratio << " (I thought this couldn't happen?)" << std::endl;
            }
        }

        if (!acceptableSquashRange) {
            // Need to make maxDisplacement smaller as a proportion of
            // the total displacement, yet ensure that the
            // displacements still sum to the total.  How?

//            std::cerr << "Adjusting df values by " << maxDf/10 << "..." << std::endl;

//            std::cerr << "now: ";
//            for (size_t i = 0; i < df.size(); ++i) {
//                df[i] += maxDf/10;
//                std::cerr << df[i] << " ";
//            }
//            std::cerr << std::endl;
            adj += maxDf/10;
            
            //...
        }
    }

    for (size_t i = 0; i < df.size(); ++i) {

        double displacement = maxDf - df[i];
        if (displacement < 0) displacement -= adj;
        else displacement += adj;

        if (i == 0 && lock) {
            if (df.size() == 1) {
                increments.push_back(duration);
                totalIncrement += duration;
            } else {
                increments.push_back(m_increment);
                totalIncrement += m_increment;
            }
            totalDisplacement -= displacement;
            continue;
        }

        double theoreticalAllotment = 0;

        if (totalDisplacement != 0) {
            theoreticalAllotment = (toAllot * displacement) / totalDisplacement;
        }
        int allotment = lrint(theoreticalAllotment);
        if (i + 1 == df.size()) allotment = toAllot;

        int increment = m_increment + allotment;

        if (increment <= 0) {
            //!!! this is a serious problem, the allocation is quite wrong if it allows increment to diverge so far from the input increment
            std::cerr << "*** WARNING: increment " << increment << " <= 0, rounding to zero" << std::endl;
            increment = 0;
            allotment = increment - m_increment;
        }

        increments.push_back(increment);
        totalIncrement += increment;

        toAllot -= allotment;
        totalDisplacement -= displacement;

        if (m_debugLevel > 2) {
            std::cerr << "df " << df[i] << ", smoothed " << df[i] << ", disp " << displacement << ", allot " << theoreticalAllotment << ", incr " << increment << ", remain " << toAllot << std::endl;
        }
    }
    
    if (m_debugLevel > 2) {
        std::cerr << "total increment: " << totalIncrement << ", left over: " << toAllot << " to allot, displacement " << totalDisplacement << std::endl;
    }

    if (totalIncrement != duration) {
        std::cerr << "*** WARNING: calculated output duration " << totalIncrement << " != expected " << duration << std::endl;
    }

    return increments;
}

void
StretchCalculator::calculateDisplacements(const std::vector<float> &df,
                                          float &maxDf,
                                          double &totalDisplacement,
                                          double &maxDisplacement,
                                          float adj) const
{
    totalDisplacement = maxDisplacement = 0;

    maxDf = 0;

    for (size_t i = 0; i < df.size(); ++i) {
        if (i == 0 || df[i] > maxDf) maxDf = df[i];
    }

    for (size_t i = 0; i < df.size(); ++i) {
        double displacement = maxDf - df[i];
        if (displacement < 0) displacement -= adj;
        else displacement += adj;
        totalDisplacement += displacement;
        if (i == 0 || displacement > maxDisplacement) {
            maxDisplacement = displacement;
        }
    }
}

}

