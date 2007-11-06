/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

#include "StretcherChannelData.h"

#include "Resampler.h"

namespace RubberBand 
{

RubberBandStretcher::Impl::ChannelData::ChannelData(size_t blockSize,
                                                    size_t outbufSize)
{
    std::set<size_t> s;
    construct(s, blockSize, outbufSize);
}

RubberBandStretcher::Impl::ChannelData::ChannelData(const std::set<size_t> &blockSizes,
                                                    size_t initialBlockSize,
                                                    size_t outbufSize)
{
    construct(blockSizes, initialBlockSize, outbufSize);
}

void
RubberBandStretcher::Impl::ChannelData::construct(const std::set<size_t> &blockSizes,
                                                  size_t initialBlockSize,
                                                  size_t outbufSize)
{
    size_t maxSize = initialBlockSize;

    if (!blockSizes.empty()) {
        // std::set is ordered by value
        std::set<size_t>::const_iterator i = blockSizes.end();
        maxSize = *--i;
    }
    if (blockSizes.find(initialBlockSize) == blockSizes.end()) {
        if (initialBlockSize > maxSize) maxSize = initialBlockSize;
    }

    size_t realSize = maxSize/2 + 1; // size of the real "half" of freq data

    std::cerr << "ChannelData::construct([" << blockSizes.size() << "], " << maxSize << ", " << outbufSize << ")" << std::endl;
    
    if (outbufSize < maxSize) outbufSize = maxSize;

    inbuf = new RingBuffer<float>(maxSize);
    outbuf = new RingBuffer<float>(outbufSize);

    mag = new double[realSize];
    phase = new double[realSize];
    prevPhase = new double[realSize];
    unwrappedPhase = new double[realSize];
    freqPeak = new size_t[realSize];

    accumulator = new float[maxSize];
    windowAccumulator = new float[maxSize];

    fltbuf = new float[maxSize];
    dblbuf = new double[maxSize];

    for (std::set<size_t>::const_iterator i = blockSizes.begin();
         i != blockSizes.end(); ++i) {
        ffts[*i] = new FFT(*i);
        ffts[*i]->initDouble();
    }
    if (blockSizes.find(initialBlockSize) == blockSizes.end()) {
        ffts[initialBlockSize] = new FFT(initialBlockSize);
        ffts[initialBlockSize]->initDouble();
    }
    fft = ffts[initialBlockSize];

    resampler = 0;
    resamplebuf = 0;
    resamplebufSize = 0;

    reset();

    for (size_t i = 0; i < realSize; ++i) {
        mag[i] = 0.0;
        phase[i] = 0.0;
        prevPhase[i] = 0.0;
        unwrappedPhase[i] = 0.0;
        freqPeak[i] = 0;
    }

    for (size_t i = 0; i < maxSize; ++i) {
        accumulator[i] = 0.f;
        windowAccumulator[i] = 0.f;
        dblbuf[i] = 0.0;
        fltbuf[i] = 0.0;
    }
}

void
RubberBandStretcher::Impl::ChannelData::setBlockSize(size_t blockSize)
{
    size_t oldSize = inbuf->getSize();
    size_t realSize = blockSize/2 + 1;

    std::cerr << "ChannelData::setBlockSize(" << blockSize << ") [from " << oldSize << "]" << std::endl;

    if (oldSize >= blockSize) {

        // no need to reallocate buffers, just reselect fft

        //!!! we can't actually do this without locking against the
        //process thread, can we?  we need to zero the mag/phase
        //buffers without interference

        if (ffts.find(blockSize) == ffts.end()) {
            //!!! this also requires a lock, but it shouldn't occur in
            //RT mode with proper initialisation
            ffts[blockSize] = new FFT(blockSize);
            ffts[blockSize]->initDouble();
        }
        
        fft = ffts[blockSize];

        for (size_t i = 0; i < realSize; ++i) {
            mag[i] = 0.0;
            phase[i] = 0.0;
            prevPhase[i] = 0.0;
            unwrappedPhase[i] = 0.0;
            freqPeak[i] = 0;
        }

        return;
    }

    //!!! at this point we need a lock in case a different client
    //thread is calling process() -- we need this lock even if we
    //aren't running in threaded mode ourselves -- if we're in RT
    //mode, then the process call should trylock and fail if the lock
    //is unavailable (since this should never normally be the case in
    //general use in RT mode)

    RingBuffer<float> *newbuf = inbuf->resized(blockSize);
    delete inbuf;
    inbuf = newbuf;

    // We don't want to preserve data in these arrays

    delete[] mag;
    delete[] phase;
    delete[] prevPhase;
    delete[] unwrappedPhase;
    delete[] freqPeak;

    mag = new double[realSize];
    phase = new double[realSize];
    prevPhase = new double[realSize];
    unwrappedPhase = new double[realSize];
    freqPeak = new size_t[realSize];

    delete[] fltbuf;
    delete[] dblbuf;

    fltbuf = new float[blockSize];
    dblbuf = new double[blockSize];

    // But we do want to preserve data in these

    float *newAcc = new float[blockSize];
    for (size_t i = 0; i < oldSize; ++i) newAcc[i] = accumulator[i];
    delete[] accumulator;
    accumulator = newAcc;

    newAcc = new float[blockSize];
    for (size_t i = 0; i < oldSize; ++i) newAcc[i] = windowAccumulator[i];
    delete[] windowAccumulator;
    windowAccumulator = newAcc;
    
    //!!! and resampler?

    for (size_t i = 0; i < realSize; ++i) {
        mag[i] = 0.0;
        phase[i] = 0.0;
        prevPhase[i] = 0.0;
        unwrappedPhase[i] = 0.0;
        freqPeak[i] = 0;
    }

    for (size_t i = 0; i < blockSize; ++i) {
        dblbuf[i] = 0.0;
        fltbuf[i] = 0.0;
    }

    for (size_t i = oldSize; i < blockSize; ++i) {
        accumulator[i] = 0.f;
        windowAccumulator[i] = 0.f;
    }

    if (ffts.find(blockSize) == ffts.end()) {
        ffts[blockSize] = new FFT(blockSize);
        ffts[blockSize]->initDouble();
    }
    
    fft = ffts[blockSize];
}

void
RubberBandStretcher::Impl::ChannelData::setOutbufSize(size_t outbufSize)
{
    size_t oldSize = outbuf->getSize();

    std::cerr << "ChannelData::setOutbufSize(" << outbufSize << ") [from " << oldSize << "]" << std::endl;

    if (oldSize < outbufSize) {

        //!!! at this point we need a lock in case a different client
        //thread is calling process()

        //!!! this doesn't do what we want -- we want a locking resize
        //that preserves the existing data
//        outbuf->resize(outbufSize);

        RingBuffer<float> *newbuf = outbuf->resized(outbufSize);
        delete outbuf;
        outbuf = newbuf;
    }
}

RubberBandStretcher::Impl::ChannelData::~ChannelData()
{
    delete resampler;
    delete[] resamplebuf;

    delete inbuf;
    delete outbuf;
    delete[] mag;
    delete[] phase;
    delete[] prevPhase;
    delete[] unwrappedPhase;
    delete[] freqPeak;
    delete[] accumulator;
    delete[] windowAccumulator;
    delete[] fltbuf;
    delete[] dblbuf;

    for (std::map<size_t, FFT *>::iterator i = ffts.begin();
         i != ffts.end(); ++i) {
        delete i->second;
    }
}

void
RubberBandStretcher::Impl::ChannelData::reset()
{
    inbuf->reset();
    outbuf->reset();

    if (resampler) resampler->reset();

    accumulatorFill = 0;
    prevIncrement = 0;
    blockCount = 0;
    inCount = 0;
    inputSize = -1;
    outCount = 0;
    draining = false;
    outputComplete = false;
}

}
