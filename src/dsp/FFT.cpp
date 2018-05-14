/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band Library
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2018 Particular Programs Ltd.

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

#include "FFT.h"
#include "system/Thread.h"
#include "base/Profiler.h"
#include "system/Allocators.h"
#include "system/VectorOps.h"
#include "system/VectorOpsComplex.h"

//#define FFT_MEASUREMENT 1

#ifdef FFT_MEASUREMENT
#include <sstream>
#endif

#ifdef HAVE_IPP
#include <ippversion.h>
#include <ipps.h>
#endif

#ifdef HAVE_FFTW3
#include <fftw3.h>
#endif

#ifdef HAVE_VDSP
#include <Accelerate/Accelerate.h>
#endif

#ifdef HAVE_MEDIALIB
#include <mlib_signal.h>
#endif

#ifdef HAVE_OPENMAX
#include <omxSP.h>
#endif

#ifdef HAVE_SFFT
extern "C" {
#include <sfft.h>
}
#endif

#ifdef USE_KISSFFT
#include "kissfft/kiss_fftr.h"
#endif

#ifndef HAVE_IPP
#ifndef HAVE_FFTW3
#ifndef USE_KISSFFT
#ifndef USE_BUILTIN_FFT
#ifndef HAVE_VDSP
#ifndef HAVE_MEDIALIB
#ifndef HAVE_OPENMAX
#ifndef HAVE_SFFT
#error No FFT implementation selected!
#endif
#endif
#endif
#endif
#endif
#endif
#endif
#endif

#include <cmath>
#include <iostream>
#include <map>
#include <cstdio>
#include <cstdlib>
#include <vector>

#ifdef FFT_MEASUREMENT
#ifndef _WIN32
#include <unistd.h>
#endif
#endif

namespace RubberBand {

class FFTImpl
{
public:
    virtual ~FFTImpl() { }

    virtual FFT::Precisions getSupportedPrecisions() const = 0;

    virtual void initFloat() = 0;
    virtual void initDouble() = 0;

    virtual void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) = 0;
    virtual void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) = 0;
    virtual void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) = 0;
    virtual void forwardMagnitude(const double *R__ realIn, double *R__ magOut) = 0;

    virtual void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) = 0;
    virtual void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) = 0;
    virtual void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) = 0;
    virtual void forwardMagnitude(const float *R__ realIn, float *R__ magOut) = 0;

    virtual void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) = 0;
    virtual void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) = 0;
    virtual void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) = 0;
    virtual void inverseCepstral(const double *R__ magIn, double *R__ cepOut) = 0;

    virtual void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) = 0;
    virtual void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) = 0;
    virtual void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) = 0;
    virtual void inverseCepstral(const float *R__ magIn, float *R__ cepOut) = 0;
};    

namespace FFTs {

#ifdef HAVE_IPP

class D_IPP : public FFTImpl
{
public:
    D_IPP(int size) :
        m_size(size), m_fspec(0), m_dspec(0)
    { 
        for (int i = 0; ; ++i) {
            if (m_size & (1 << i)) {
                m_order = i;
                break;
            }
        }
    }

    ~D_IPP() {
        if (m_fspec) {
#if (IPP_VERSION_MAJOR >= 9)
            ippsFree(m_fspecbuf);
#else
            ippsFFTFree_R_32f(m_fspec);
#endif
            ippsFree(m_fbuf);
            ippsFree(m_fpacked);
            ippsFree(m_fspare);
        }
        if (m_dspec) {
#if (IPP_VERSION_MAJOR >= 9)
            ippsFree(m_dspecbuf);
#else
            ippsFFTFree_R_64f(m_dspec);
#endif
            ippsFree(m_dbuf);
            ippsFree(m_dpacked);
            ippsFree(m_dspare);
        }
    }

    FFT::Precisions
    getSupportedPrecisions() const {
        return FFT::SinglePrecision | FFT::DoublePrecision;
    }

    //!!! rv check

    void initFloat() {
        if (m_fspec) return;
#if (IPP_VERSION_MAJOR >= 9)
        int specSize, specBufferSize, bufferSize;
        ippsFFTGetSize_R_32f(m_order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                             &specSize, &specBufferSize, &bufferSize);
        m_fspecbuf = ippsMalloc_8u(specSize);
        Ipp8u *tmp = ippsMalloc_8u(specBufferSize);
        m_fbuf = ippsMalloc_8u(bufferSize);
        m_fpacked = ippsMalloc_32f(m_size + 2);
        m_fspare = ippsMalloc_32f(m_size / 2 + 1);
        ippsFFTInit_R_32f(&m_fspec,
                          m_order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                          m_fspecbuf, tmp);
        ippsFree(tmp);
#else
        int specSize, specBufferSize, bufferSize;
        ippsFFTGetSize_R_32f(m_order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                             &specSize, &specBufferSize, &bufferSize);
        m_fbuf = ippsMalloc_8u(bufferSize);
        m_fpacked = ippsMalloc_32f(m_size + 2);
        m_fspare = ippsMalloc_32f(m_size / 2 + 1);
        ippsFFTInitAlloc_R_32f(&m_fspec, m_order, IPP_FFT_NODIV_BY_ANY, 
                               ippAlgHintFast);
#endif
    }

    void initDouble() {
        if (m_dspec) return;
#if (IPP_VERSION_MAJOR >= 9)
        int specSize, specBufferSize, bufferSize;
        ippsFFTGetSize_R_64f(m_order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                             &specSize, &specBufferSize, &bufferSize);
        m_dspecbuf = ippsMalloc_8u(specSize);
        Ipp8u *tmp = ippsMalloc_8u(specBufferSize);
        m_dbuf = ippsMalloc_8u(bufferSize);
        m_dpacked = ippsMalloc_64f(m_size + 2);
        m_dspare = ippsMalloc_64f(m_size / 2 + 1);
        ippsFFTInit_R_64f(&m_dspec,
                          m_order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                          m_dspecbuf, tmp);
        ippsFree(tmp);
#else
        int specSize, specBufferSize, bufferSize;
        ippsFFTGetSize_R_64f(m_order, IPP_FFT_NODIV_BY_ANY, ippAlgHintFast,
                             &specSize, &specBufferSize, &bufferSize);
        m_dbuf = ippsMalloc_8u(bufferSize);
        m_dpacked = ippsMalloc_64f(m_size + 2);
        m_dspare = ippsMalloc_64f(m_size / 2 + 1);
        ippsFFTInitAlloc_R_64f(&m_dspec, m_order, IPP_FFT_NODIV_BY_ANY, 
                               ippAlgHintFast);
#endif
    }

    void packFloat(const float *R__ re, const float *R__ im) {
        Profiler profiler("D_IPP::packFloat");
        int index = 0;
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_fpacked[index++] = re[i];
            index++;
        }
        index = 0;
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_fpacked[index++] = im[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_fpacked[index++] = 0.f;
            }
        }
    }

    void packDouble(const double *R__ re, const double *R__ im) {
        Profiler profiler("D_IPP::packDouble");
        int index = 0;
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_dpacked[index++] = re[i];
            index++;
        }
        index = 0;
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_dpacked[index++] = im[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_dpacked[index++] = 0.0;
            }
        }
    }

    void unpackFloat(float *re, float *R__ im) { // re may be equal to m_fpacked
        Profiler profiler("D_IPP::unpackFloat");
        int index = 0;
        const int hs = m_size/2;
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                im[i] = m_fpacked[index++];
            }
        }
        index = 0;
        for (int i = 0; i <= hs; ++i) {
            re[i] = m_fpacked[index++];
            index++;
        }
    }        

    void unpackDouble(double *re, double *R__ im) { // re may be equal to m_dpacked
        Profiler profiler("D_IPP::unpackDouble");
        int index = 0;
        const int hs = m_size/2;
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                im[i] = m_dpacked[index++];
            }
        }
        index = 0;
        for (int i = 0; i <= hs; ++i) {
            re[i] = m_dpacked[index++];
            index++;
        }
    }        

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        Profiler profiler("D_IPP::forward [d]");
        if (!m_dspec) initDouble();
        ippsFFTFwd_RToCCS_64f(realIn, m_dpacked, m_dspec, m_dbuf);
        unpackDouble(realOut, imagOut);
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        Profiler profiler("D_IPP::forwardInterleaved [d]");
        if (!m_dspec) initDouble();
        ippsFFTFwd_RToCCS_64f(realIn, complexOut, m_dspec, m_dbuf);
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        Profiler profiler("D_IPP::forwardPolar [d]");
        if (!m_dspec) initDouble();
        ippsFFTFwd_RToCCS_64f(realIn, m_dpacked, m_dspec, m_dbuf);
        unpackDouble(m_dpacked, m_dspare);
        Profiler profiler2("D_IPP::forwardPolar [d] conv");
        ippsCartToPolar_64f(m_dpacked, m_dspare, magOut, phaseOut, m_size/2+1);
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        Profiler profiler("D_IPP::forwardMagnitude [d]");
        if (!m_dspec) initDouble();
        ippsFFTFwd_RToCCS_64f(realIn, m_dpacked, m_dspec, m_dbuf);
        unpackDouble(m_dpacked, m_dspare);
        ippsMagnitude_64f(m_dpacked, m_dspare, magOut, m_size/2+1);
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        Profiler profiler("D_IPP::forward [f]");
        if (!m_fspec) initFloat();
        ippsFFTFwd_RToCCS_32f(realIn, m_fpacked, m_fspec, m_fbuf);
        unpackFloat(realOut, imagOut);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        Profiler profiler("D_IPP::forwardInterleaved [f]");
        if (!m_fspec) initFloat();
        ippsFFTFwd_RToCCS_32f(realIn, complexOut, m_fspec, m_fbuf);
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        Profiler profiler("D_IPP::forwardPolar [f]");
        if (!m_fspec) initFloat();
        ippsFFTFwd_RToCCS_32f(realIn, m_fpacked, m_fspec, m_fbuf);
        unpackFloat(m_fpacked, m_fspare);
        Profiler profiler2("D_IPP::forwardPolar [f] conv");
        ippsCartToPolar_32f(m_fpacked, m_fspare, magOut, phaseOut, m_size/2+1);
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        Profiler profiler("D_IPP::forwardMagnitude [f]");
        if (!m_fspec) initFloat();
        ippsFFTFwd_RToCCS_32f(realIn, m_fpacked, m_fspec, m_fbuf);
        unpackFloat(m_fpacked, m_fspare);
        ippsMagnitude_32f(m_fpacked, m_fspare, magOut, m_size/2+1);
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        Profiler profiler("D_IPP::inverse [d]");
        if (!m_dspec) initDouble();
        packDouble(realIn, imagIn);
        ippsFFTInv_CCSToR_64f(m_dpacked, realOut, m_dspec, m_dbuf);
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        Profiler profiler("D_IPP::inverse [d]");
        if (!m_dspec) initDouble();
        ippsFFTInv_CCSToR_64f(complexIn, realOut, m_dspec, m_dbuf);
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        Profiler profiler("D_IPP::inversePolar [d]");
        if (!m_dspec) initDouble();
        ippsPolarToCart_64f(magIn, phaseIn, realOut, m_dspare, m_size/2+1);
        Profiler profiler2("D_IPP::inversePolar [d] postconv");
        packDouble(realOut, m_dspare); // to m_dpacked
        ippsFFTInv_CCSToR_64f(m_dpacked, realOut, m_dspec, m_dbuf);
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        Profiler profiler("D_IPP::inverseCepstral [d]");
        if (!m_dspec) initDouble();
        const int hs1 = m_size/2 + 1;
        ippsCopy_64f(magIn, m_dspare, hs1);
        ippsAddC_64f_I(0.000001, m_dspare, hs1);
        ippsLn_64f_I(m_dspare, hs1);
        packDouble(m_dspare, 0);
        ippsFFTInv_CCSToR_64f(m_dpacked, cepOut, m_dspec, m_dbuf);
    }
    
    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        Profiler profiler("D_IPP::inverse [f]");
        if (!m_fspec) initFloat();
        packFloat(realIn, imagIn);
        ippsFFTInv_CCSToR_32f(m_fpacked, realOut, m_fspec, m_fbuf);
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        Profiler profiler("D_IPP::inverse [f]");
        if (!m_fspec) initFloat();
        ippsFFTInv_CCSToR_32f(complexIn, realOut, m_fspec, m_fbuf);
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        Profiler profiler("D_IPP::inversePolar [f]");
        if (!m_fspec) initFloat();
        ippsPolarToCart_32f(magIn, phaseIn, realOut, m_fspare, m_size/2+1);
        Profiler profiler2("D_IPP::inversePolar [f] postconv");
        packFloat(realOut, m_fspare); // to m_fpacked
        ippsFFTInv_CCSToR_32f(m_fpacked, realOut, m_fspec, m_fbuf);
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        Profiler profiler("D_IPP::inverseCepstral [f]");
        if (!m_fspec) initFloat();
        const int hs1 = m_size/2 + 1;
        ippsCopy_32f(magIn, m_fspare, hs1);
        ippsAddC_32f_I(0.000001f, m_fspare, hs1);
        ippsLn_32f_I(m_fspare, hs1);
        packFloat(m_fspare, 0);
        ippsFFTInv_CCSToR_32f(m_fpacked, cepOut, m_fspec, m_fbuf);
    }

private:
    const int m_size;
    int m_order;
    IppsFFTSpec_R_32f *m_fspec;
    IppsFFTSpec_R_64f *m_dspec;
    Ipp8u *m_fspecbuf;
    Ipp8u *m_dspecbuf;
    Ipp8u *m_fbuf;
    Ipp8u *m_dbuf;
    float *m_fpacked;
    float *m_fspare;
    double *m_dpacked;
    double *m_dspare;
};

#endif /* HAVE_IPP */

#ifdef HAVE_VDSP

class D_VDSP : public FFTImpl
{
public:
    D_VDSP(int size) :
        m_size(size), m_fspec(0), m_dspec(0),
        m_fpacked(0), m_fspare(0),
        m_dpacked(0), m_dspare(0)
    { 
        for (int i = 0; ; ++i) {
            if (m_size & (1 << i)) {
                m_order = i;
                break;
            }
        }
    }

    ~D_VDSP() {
        if (m_fspec) {
            vDSP_destroy_fftsetup(m_fspec);
            deallocate(m_fspare);
            deallocate(m_fspare2);
            deallocate(m_fbuf->realp);
            deallocate(m_fbuf->imagp);
            delete m_fbuf;
            deallocate(m_fpacked->realp);
            deallocate(m_fpacked->imagp);
            delete m_fpacked;
        }
        if (m_dspec) {
            vDSP_destroy_fftsetupD(m_dspec);
            deallocate(m_dspare);
            deallocate(m_dspare2);
            deallocate(m_dbuf->realp);
            deallocate(m_dbuf->imagp);
            delete m_dbuf;
            deallocate(m_dpacked->realp);
            deallocate(m_dpacked->imagp);
            delete m_dpacked;
        }
    }

    FFT::Precisions
    getSupportedPrecisions() const {
        return FFT::SinglePrecision | FFT::DoublePrecision;
    }

    //!!! rv check

    void initFloat() {
        if (m_fspec) return;
        m_fspec = vDSP_create_fftsetup(m_order, FFT_RADIX2);
        m_fbuf = new DSPSplitComplex;
        //!!! "If possible, tempBuffer->realp and tempBuffer->imagp should be 32-byte aligned for best performance."
        m_fbuf->realp = allocate<float>(m_size);
        m_fbuf->imagp = allocate<float>(m_size);
        m_fpacked = new DSPSplitComplex;
        m_fpacked->realp = allocate<float>(m_size / 2 + 1);
        m_fpacked->imagp = allocate<float>(m_size / 2 + 1);
        m_fspare = allocate<float>(m_size + 2);
        m_fspare2 = allocate<float>(m_size + 2);
    }

    void initDouble() {
        if (m_dspec) return;
        m_dspec = vDSP_create_fftsetupD(m_order, FFT_RADIX2);
        m_dbuf = new DSPDoubleSplitComplex;
        //!!! "If possible, tempBuffer->realp and tempBuffer->imagp should be 32-byte aligned for best performance."
        m_dbuf->realp = allocate<double>(m_size);
        m_dbuf->imagp = allocate<double>(m_size);
        m_dpacked = new DSPDoubleSplitComplex;
        m_dpacked->realp = allocate<double>(m_size / 2 + 1);
        m_dpacked->imagp = allocate<double>(m_size / 2 + 1);
        m_dspare = allocate<double>(m_size + 2);
        m_dspare2 = allocate<double>(m_size + 2);
    }

    void packReal(const float *R__ const re) {
        // Pack input for forward transform 
        vDSP_ctoz((DSPComplex *)re, 2, m_fpacked, 1, m_size/2);
    }
    void packComplex(const float *R__ const re, const float *R__ const im) {
        // Pack input for inverse transform 
        if (re) v_copy(m_fpacked->realp, re, m_size/2 + 1);
        else v_zero(m_fpacked->realp, m_size/2 + 1);
        if (im) v_copy(m_fpacked->imagp, im, m_size/2 + 1);
        else v_zero(m_fpacked->imagp, m_size/2 + 1);
        fnyq();
    }

    void unpackReal(float *R__ const re) {
        // Unpack output for inverse transform
        vDSP_ztoc(m_fpacked, 1, (DSPComplex *)re, 2, m_size/2);
    }
    void unpackComplex(float *R__ const re, float *R__ const im) {
        // Unpack output for forward transform
        // vDSP forward FFTs are scaled 2x (for some reason)
        float two = 2.f;
        vDSP_vsdiv(m_fpacked->realp, 1, &two, re, 1, m_size/2 + 1);
        vDSP_vsdiv(m_fpacked->imagp, 1, &two, im, 1, m_size/2 + 1);
    }
    void unpackComplex(float *R__ const cplx) {
        // Unpack output for forward transform
        // vDSP forward FFTs are scaled 2x (for some reason)
        const int hs1 = m_size/2 + 1;
        for (int i = 0; i < hs1; ++i) {
            cplx[i*2] = m_fpacked->realp[i] / 2.f;
            cplx[i*2+1] = m_fpacked->imagp[i] / 2.f;
        }
    }

    void packReal(const double *R__ const re) {
        // Pack input for forward transform
        vDSP_ctozD((DSPDoubleComplex *)re, 2, m_dpacked, 1, m_size/2);
    }
    void packComplex(const double *R__ const re, const double *R__ const im) {
        // Pack input for inverse transform
        if (re) v_copy(m_dpacked->realp, re, m_size/2 + 1);
        else v_zero(m_dpacked->realp, m_size/2 + 1);
        if (im) v_copy(m_dpacked->imagp, im, m_size/2 + 1);
        else v_zero(m_dpacked->imagp, m_size/2 + 1);
        dnyq();
    }

    void unpackReal(double *R__ const re) {
        // Unpack output for inverse transform
        vDSP_ztocD(m_dpacked, 1, (DSPDoubleComplex *)re, 2, m_size/2);
    }
    void unpackComplex(double *R__ const re, double *R__ const im) {
        // Unpack output for forward transform
        // vDSP forward FFTs are scaled 2x (for some reason)
        double two = 2.0;
        vDSP_vsdivD(m_dpacked->realp, 1, &two, re, 1, m_size/2 + 1);
        vDSP_vsdivD(m_dpacked->imagp, 1, &two, im, 1, m_size/2 + 1);
    }
    void unpackComplex(double *R__ const cplx) {
        // Unpack output for forward transform
        // vDSP forward FFTs are scaled 2x (for some reason)
        const int hs1 = m_size/2 + 1;
        for (int i = 0; i < hs1; ++i) {
            cplx[i*2] = m_dpacked->realp[i] / 2.0;
            cplx[i*2+1] = m_dpacked->imagp[i] / 2.0;
        }
    }

    void fdenyq() {
        // for fft result in packed form, unpack the DC and Nyquist bins
        const int hs = m_size/2;
        m_fpacked->realp[hs] = m_fpacked->imagp[0];
        m_fpacked->imagp[hs] = 0.f;
        m_fpacked->imagp[0] = 0.f;
    }
    void ddenyq() {
        // for fft result in packed form, unpack the DC and Nyquist bins
        const int hs = m_size/2;
        m_dpacked->realp[hs] = m_dpacked->imagp[0];
        m_dpacked->imagp[hs] = 0.;
        m_dpacked->imagp[0] = 0.;
    }

    void fnyq() {
        // for ifft input in packed form, pack the DC and Nyquist bins
        const int hs = m_size/2;
        m_fpacked->imagp[0] = m_fpacked->realp[hs];
        m_fpacked->realp[hs] = 0.f;
        m_fpacked->imagp[hs] = 0.f;
    }
    void dnyq() {
        // for ifft input in packed form, pack the DC and Nyquist bins
        const int hs = m_size/2;
        m_dpacked->imagp[0] = m_dpacked->realp[hs];
        m_dpacked->realp[hs] = 0.;
        m_dpacked->imagp[hs] = 0.;
    }

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        Profiler profiler("D_VDSP::forward [d]");
        if (!m_dspec) initDouble();
        packReal(realIn);
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_FORWARD);
        ddenyq();
        unpackComplex(realOut, imagOut);
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        Profiler profiler("D_VDSP::forward [d]");
        if (!m_dspec) initDouble();
        packReal(realIn);
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_FORWARD);
        ddenyq();
        unpackComplex(complexOut);
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        Profiler profiler("D_VDSP::forwardPolar [d]");
        if (!m_dspec) initDouble();
        const int hs1 = m_size/2+1;
        packReal(realIn);
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_FORWARD);
        ddenyq();
        // vDSP forward FFTs are scaled 2x (for some reason)
        for (int i = 0; i < hs1; ++i) m_dpacked->realp[i] /= 2.0;
        for (int i = 0; i < hs1; ++i) m_dpacked->imagp[i] /= 2.0;
        v_cartesian_to_polar(magOut, phaseOut,
                             m_dpacked->realp, m_dpacked->imagp, hs1);
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        Profiler profiler("D_VDSP::forwardMagnitude [d]");
        if (!m_dspec) initDouble();
        packReal(realIn);
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_FORWARD);
        ddenyq();
        const int hs1 = m_size/2+1;
        vDSP_zvmagsD(m_dpacked, 1, m_dspare, 1, hs1);
        vvsqrt(m_dspare2, m_dspare, &hs1);
        // vDSP forward FFTs are scaled 2x (for some reason)
        double two = 2.0;
        vDSP_vsdivD(m_dspare2, 1, &two, magOut, 1, hs1);
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        Profiler profiler("D_VDSP::forward [f]");
        if (!m_fspec) initFloat();
        packReal(realIn);
        vDSP_fft_zript(m_fspec, m_fpacked, 1, m_fbuf, m_order, FFT_FORWARD);
        fdenyq();
        unpackComplex(realOut, imagOut);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        Profiler profiler("D_VDSP::forward [f]");
        if (!m_fspec) initFloat();
        packReal(realIn);
        vDSP_fft_zript(m_fspec, m_fpacked, 1, m_fbuf, m_order, FFT_FORWARD);
        fdenyq();
        unpackComplex(complexOut);
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        Profiler profiler("D_VDSP::forwardPolar [f]");
        if (!m_fspec) initFloat();
        const int hs1 = m_size/2+1;
        packReal(realIn);
        vDSP_fft_zript(m_fspec, m_fpacked, 1, m_fbuf, m_order, FFT_FORWARD);
        fdenyq();
        // vDSP forward FFTs are scaled 2x (for some reason)
        for (int i = 0; i < hs1; ++i) m_fpacked->realp[i] /= 2.f;
        for (int i = 0; i < hs1; ++i) m_fpacked->imagp[i] /= 2.f;
        v_cartesian_to_polar(magOut, phaseOut,
                             m_fpacked->realp, m_fpacked->imagp, hs1);
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        Profiler profiler("D_VDSP::forwardMagnitude [f]");
        if (!m_fspec) initFloat();
        packReal(realIn);
        vDSP_fft_zript(m_fspec, m_fpacked, 1, m_fbuf, m_order, FFT_FORWARD);
        fdenyq();
        const int hs1 = m_size/2 + 1;
        vDSP_zvmags(m_fpacked, 1, m_fspare, 1, hs1);
        vvsqrtf(m_fspare2, m_fspare, &hs1);
        // vDSP forward FFTs are scaled 2x (for some reason)
        float two = 2.f;
        vDSP_vsdiv(m_fspare2, 1, &two, magOut, 1, hs1);
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        Profiler profiler("D_VDSP::inverse [d]");
        if (!m_dspec) initDouble();
        packComplex(realIn, imagIn);
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_INVERSE);
        unpackReal(realOut);
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        Profiler profiler("D_VDSP::inverseInterleaved [d]");
        if (!m_dspec) initDouble();
        double *d[2] = { m_dpacked->realp, m_dpacked->imagp };
        v_deinterleave(d, complexIn, 2, m_size/2 + 1);
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_INVERSE);
        unpackReal(realOut);
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        Profiler profiler("D_VDSP::inversePolar [d]");
        if (!m_dspec) initDouble();
        const int hs1 = m_size/2+1;
        vvsincos(m_dpacked->imagp, m_dpacked->realp, phaseIn, &hs1);
        double *const rp = m_dpacked->realp;
        double *const ip = m_dpacked->imagp;
        for (int i = 0; i < hs1; ++i) rp[i] *= magIn[i];
        for (int i = 0; i < hs1; ++i) ip[i] *= magIn[i];
        dnyq();
        vDSP_fft_zriptD(m_dspec, m_dpacked, 1, m_dbuf, m_order, FFT_INVERSE);
        unpackReal(realOut);
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        Profiler profiler("D_VDSP::inverseCepstral [d]");
        if (!m_dspec) initDouble();
        const int hs1 = m_size/2 + 1;
        v_copy(m_dspare, magIn, hs1);
        for (int i = 0; i < hs1; ++i) m_dspare[i] += 0.000001;
        vvlog(m_dspare2, m_dspare, &hs1);
        inverse(m_dspare2, 0, cepOut);
    }
    
    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        Profiler profiler("D_VDSP::inverse [f]");
        if (!m_fspec) initFloat();
        packComplex(realIn, imagIn);
        vDSP_fft_zript(m_fspec, m_fpacked, 1, m_fbuf, m_order, FFT_INVERSE);
        unpackReal(realOut);
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        Profiler profiler("D_VDSP::inverseInterleaved [f]");
        if (!m_fspec) initFloat();
        float *f[2] = { m_fpacked->realp, m_fpacked->imagp };
        v_deinterleave(f, complexIn, 2, m_size/2 + 1);
        vDSP_fft_zript(m_fspec, m_fpacked, 1, m_fbuf, m_order, FFT_INVERSE);
        unpackReal(realOut);
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        Profiler profiler("D_VDSP::inversePolar [f]");
        if (!m_fspec) initFloat();

        const int hs1 = m_size/2+1;
        vvsincosf(m_fpacked->imagp, m_fpacked->realp, phaseIn, &hs1);
        float *const rp = m_fpacked->realp;
        float *const ip = m_fpacked->imagp;
        for (int i = 0; i < hs1; ++i) rp[i] *= magIn[i];
        for (int i = 0; i < hs1; ++i) ip[i] *= magIn[i];
        fnyq();
        vDSP_fft_zript(m_fspec, m_fpacked, 1, m_fbuf, m_order, FFT_INVERSE);
        unpackReal(realOut);
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        Profiler profiler("D_VDSP::inverseCepstral [f]");
        if (!m_fspec) initFloat();
        const int hs1 = m_size/2 + 1;
        v_copy(m_fspare, magIn, hs1);
        for (int i = 0; i < hs1; ++i) m_fspare[i] += 0.000001f;
        vvlogf(m_fspare2, m_fspare, &hs1);
        inverse(m_fspare2, 0, cepOut);
    }

private:
    const int m_size;
    int m_order;
    FFTSetup m_fspec;
    FFTSetupD m_dspec;
    DSPSplitComplex *m_fbuf;
    DSPDoubleSplitComplex *m_dbuf;
    DSPSplitComplex *m_fpacked;
    float *m_fspare;
    float *m_fspare2;
    DSPDoubleSplitComplex *m_dpacked;
    double *m_dspare;
    double *m_dspare2;
};

#endif /* HAVE_VDSP */

#ifdef HAVE_MEDIALIB

class D_MEDIALIB : public FFTImpl
{
public:
    D_MEDIALIB(int size) :
        m_size(size),
        m_dpacked(0), m_fpacked(0)
    { 
        for (int i = 0; ; ++i) {
            if (m_size & (1 << i)) {
                m_order = i;
                break;
            }
        }
    }

    ~D_MEDIALIB() {
        if (m_dpacked) {
            deallocate(m_dpacked);
        }
        if (m_fpacked) {
            deallocate(m_fpacked);
        }
    }

    FFT::Precisions
    getSupportedPrecisions() const {
        return FFT::SinglePrecision | FFT::DoublePrecision;
    }

    //!!! rv check

    void initFloat() {
        m_fpacked = allocate<float>(m_size*2);
    }

    void initDouble() {
        m_dpacked = allocate<double>(m_size*2);
    }

    void packFloatConjugates() {
        const int hs = m_size / 2;
        for (int i = 1; i <= hs; ++i) {
            m_fpacked[(m_size-i)*2] = m_fpacked[2*i];
            m_fpacked[(m_size-i)*2 + 1] = -m_fpacked[2*i + 1];
        }
    }

    void packDoubleConjugates() {
        const int hs = m_size / 2;
        for (int i = 1; i <= hs; ++i) {
            m_dpacked[(m_size-i)*2] = m_dpacked[2*i];
            m_dpacked[(m_size-i)*2 + 1] = -m_dpacked[2*i + 1];
        }
    }

    void packFloat(const float *R__ re, const float *R__ im) {
        int index = 0;
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_fpacked[index++] = re[i];
            index++;
        }
        index = 0;
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_fpacked[index++] = im[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_fpacked[index++] = 0.f;
            }
        }
        packFloatConjugates();
    }

    void packDouble(const double *R__ re, const double *R__ im) {
        int index = 0;
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_dpacked[index++] = re[i];
            index++;
        }
        index = 0;
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_dpacked[index++] = im[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_dpacked[index++] = 0.0;
            }
        }
        packDoubleConjugates();
    }

    void unpackFloat(float *re, float *R__ im) { // re may be equal to m_fpacked
        int index = 0;
        const int hs = m_size/2;
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                im[i] = m_fpacked[index++];
            }
        }
        index = 0;
        for (int i = 0; i <= hs; ++i) {
            re[i] = m_fpacked[index++];
            index++;
        }
    }        

    void unpackDouble(double *re, double *R__ im) { // re may be equal to m_dpacked
        int index = 0;
        const int hs = m_size/2;
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                im[i] = m_dpacked[index++];
            }
        }
        index = 0;
        for (int i = 0; i <= hs; ++i) {
            re[i] = m_dpacked[index++];
            index++;
        }
    }

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        Profiler profiler("D_MEDIALIB::forward [d]");
        if (!m_dpacked) initDouble();
        mlib_SignalFFT_1_D64C_D64(m_dpacked, realIn, m_order);
        unpackDouble(realOut, imagOut);
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        Profiler profiler("D_MEDIALIB::forwardInterleaved [d]");
        if (!m_dpacked) initDouble();
        // mlib FFT gives the whole redundant complex result
        mlib_SignalFFT_1_D64C_D64(m_dpacked, realIn, m_order);
        v_copy(complexOut, m_dpacked, m_size + 2);
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        Profiler profiler("D_MEDIALIB::forwardPolar [d]");
        if (!m_dpacked) initDouble();
        mlib_SignalFFT_1_D64C_D64(m_dpacked, realIn, m_order);
        const int hs = m_size/2;
        int index = 0;
        for (int i = 0; i <= hs; ++i) {
            int reali = index;
            ++index;
            magOut[i] = sqrt(m_dpacked[reali] * m_dpacked[reali] +
                             m_dpacked[index] * m_dpacked[index]);
            phaseOut[i] = atan2(m_dpacked[index], m_dpacked[reali]) ;
            ++index;
        }
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        Profiler profiler("D_MEDIALIB::forwardMagnitude [d]");
        if (!m_dpacked) initDouble();
        mlib_SignalFFT_1_D64C_D64(m_dpacked, realIn, m_order);
        const int hs = m_size/2;
        int index = 0;
        for (int i = 0; i <= hs; ++i) {
            int reali = index;
            ++index;
            magOut[i] = sqrt(m_dpacked[reali] * m_dpacked[reali] +
                             m_dpacked[index] * m_dpacked[index]);
            ++index;
        }
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        Profiler profiler("D_MEDIALIB::forward [f]");
        if (!m_fpacked) initFloat();
        mlib_SignalFFT_1_F32C_F32(m_fpacked, realIn, m_order);
        unpackFloat(realOut, imagOut);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        Profiler profiler("D_MEDIALIB::forwardInterleaved [f]");
        if (!m_fpacked) initFloat();
        // mlib FFT gives the whole redundant complex result
        mlib_SignalFFT_1_F32C_F32(m_fpacked, realIn, m_order);
        v_copy(complexOut, m_fpacked, m_size + 2);
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        Profiler profiler("D_MEDIALIB::forwardPolar [f]");
        if (!m_fpacked) initFloat();
        mlib_SignalFFT_1_F32C_F32(m_fpacked, realIn, m_order);
        const int hs = m_size/2;
        int index = 0;
        for (int i = 0; i <= hs; ++i) {
            int reali = index;
            ++index;
            magOut[i] = sqrtf(m_fpacked[reali] * m_fpacked[reali] +
                              m_fpacked[index] * m_fpacked[index]);
            phaseOut[i] = atan2f(m_fpacked[index], m_fpacked[reali]);
            ++index;
        }
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        Profiler profiler("D_MEDIALIB::forwardMagnitude [f]");
        if (!m_fpacked) initFloat();
        mlib_SignalFFT_1_F32C_F32(m_fpacked, realIn, m_order);
        const int hs = m_size/2;
        int index = 0;
        for (int i = 0; i <= hs; ++i) {
            int reali = index;
            ++index;
            magOut[i] = sqrtf(m_fpacked[reali] * m_fpacked[reali] +
                              m_fpacked[index] * m_fpacked[index]);
            ++index;
        }
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        Profiler profiler("D_MEDIALIB::inverse [d]");
        if (!m_dpacked) initDouble();
        packDouble(realIn, imagIn);
        mlib_SignalIFFT_2_D64_D64C(realOut, m_dpacked, m_order);
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        Profiler profiler("D_MEDIALIB::inverseInterleaved [d]");
        if (!m_dpacked) initDouble();
        v_copy(m_dpacked, complexIn, m_size + 2);
        packDoubleConjugates();
        mlib_SignalIFFT_2_D64_D64C(realOut, m_dpacked, m_order);
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        Profiler profiler("D_MEDIALIB::inversePolar [d]");
        if (!m_dpacked) initDouble();
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            double real = magIn[i] * cos(phaseIn[i]);
            double imag = magIn[i] * sin(phaseIn[i]);
            m_dpacked[i*2] = real;
            m_dpacked[i*2 + 1] = imag;
        }
        packDoubleConjugates();
        mlib_SignalIFFT_2_D64_D64C(realOut, m_dpacked, m_order);
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        Profiler profiler("D_MEDIALIB::inverseCepstral [d]");
        if (!m_dpacked) initDouble();
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_dpacked[i*2] = log(magIn[i] + 0.000001);
            m_dpacked[i*2 + 1] = 0.0;
        }
        packDoubleConjugates();
        mlib_SignalIFFT_2_D64_D64C(cepOut, m_dpacked, m_order);
    }
    
    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        Profiler profiler("D_MEDIALIB::inverse [f]");
        if (!m_fpacked) initFloat();
        packFloat(realIn, imagIn);
        mlib_SignalIFFT_2_F32_F32C(realOut, m_fpacked, m_order);
    }
    
    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        Profiler profiler("D_MEDIALIB::inverseInterleaved [f]");
        if (!m_fpacked) initFloat();
        v_convert(m_fpacked, complexIn, m_size + 2);
        packFloatConjugates();
        mlib_SignalIFFT_2_F32_F32C(realOut, m_fpacked, m_order);
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        Profiler profiler("D_MEDIALIB::inversePolar [f]");
        if (!m_fpacked) initFloat();
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            double real = magIn[i] * cos(phaseIn[i]);
            double imag = magIn[i] * sin(phaseIn[i]);
            m_fpacked[i*2] = real;
            m_fpacked[i*2 + 1] = imag;
        }
        packFloatConjugates();
        mlib_SignalIFFT_2_F32_F32C(realOut, m_fpacked, m_order);
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        Profiler profiler("D_MEDIALIB::inverseCepstral [f]");
        if (!m_fpacked) initFloat();
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_fpacked[i*2] = logf(magIn[i] + 0.000001);
            m_fpacked[i*2 + 1] = 0.f;
        }
        packFloatConjugates();
        mlib_SignalIFFT_2_F32_F32C(cepOut, m_fpacked, m_order);
    }

private:
    const int m_size;
    int m_order;
    double *m_dpacked;
    float *m_fpacked;
};

#endif /* HAVE_MEDIALIB */

#ifdef HAVE_OPENMAX

class D_OPENMAX : public FFTImpl
{
    // Convert a signed 32-bit integer to a float in the range [-1,1)
    static inline float i2f(OMX_S32 i)
    {
        return float(i) / float(OMX_MAX_S32);
    }

    // Convert a signed 32-bit integer to a double in the range [-1,1)
    static inline double i2d(OMX_S32 i)
    {
        return double(i) / double(OMX_MAX_S32);
    }

    // Convert a float in the range [-1,1) to a signed 32-bit integer
    static inline OMX_S32 f2i(float f)
    {
        return OMX_S32(f * OMX_MAX_S32);
    }

    // Convert a double in the range [-1,1) to a signed 32-bit integer
    static inline OMX_S32 d2i(double d)
    {
        return OMX_S32(d * OMX_MAX_S32);
    }

public:
    D_OPENMAX(int size) :
        m_size(size),
        m_packed(0)
    { 
        for (int i = 0; ; ++i) {
            if (m_size & (1 << i)) {
                m_order = i;
                break;
            }
        }
    }

    ~D_OPENMAX() {
        if (m_packed) {
            deallocate(m_packed);
            deallocate(m_buf);
            deallocate(m_fbuf);
            deallocate(m_spec);
        }
    }

    FFT::Precisions
    getSupportedPrecisions() const {
        return FFT::SinglePrecision;
    }

    //!!! rv check

    // The OpenMAX implementation uses a fixed-point representation in
    // 32-bit signed integers, with a downward scaling factor (0-32
    // bits) supplied as an argument to the FFT function.

    void initFloat() {
        initDouble();
    }

    void initDouble() {
        if (!m_packed) {
            m_buf = allocate<OMX_S32>(m_size);
            m_packed = allocate<OMX_S32>(m_size*2 + 2);
            m_fbuf = allocate<float>(m_size*2 + 2);
            OMX_INT sz = 0;
            omxSP_FFTGetBufSize_R_S32(m_order, &sz);
            m_spec = (OMXFFTSpec_R_S32 *)allocate<char>(sz);
            omxSP_FFTInit_R_S32(m_spec, m_order);
        }
    }

    void packFloat(const float *R__ re) {
        // prepare fixed point input for forward transform
        for (int i = 0; i < m_size; ++i) {
            m_buf[i] = f2i(re[i]);
        }
    }

    void packDouble(const double *R__ re) {
        // prepare fixed point input for forward transform
        for (int i = 0; i < m_size; ++i) {
            m_buf[i] = d2i(re[i]);
        }
    }

    void unpackFloat(float *R__ re, float *R__ im) {
        // convert fixed point output for forward transform
        int index = 0;
        const int hs = m_size/2;
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                im[i] = i2f(m_packed[index++]);
            }
            v_scale(im, m_size, hs + 1);
        }
        index = 0;
        for (int i = 0; i <= hs; ++i) {
            re[i] = i2f(m_packed[index++]);
            index++;
        }
        v_scale(re, m_size, hs + 1);
    }        

    void unpackDouble(double *R__ re, double *R__ im) {
        // convert fixed point output for forward transform
        int index = 0;
        const int hs = m_size/2;
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                im[i] = i2d(m_packed[index++]);
            }
            v_scale(im, m_size, hs + 1);
        }
        index = 0;
        for (int i = 0; i <= hs; ++i) {
            re[i] = i2d(m_packed[index++]);
            index++;
        }
        v_scale(re, m_size, hs + 1);
    }

    void unpackFloatInterleaved(float *R__ cplx) {
        // convert fixed point output for forward transform
        for (int i = 0; i < m_size + 2; ++i) {
            cplx[i] = i2f(m_packed[i]);
        }            
        v_scale(cplx, m_size, m_size + 2);
    }

    void unpackDoubleInterleaved(double *R__ cplx) {
        // convert fixed point output for forward transform
        for (int i = 0; i < m_size + 2; ++i) {
            cplx[i] = i2d(m_packed[i]);
        }            
        v_scale(cplx, m_size, m_size + 2);
    }

    void packFloat(const float *R__ re, const float *R__ im) {
        // prepare fixed point input for inverse transform
        int index = 0;
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_packed[index++] = f2i(re[i]);
            index++;
        }
        index = 0;
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_packed[index++] = f2i(im[i]);
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_packed[index++] = 0;
            }
        }
    }

    void packDouble(const double *R__ re, const double *R__ im) {
        // prepare fixed point input for inverse transform
        int index = 0;
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_packed[index++] = d2i(re[i]);
            index++;
        }
        index = 0;
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_packed[index++] = d2i(im[i]);
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                index++;
                m_packed[index++] = 0;
            }
        }
    }

    void convertFloat(const float *R__ f) {
        // convert interleaved input for inverse interleaved transform
        const int n = m_size + 2;
        for (int i = 0; i < n; ++i) {
            m_packed[i] = f2i(f[i]);
        }
    }        

    void convertDouble(const double *R__ d) {
        // convert interleaved input for inverse interleaved transform
        const int n = m_size + 2;
        for (int i = 0; i < n; ++i) {
            m_packed[i] = d2i(d[i]);
        }
    }        

    void unpackFloat(float *R__ re) {
        // convert fixed point output for inverse transform
        for (int i = 0; i < m_size; ++i) {
            re[i] = i2f(m_buf[i]) * m_size;
        }
    }

    void unpackDouble(double *R__ re) {
        // convert fixed point output for inverse transform
        for (int i = 0; i < m_size; ++i) {
            re[i] = i2d(m_buf[i]) * m_size;
        }
    }

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        Profiler profiler("D_OPENMAX::forward [d]");
        if (!m_packed) initDouble();
        packDouble(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        unpackDouble(realOut, imagOut);
    }
    
    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        Profiler profiler("D_OPENMAX::forwardInterleaved [d]");
        if (!m_packed) initDouble();
        packDouble(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        unpackDoubleInterleaved(complexOut);
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        Profiler profiler("D_OPENMAX::forwardPolar [d]");
        if (!m_packed) initDouble();
        packDouble(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        unpackDouble(magOut, phaseOut); // temporarily
        // at this point we actually have real/imag in the mag/phase arrays
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            double real = magOut[i];
            double imag = phaseOut[i];
            c_magphase(magOut + i, phaseOut + i, real, imag);
        }
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        Profiler profiler("D_OPENMAX::forwardMagnitude [d]");
        if (!m_packed) initDouble();
        packDouble(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            int reali = i * 2;
            int imagi = reali + 1;
            double real = i2d(m_packed[reali]) * m_size;
            double imag = i2d(m_packed[imagi]) * m_size;
            magOut[i] = sqrt(real * real + imag * imag);
        }
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        Profiler profiler("D_OPENMAX::forward [f]");
        if (!m_packed) initFloat();
        packFloat(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        unpackFloat(realOut, imagOut);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        Profiler profiler("D_OPENMAX::forwardInterleaved [f]");
        if (!m_packed) initFloat();
        packFloat(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        unpackFloatInterleaved(complexOut);
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        Profiler profiler("D_OPENMAX::forwardPolar [f]");
        if (!m_packed) initFloat();

        packFloat(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        unpackFloat(magOut, phaseOut); // temporarily
        // at this point we actually have real/imag in the mag/phase arrays
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            float real = magOut[i];
            float imag = phaseOut[i];
            c_magphase(magOut + i, phaseOut + i, real, imag);
        }
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        Profiler profiler("D_OPENMAX::forwardMagnitude [f]");
        if (!m_packed) initFloat();
        packFloat(realIn);
        omxSP_FFTFwd_RToCCS_S32_Sfs(m_buf, m_packed, m_spec, m_order);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            int reali = i * 2;
            int imagi = reali + 1;
            float real = i2f(m_packed[reali]) * m_size;
            float imag = i2f(m_packed[imagi]) * m_size;
            magOut[i] = sqrtf(real * real + imag * imag);
        }
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        Profiler profiler("D_OPENMAX::inverse [d]");
        if (!m_packed) initDouble();
        packDouble(realIn, imagIn);
        omxSP_FFTInv_CCSToR_S32_Sfs(m_packed, m_buf, m_spec, 0);
        unpackDouble(realOut);
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        Profiler profiler("D_OPENMAX::inverseInterleaved [d]");
        if (!m_packed) initDouble();
        convertDouble(complexIn);
        omxSP_FFTInv_CCSToR_S32_Sfs(m_packed, m_buf, m_spec, 0);
        unpackDouble(realOut);
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        Profiler profiler("D_OPENMAX::inversePolar [d]");
        if (!m_packed) initDouble();
        int index = 0;
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            double real, imag;
            c_phasor(&real, &imag, phaseIn[i]);
            m_fbuf[index++] = float(real);
            m_fbuf[index++] = float(imag);
        }
        convertFloat(m_fbuf);
        omxSP_FFTInv_CCSToR_S32_Sfs(m_packed, m_buf, m_spec, 0);
        unpackDouble(realOut);
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        Profiler profiler("D_OPENMAX::inverseCepstral [d]");
        if (!m_packed) initDouble();
        //!!! implement
    }
    
    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        Profiler profiler("D_OPENMAX::inverse [f]");
        if (!m_packed) initFloat();
        packFloat(realIn, imagIn);
        omxSP_FFTInv_CCSToR_S32_Sfs(m_packed, m_buf, m_spec, 0);
        unpackFloat(realOut);
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        Profiler profiler("D_OPENMAX::inverse [f]");
        if (!m_packed) initFloat();
        convertFloat(complexIn);
        omxSP_FFTInv_CCSToR_S32_Sfs(m_packed, m_buf, m_spec, 0);
        unpackFloat(realOut);
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        Profiler profiler("D_OPENMAX::inversePolar [f]");
        if (!m_packed) initFloat();
        const int hs = m_size/2;
        v_polar_to_cartesian_interleaved(m_fbuf, magIn, phaseIn, hs+1);
        convertFloat(m_fbuf);
        omxSP_FFTInv_CCSToR_S32_Sfs(m_packed, m_buf, m_spec, 0);
        unpackFloat(realOut);
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        Profiler profiler("D_OPENMAX::inverseCepstral [f]");
        if (!m_packed) initFloat();
        //!!! implement
    }

private:
    const int m_size;
    int m_order;
    OMX_S32 *m_packed;
    OMX_S32 *m_buf;
    float *m_fbuf;
    OMXFFTSpec_R_S32 *m_spec;

};

#endif /* HAVE_OPENMAX */

#ifdef HAVE_FFTW3

/*
 Define FFTW_DOUBLE_ONLY to make all uses of FFTW functions be
 double-precision (so "float" FFTs are calculated by casting to
 doubles and using the double-precision FFTW function).

 Define FFTW_SINGLE_ONLY to make all uses of FFTW functions be
 single-precision (so "double" FFTs are calculated by casting to
 floats and using the single-precision FFTW function).

 Neither of these flags is desirable for either performance or
 precision. The main reason to define either flag is to avoid linking
 against both fftw3 and fftw3f libraries.
*/

//#define FFTW_DOUBLE_ONLY 1
//#define FFTW_SINGLE_ONLY 1

#if defined(FFTW_DOUBLE_ONLY) && defined(FFTW_SINGLE_ONLY)
// Can't meaningfully define both
#error Can only define one of FFTW_DOUBLE_ONLY and FFTW_SINGLE_ONLY
#endif

#if defined(FFTW_FLOAT_ONLY)
#warning FFTW_FLOAT_ONLY is deprecated, use FFTW_SINGLE_ONLY instead
#define FFTW_SINGLE_ONLY 1
#endif

#ifdef FFTW_DOUBLE_ONLY
#define fft_float_type double
#define fftwf_complex fftw_complex
#define fftwf_plan fftw_plan
#define fftwf_plan_dft_r2c_1d fftw_plan_dft_r2c_1d
#define fftwf_plan_dft_c2r_1d fftw_plan_dft_c2r_1d
#define fftwf_destroy_plan fftw_destroy_plan
#define fftwf_malloc fftw_malloc
#define fftwf_free fftw_free
#define fftwf_execute fftw_execute
#define atan2f atan2
#define sqrtf sqrt
#define cosf cos
#define sinf sin
#else
#define fft_float_type float
#endif /* FFTW_DOUBLE_ONLY */

#ifdef FFTW_SINGLE_ONLY
#define fft_double_type float
#define fftw_complex fftwf_complex
#define fftw_plan fftwf_plan
#define fftw_plan_dft_r2c_1d fftwf_plan_dft_r2c_1d
#define fftw_plan_dft_c2r_1d fftwf_plan_dft_c2r_1d
#define fftw_destroy_plan fftwf_destroy_plan
#define fftw_malloc fftwf_malloc
#define fftw_free fftwf_free
#define fftw_execute fftwf_execute
#define atan2 atan2f
#define sqrt sqrtf
#define cos cosf
#define sin sinf
#else
#define fft_double_type double
#endif /* FFTW_SINGLE_ONLY */

class D_FFTW : public FFTImpl
{
public:
    D_FFTW(int size) :
        m_fplanf(0), m_dplanf(0), m_size(size)
    {
    }

    ~D_FFTW() {
        if (m_fplanf) {
#ifndef NO_THREADING
            m_commonMutex.lock();
#endif
            bool save = false;
            if (m_extantf > 0 && --m_extantf == 0) save = true;
            (void)save;
#ifndef FFTW_DOUBLE_ONLY
            if (save) saveWisdom('f');
#endif
            fftwf_destroy_plan(m_fplanf);
            fftwf_destroy_plan(m_fplani);
            fftwf_free(m_fbuf);
            fftwf_free(m_fpacked);
#ifndef NO_THREADING
            m_commonMutex.unlock();
#endif
        }
        if (m_dplanf) {
#ifndef NO_THREADING
            m_commonMutex.lock();
#endif
            bool save = false;
            if (m_extantd > 0 && --m_extantd == 0) save = true;
            (void)save;
#ifndef FFTW_SINGLE_ONLY
            if (save) saveWisdom('d');
#endif
            fftw_destroy_plan(m_dplanf);
            fftw_destroy_plan(m_dplani);
            fftw_free(m_dbuf);
            fftw_free(m_dpacked);
#ifndef NO_THREADING
            m_commonMutex.unlock();
#endif
        }
#ifndef NO_THREADING
        m_commonMutex.lock();
#endif
        if (m_extantf <= 0 && m_extantd <= 0) {
            fftw_cleanup();
        }
#ifndef NO_THREADING
        m_commonMutex.unlock();
#endif
    }

    FFT::Precisions
    getSupportedPrecisions() const {
#ifdef FFTW_SINGLE_ONLY
        return FFT::SinglePrecision;
#else
#ifdef FFTW_DOUBLE_ONLY
        return FFT::DoublePrecision;
#else
        return FFT::SinglePrecision | FFT::DoublePrecision;
#endif
#endif
    }

    void initFloat() {
        if (m_fplanf) return;
        bool load = false;
#ifndef NO_THREADING
        m_commonMutex.lock();
#endif
        if (m_extantf++ == 0) load = true;
#ifdef FFTW_DOUBLE_ONLY
        if (load) loadWisdom('d');
#else
        if (load) loadWisdom('f');
#endif
        m_fbuf = (fft_float_type *)fftw_malloc(m_size * sizeof(fft_float_type));
        m_fpacked = (fftwf_complex *)fftw_malloc
            ((m_size/2 + 1) * sizeof(fftwf_complex));
        m_fplanf = fftwf_plan_dft_r2c_1d
            (m_size, m_fbuf, m_fpacked, FFTW_MEASURE);
        m_fplani = fftwf_plan_dft_c2r_1d
            (m_size, m_fpacked, m_fbuf, FFTW_MEASURE);
#ifndef NO_THREADING
        m_commonMutex.unlock();
#endif
    }

    void initDouble() {
        if (m_dplanf) return;
        bool load = false;
#ifndef NO_THREADING
        m_commonMutex.lock();
#endif
        if (m_extantd++ == 0) load = true;
#ifdef FFTW_SINGLE_ONLY
        if (load) loadWisdom('f');
#else
        if (load) loadWisdom('d');
#endif
        m_dbuf = (fft_double_type *)fftw_malloc(m_size * sizeof(fft_double_type));
        m_dpacked = (fftw_complex *)fftw_malloc
            ((m_size/2 + 1) * sizeof(fftw_complex));
        m_dplanf = fftw_plan_dft_r2c_1d
            (m_size, m_dbuf, m_dpacked, FFTW_MEASURE);
        m_dplani = fftw_plan_dft_c2r_1d
            (m_size, m_dpacked, m_dbuf, FFTW_MEASURE);
#ifndef NO_THREADING
        m_commonMutex.unlock();
#endif
    }

    void loadWisdom(char type) { wisdom(false, type); }
    void saveWisdom(char type) { wisdom(true, type); }

    void wisdom(bool save, char type) {

#ifdef FFTW_DOUBLE_ONLY
        if (type == 'f') return;
#endif
#ifdef FFTW_SINGLE_ONLY
        if (type == 'd') return;
#endif

        const char *home = getenv("HOME");
        if (!home) return;

        char fn[256];
        snprintf(fn, 256, "%s/%s.%c", home, ".rubberband.wisdom", type);

        FILE *f = fopen(fn, save ? "wb" : "rb");
        if (!f) return;

        if (save) {
            switch (type) {
#ifdef FFTW_DOUBLE_ONLY
            case 'f': break;
#else
            case 'f': fftwf_export_wisdom_to_file(f); break;
#endif
#ifdef FFTW_SINGLE_ONLY
            case 'd': break;
#else
            case 'd': fftw_export_wisdom_to_file(f); break;
#endif
            default: break;
            }
        } else {
            switch (type) {
#ifdef FFTW_DOUBLE_ONLY
            case 'f': break;
#else
            case 'f': fftwf_import_wisdom_from_file(f); break;
#endif
#ifdef FFTW_SINGLE_ONLY
            case 'd': break;
#else
            case 'd': fftw_import_wisdom_from_file(f); break;
#endif
            default: break;
            }
        }

        fclose(f);
    }

    void packFloat(const float *R__ re, const float *R__ im) {
        const int hs = m_size/2;
        fftwf_complex *const R__ fpacked = m_fpacked; 
        for (int i = 0; i <= hs; ++i) {
            fpacked[i][0] = re[i];
        }
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                fpacked[i][1] = im[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                fpacked[i][1] = 0.f;
            }
        }                
    }

    void packDouble(const double *R__ re, const double *R__ im) {
        const int hs = m_size/2;
        fftw_complex *const R__ dpacked = m_dpacked; 
        for (int i = 0; i <= hs; ++i) {
            dpacked[i][0] = re[i];
        }
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                dpacked[i][1] = im[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                dpacked[i][1] = 0.0;
            }
        }
    }

    void unpackFloat(float *R__ re, float *R__ im) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            re[i] = m_fpacked[i][0];
        }
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                im[i] = m_fpacked[i][1];
            }
        }
    }        

    void unpackDouble(double *R__ re, double *R__ im) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            re[i] = m_dpacked[i][0];
        }
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                im[i] = m_dpacked[i][1];
            }
        }
    }        

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        if (!m_dplanf) initDouble();
        const int sz = m_size;
        fft_double_type *const R__ dbuf = m_dbuf;
#ifndef FFTW_SINGLE_ONLY
        if (realIn != dbuf) 
#endif
            for (int i = 0; i < sz; ++i) {
                dbuf[i] = realIn[i];
            }
        fftw_execute(m_dplanf);
        unpackDouble(realOut, imagOut);
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        if (!m_dplanf) initDouble();
        const int sz = m_size;
        fft_double_type *const R__ dbuf = m_dbuf;
#ifndef FFTW_SINGLE_ONLY
        if (realIn != dbuf) 
#endif
            for (int i = 0; i < sz; ++i) {
                dbuf[i] = realIn[i];
            }
        fftw_execute(m_dplanf);
        v_convert(complexOut, (fft_double_type *)m_dpacked, sz + 2);
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        if (!m_dplanf) initDouble();
        fft_double_type *const R__ dbuf = m_dbuf;
        const int sz = m_size;
#ifndef FFTW_SINGLE_ONLY
        if (realIn != dbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                dbuf[i] = realIn[i];
            }
        fftw_execute(m_dplanf);
        v_cartesian_interleaved_to_polar(magOut, phaseOut,
                                         (double *)m_dpacked, m_size/2+1);
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        if (!m_dplanf) initDouble();
        fft_double_type *const R__ dbuf = m_dbuf;
        const int sz = m_size;
#ifndef FFTW_SINGLE_ONLY
        if (realIn != m_dbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                dbuf[i] = realIn[i];
            }
        fftw_execute(m_dplanf);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(m_dpacked[i][0] * m_dpacked[i][0] +
                             m_dpacked[i][1] * m_dpacked[i][1]);
        }
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        if (!m_fplanf) initFloat();
        fft_float_type *const R__ fbuf = m_fbuf;
        const int sz = m_size;
#ifndef FFTW_DOUBLE_ONLY
        if (realIn != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                fbuf[i] = realIn[i];
            }
        fftwf_execute(m_fplanf);
        unpackFloat(realOut, imagOut);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        if (!m_fplanf) initFloat();
        fft_float_type *const R__ fbuf = m_fbuf;
        const int sz = m_size;
#ifndef FFTW_DOUBLE_ONLY
        if (realIn != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                fbuf[i] = realIn[i];
            }
        fftwf_execute(m_fplanf);
        v_convert(complexOut, (fft_float_type *)m_fpacked, sz + 2);
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        if (!m_fplanf) initFloat();
        fft_float_type *const R__ fbuf = m_fbuf;
        const int sz = m_size;
#ifndef FFTW_DOUBLE_ONLY
        if (realIn != fbuf) 
#endif
            for (int i = 0; i < sz; ++i) {
                fbuf[i] = realIn[i];
            }
        fftwf_execute(m_fplanf);
        v_cartesian_interleaved_to_polar(magOut, phaseOut,
                                         (float *)m_fpacked, m_size/2+1);
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        if (!m_fplanf) initFloat();
        fft_float_type *const R__ fbuf = m_fbuf;
        const int sz = m_size;
#ifndef FFTW_DOUBLE_ONLY
        if (realIn != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                fbuf[i] = realIn[i];
            }
        fftwf_execute(m_fplanf);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrtf(m_fpacked[i][0] * m_fpacked[i][0] +
                              m_fpacked[i][1] * m_fpacked[i][1]);
        }
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        if (!m_dplanf) initDouble();
        packDouble(realIn, imagIn);
        fftw_execute(m_dplani);
        const int sz = m_size;
        fft_double_type *const R__ dbuf = m_dbuf;
#ifndef FFTW_SINGLE_ONLY
        if (realOut != dbuf) 
#endif
            for (int i = 0; i < sz; ++i) {
                realOut[i] = dbuf[i];
            }
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        if (!m_dplanf) initDouble();
        v_convert((double *)m_dpacked, complexIn, m_size + 2);
        fftw_execute(m_dplani);
        const int sz = m_size;
        fft_double_type *const R__ dbuf = m_dbuf;
#ifndef FFTW_SINGLE_ONLY
        if (realOut != dbuf) 
#endif
            for (int i = 0; i < sz; ++i) {
                realOut[i] = dbuf[i];
            }
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        if (!m_dplanf) initDouble();
        const int hs = m_size/2;
        fftw_complex *const R__ dpacked = m_dpacked;
        for (int i = 0; i <= hs; ++i) {
            dpacked[i][0] = magIn[i] * cos(phaseIn[i]);
        }
        for (int i = 0; i <= hs; ++i) {
            dpacked[i][1] = magIn[i] * sin(phaseIn[i]);
        }
        fftw_execute(m_dplani);
        const int sz = m_size;
        fft_double_type *const R__ dbuf = m_dbuf;
#ifndef FFTW_SINGLE_ONLY
        if (realOut != dbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                realOut[i] = dbuf[i];
            }
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        if (!m_dplanf) initDouble();
        fft_double_type *const R__ dbuf = m_dbuf;
        fftw_complex *const R__ dpacked = m_dpacked;
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            dpacked[i][0] = log(magIn[i] + 0.000001);
        }
        for (int i = 0; i <= hs; ++i) {
            dpacked[i][1] = 0.0;
        }
        fftw_execute(m_dplani);
        const int sz = m_size;
#ifndef FFTW_SINGLE_ONLY
        if (cepOut != dbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                cepOut[i] = dbuf[i];
            }
    }

    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        if (!m_fplanf) initFloat();
        packFloat(realIn, imagIn);
        fftwf_execute(m_fplani);
        const int sz = m_size;
        fft_float_type *const R__ fbuf = m_fbuf;
#ifndef FFTW_DOUBLE_ONLY
        if (realOut != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                realOut[i] = fbuf[i];
            }
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        if (!m_fplanf) initFloat();
        v_copy((float *)m_fpacked, complexIn, m_size + 2);
        fftwf_execute(m_fplani);
        const int sz = m_size;
        fft_float_type *const R__ fbuf = m_fbuf;
#ifndef FFTW_DOUBLE_ONLY
        if (realOut != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                realOut[i] = fbuf[i];
            }
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        if (!m_fplanf) initFloat();
        const int hs = m_size/2;
        fftwf_complex *const R__ fpacked = m_fpacked;
        for (int i = 0; i <= hs; ++i) {
            fpacked[i][0] = magIn[i] * cosf(phaseIn[i]);
        }
        for (int i = 0; i <= hs; ++i) {
            fpacked[i][1] = magIn[i] * sinf(phaseIn[i]);
        }
        fftwf_execute(m_fplani);
        const int sz = m_size;
        fft_float_type *const R__ fbuf = m_fbuf;
#ifndef FFTW_DOUBLE_ONLY
        if (realOut != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                realOut[i] = fbuf[i];
            }
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        if (!m_fplanf) initFloat();
        const int hs = m_size/2;
        fftwf_complex *const R__ fpacked = m_fpacked;
        for (int i = 0; i <= hs; ++i) {
            fpacked[i][0] = logf(magIn[i] + 0.000001f);
        }
        for (int i = 0; i <= hs; ++i) {
            fpacked[i][1] = 0.f;
        }
        fftwf_execute(m_fplani);
        const int sz = m_size;
        fft_float_type *const R__ fbuf = m_fbuf;
#ifndef FFTW_DOUBLE_ONLY
        if (cepOut != fbuf)
#endif
            for (int i = 0; i < sz; ++i) {
                cepOut[i] = fbuf[i];
            }
    }

private:
    fftwf_plan m_fplanf;
    fftwf_plan m_fplani;
#ifdef FFTW_DOUBLE_ONLY
    double *m_fbuf;
#else
    float *m_fbuf;
#endif
    fftwf_complex *m_fpacked;
    fftw_plan m_dplanf;
    fftw_plan m_dplani;
#ifdef FFTW_SINGLE_ONLY
    float *m_dbuf;
#else
    double *m_dbuf;
#endif
    fftw_complex *m_dpacked;
    const int m_size;
    static int m_extantf;
    static int m_extantd;
#ifndef NO_THREADING
    static Mutex m_commonMutex;
#endif
};

int
D_FFTW::m_extantf = 0;

int
D_FFTW::m_extantd = 0;

#ifndef NO_THREADING
Mutex
D_FFTW::m_commonMutex;
#endif

#endif /* HAVE_FFTW3 */

#ifdef HAVE_SFFT

/*
 Define SFFT_DOUBLE_ONLY to make all uses of SFFT functions be
 double-precision (so "float" FFTs are calculated by casting to
 doubles and using the double-precision SFFT function).

 Define SFFT_SINGLE_ONLY to make all uses of SFFT functions be
 single-precision (so "double" FFTs are calculated by casting to
 floats and using the single-precision SFFT function).

 Neither of these flags is desirable for either performance or
 precision.
*/

//#define SFFT_DOUBLE_ONLY 1
//#define SFFT_SINGLE_ONLY 1

#if defined(SFFT_DOUBLE_ONLY) && defined(SFFT_SINGLE_ONLY)
// Can't meaningfully define both
#error Can only define one of SFFT_DOUBLE_ONLY and SFFT_SINGLE_ONLY
#endif

#ifdef SFFT_DOUBLE_ONLY
#define fft_float_type double
#define FLAG_SFFT_FLOAT SFFT_DOUBLE
#define atan2f atan2
#define sqrtf sqrt
#define cosf cos
#define sinf sin
#define logf log
#else
#define FLAG_SFFT_FLOAT SFFT_FLOAT
#define fft_float_type float
#endif /* SFFT_DOUBLE_ONLY */

#ifdef SFFT_SINGLE_ONLY
#define fft_double_type float
#define FLAG_SFFT_DOUBLE SFFT_FLOAT
#define atan2 atan2f
#define sqrt sqrtf
#define cos cosf
#define sin sinf
#define log logf
#else
#define FLAG_SFFT_DOUBLE SFFT_DOUBLE
#define fft_double_type double
#endif /* SFFT_SINGLE_ONLY */

class D_SFFT : public FFTImpl
{
public:
    D_SFFT(int size) :
        m_fplanf(0), m_fplani(0), m_dplanf(0), m_dplani(0), m_size(size)
    {
    }

    ~D_SFFT() {
        if (m_fplanf) {
            sfft_free(m_fplanf);
            sfft_free(m_fplani);
            deallocate(m_fbuf);
            deallocate(m_fresult);
        }
        if (m_dplanf) {
            sfft_free(m_dplanf);
            sfft_free(m_dplani);
            deallocate(m_dbuf);
            deallocate(m_dresult);
        }
    }

    FFT::Precisions
    getSupportedPrecisions() const {
#ifdef SFFT_SINGLE_ONLY
        return FFT::SinglePrecision;
#else
#ifdef SFFT_DOUBLE_ONLY
        return FFT::DoublePrecision;
#else
        return FFT::SinglePrecision | FFT::DoublePrecision;
#endif
#endif
    }

    void initFloat() {
        if (m_fplanf) return;
        m_fbuf = allocate<fft_float_type>(2 * m_size);
        m_fresult = allocate<fft_float_type>(2 * m_size);
        m_fplanf = sfft_init(m_size, SFFT_FORWARD | FLAG_SFFT_FLOAT);
        m_fplani = sfft_init(m_size, SFFT_BACKWARD | FLAG_SFFT_FLOAT);
        if (!m_fplanf || !m_fplani) {
            if (!m_fplanf) {
                std::cerr << "D_SFFT: Failed to construct forward float transform for size " << m_size << " (check SFFT library's target configuration)" << std::endl;
            } else {
                std::cerr << "D_SFFT: Failed to construct inverse float transform for size " << m_size << " (check SFFT library's target configuration)" << std::endl;
            }
#ifndef NO_EXCEPTIONS
            throw FFT::InternalError;
#else
            abort();
#endif
        }
    }

    void initDouble() {
        if (m_dplanf) return;
        m_dbuf = allocate<fft_double_type>(2 * m_size);
        m_dresult = allocate<fft_double_type>(2 * m_size);
        m_dplanf = sfft_init(m_size, SFFT_FORWARD | FLAG_SFFT_DOUBLE);
        m_dplani = sfft_init(m_size, SFFT_BACKWARD | FLAG_SFFT_DOUBLE);
        if (!m_dplanf || !m_dplani) {
            if (!m_dplanf) {
                std::cerr << "D_SFFT: Failed to construct forward double transform for size " << m_size << " (check SFFT library's target configuration)" << std::endl;
            } else {
                std::cerr << "D_SFFT: Failed to construct inverse double transform for size " << m_size << " (check SFFT library's target configuration)" << std::endl;
            }
#ifndef NO_EXCEPTIONS
            throw FFT::InternalError;
#else
            abort();
#endif
        }
    }

    void packFloat(const float *R__ re, const float *R__ im, fft_float_type *target, int n) {
        for (int i = 0; i < n; ++i) target[i*2] = re[i];
        if (im) {
            for (int i = 0; i < n; ++i) target[i*2+1] = im[i]; 
        } else {
            for (int i = 0; i < n; ++i) target[i*2+1] = 0.f;
        }                
    }

    void packDouble(const double *R__ re, const double *R__ im, fft_double_type *target, int n) {
        for (int i = 0; i < n; ++i) target[i*2] = re[i];
        if (im) {
            for (int i = 0; i < n; ++i) target[i*2+1] = im[i];
        } else {
            for (int i = 0; i < n; ++i) target[i*2+1] = 0.0;
        }                
    }

    void unpackFloat(const fft_float_type *source, float *R__ re, float *R__ im, int n) {
        for (int i = 0; i < n; ++i) re[i] = source[i*2];
        if (im) {
            for (int i = 0; i < n; ++i) im[i] = source[i*2+1];
        }
    }        

    void unpackDouble(const fft_double_type *source, double *R__ re, double *R__ im, int n) {
        for (int i = 0; i < n; ++i) re[i] = source[i*2];
        if (im) {
            for (int i = 0; i < n; ++i) im[i] = source[i*2+1];
        }
    }        

    template<typename T>
    void mirror(T *R__ cplx, int n) {
        for (int i = 1; i <= n/2; ++i) {
            int j = n-i;
            cplx[j*2] = cplx[i*2];
            cplx[j*2+1] = -cplx[i*2+1];
        }
    }

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        if (!m_dplanf) initDouble();
        packDouble(realIn, 0, m_dbuf, m_size);
        sfft_execute(m_dplanf, m_dbuf, m_dresult);
        unpackDouble(m_dresult, realOut, imagOut, m_size/2+1);
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        if (!m_dplanf) initDouble();
        packDouble(realIn, 0, m_dbuf, m_size);
        sfft_execute(m_dplanf, m_dbuf, m_dresult);
        v_convert(complexOut, m_dresult, m_size+2); // i.e. m_size/2+1 complex
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        if (!m_dplanf) initDouble();
        packDouble(realIn, 0, m_dbuf, m_size);
        sfft_execute(m_dplanf, m_dbuf, m_dresult);
        v_cartesian_interleaved_to_polar(magOut, phaseOut,
                                         m_dresult, m_size/2+1);
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        if (!m_dplanf) initDouble();
        packDouble(realIn, 0, m_dbuf, m_size);
        sfft_execute(m_dplanf, m_dbuf, m_dresult);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(m_dresult[i*2] * m_dresult[i*2] +
                             m_dresult[i*2+1] * m_dresult[i*2+1]);
        }
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        if (!m_fplanf) initFloat();
        packFloat(realIn, 0, m_fbuf, m_size);
        sfft_execute(m_fplanf, m_fbuf, m_fresult);
        unpackFloat(m_fresult, realOut, imagOut, m_size/2+1);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        if (!m_fplanf) initFloat();
        packFloat(realIn, 0, m_fbuf, m_size);
        sfft_execute(m_fplanf, m_fbuf, m_fresult);
        v_convert(complexOut, m_fresult, m_size+2); // i.e. m_size/2+1 complex
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        if (!m_fplanf) initFloat();
        packFloat(realIn, 0, m_fbuf, m_size);
        sfft_execute(m_fplanf, m_fbuf, m_fresult);
        v_cartesian_interleaved_to_polar(magOut, phaseOut,
                                         m_fresult, m_size/2+1);
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        if (!m_fplanf) initFloat();
        packFloat(realIn, 0, m_fbuf, m_size);
        sfft_execute(m_fplanf, m_fbuf, m_fresult);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrtf(m_fresult[i*2] * m_fresult[i*2] +
                              m_fresult[i*2+1] * m_fresult[i*2+1]);
        }
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        if (!m_dplanf) initDouble();
        packDouble(realIn, imagIn, m_dbuf, m_size/2+1);
        mirror(m_dbuf, m_size);
        sfft_execute(m_dplani, m_dbuf, m_dresult);
        for (int i = 0; i < m_size; ++i) {
            realOut[i] = m_dresult[i*2];
        }
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        if (!m_dplanf) initDouble();
        v_convert((double *)m_dbuf, complexIn, m_size + 2);
        mirror(m_dbuf, m_size);
        sfft_execute(m_dplani, m_dbuf, m_dresult);
        for (int i = 0; i < m_size; ++i) {
            realOut[i] = m_dresult[i*2];
        }
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        if (!m_dplanf) initDouble();
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_dbuf[i*2] = magIn[i] * cos(phaseIn[i]);
            m_dbuf[i*2+1] = magIn[i] * sin(phaseIn[i]);
        }
        mirror(m_dbuf, m_size);
        sfft_execute(m_dplani, m_dbuf, m_dresult);
        for (int i = 0; i < m_size; ++i) {
            realOut[i] = m_dresult[i*2];
        }
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        if (!m_dplanf) initDouble();
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_dbuf[i*2] = log(magIn[i] + 0.000001);
            m_dbuf[i*2+1] = 0.0;
        }
        mirror(m_dbuf, m_size);
        sfft_execute(m_dplani, m_dbuf, m_dresult);
        for (int i = 0; i < m_size; ++i) {
            cepOut[i] = m_dresult[i*2];
        }
    }

    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        if (!m_fplanf) initFloat();
        packFloat(realIn, imagIn, m_fbuf, m_size/2+1);
        mirror(m_fbuf, m_size);
        sfft_execute(m_fplani, m_fbuf, m_fresult);
        for (int i = 0; i < m_size; ++i) {
            realOut[i] = m_fresult[i*2];
        }
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        if (!m_fplanf) initFloat();
        v_convert((float *)m_fbuf, complexIn, m_size + 2);
        mirror(m_fbuf, m_size);
        sfft_execute(m_fplani, m_fbuf, m_fresult);
        for (int i = 0; i < m_size; ++i) {
            realOut[i] = m_fresult[i*2];
        }
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        if (!m_fplanf) initFloat();
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_fbuf[i*2] = magIn[i] * cosf(phaseIn[i]);
            m_fbuf[i*2+1] = magIn[i] * sinf(phaseIn[i]);
        }
        mirror(m_fbuf, m_size);
        sfft_execute(m_fplani, m_fbuf, m_fresult);
        for (int i = 0; i < m_size; ++i) {
            realOut[i] = m_fresult[i*2];
        }
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        if (!m_fplanf) initFloat();
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_fbuf[i*2] = logf(magIn[i] + 0.00001);
            m_fbuf[i*2+1] = 0.0f;
        }
        sfft_execute(m_fplani, m_fbuf, m_fresult);
        for (int i = 0; i < m_size; ++i) {
            cepOut[i] = m_fresult[i*2];
        }
    }

private:
    sfft_plan_t *m_fplanf;
    sfft_plan_t *m_fplani;
    fft_float_type *m_fbuf;
    fft_float_type *m_fresult;

    sfft_plan_t *m_dplanf;
    sfft_plan_t *m_dplani;
    fft_double_type *m_dbuf;
    fft_double_type *m_dresult;

    const int m_size;
};

#endif /* HAVE_SFFT */

#ifdef USE_KISSFFT

class D_KISSFFT : public FFTImpl
{
public:
    D_KISSFFT(int size) :
        m_size(size),
        m_fplanf(0),  
        m_fplani(0)
    {
#ifdef FIXED_POINT
#error KISSFFT is not configured for float values
#endif
        if (sizeof(kiss_fft_scalar) != sizeof(float)) {
            std::cerr << "ERROR: KISSFFT is not configured for float values"
                      << std::endl;
        }

        m_fbuf = new kiss_fft_scalar[m_size + 2];
        m_fpacked = new kiss_fft_cpx[m_size + 2];
        m_fplanf = kiss_fftr_alloc(m_size, 0, NULL, NULL);
        m_fplani = kiss_fftr_alloc(m_size, 1, NULL, NULL);
    }

    ~D_KISSFFT() {
        kiss_fftr_free(m_fplanf);
        kiss_fftr_free(m_fplani);
        kiss_fft_cleanup();

        delete[] m_fbuf;
        delete[] m_fpacked;
    }

    FFT::Precisions
    getSupportedPrecisions() const {
        return FFT::SinglePrecision;
    }

    void initFloat() { }
    void initDouble() { }

    void packFloat(const float *R__ re, const float *R__ im) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_fpacked[i].r = re[i];
        }
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                m_fpacked[i].i = im[i];
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                m_fpacked[i].i = 0.f;
            }
        }
    }

    void unpackFloat(float *R__ re, float *R__ im) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            re[i] = m_fpacked[i].r;
        }
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                im[i] = m_fpacked[i].i;
            }
        }
    }        

    void packDouble(const double *R__ re, const double *R__ im) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            m_fpacked[i].r = float(re[i]);
        }
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                m_fpacked[i].i = float(im[i]);
            }
        } else {
            for (int i = 0; i <= hs; ++i) {
                m_fpacked[i].i = 0.f;
            }
        }
    }

    void unpackDouble(double *R__ re, double *R__ im) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            re[i] = double(m_fpacked[i].r);
        }
        if (im) {
            for (int i = 0; i <= hs; ++i) {
                im[i] = double(m_fpacked[i].i);
            }
        }
    }        

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {

        v_convert(m_fbuf, realIn, m_size);
        kiss_fftr(m_fplanf, m_fbuf, m_fpacked);
        unpackDouble(realOut, imagOut);
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {

        v_convert(m_fbuf, realIn, m_size);
        kiss_fftr(m_fplanf, m_fbuf, m_fpacked);
        v_convert(complexOut, (float *)m_fpacked, m_size + 2);
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {

        for (int i = 0; i < m_size; ++i) {
            m_fbuf[i] = float(realIn[i]);
        }

        kiss_fftr(m_fplanf, m_fbuf, m_fpacked);

        const int hs = m_size/2;

        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(double(m_fpacked[i].r) * double(m_fpacked[i].r) +
                             double(m_fpacked[i].i) * double(m_fpacked[i].i));
        }

        for (int i = 0; i <= hs; ++i) {
            phaseOut[i] = atan2(double(m_fpacked[i].i), double(m_fpacked[i].r));
        }
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {

        for (int i = 0; i < m_size; ++i) {
            m_fbuf[i] = float(realIn[i]);
        }

        kiss_fftr(m_fplanf, m_fbuf, m_fpacked);

        const int hs = m_size/2;

        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(double(m_fpacked[i].r) * double(m_fpacked[i].r) +
                             double(m_fpacked[i].i) * double(m_fpacked[i].i));
        }
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {

        kiss_fftr(m_fplanf, realIn, m_fpacked);
        unpackFloat(realOut, imagOut);
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {

        kiss_fftr(m_fplanf, realIn, (kiss_fft_cpx *)complexOut);
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {

        kiss_fftr(m_fplanf, realIn, m_fpacked);

        const int hs = m_size/2;

        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrtf(m_fpacked[i].r * m_fpacked[i].r +
                              m_fpacked[i].i * m_fpacked[i].i);
        }

        for (int i = 0; i <= hs; ++i) {
            phaseOut[i] = atan2f(m_fpacked[i].i, m_fpacked[i].r);
        }
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {

        kiss_fftr(m_fplanf, realIn, m_fpacked);

        const int hs = m_size/2;

        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrtf(m_fpacked[i].r * m_fpacked[i].r +
                              m_fpacked[i].i * m_fpacked[i].i);
        }
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {

        packDouble(realIn, imagIn);

        kiss_fftri(m_fplani, m_fpacked, m_fbuf);

        for (int i = 0; i < m_size; ++i) {
            realOut[i] = m_fbuf[i];
        }
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {

        v_convert((float *)m_fpacked, complexIn, m_size + 2);

        kiss_fftri(m_fplani, m_fpacked, m_fbuf);

        for (int i = 0; i < m_size; ++i) {
            realOut[i] = m_fbuf[i];
        }
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {

        const int hs = m_size/2;

        for (int i = 0; i <= hs; ++i) {
            m_fpacked[i].r = float(magIn[i] * cos(phaseIn[i]));
            m_fpacked[i].i = float(magIn[i] * sin(phaseIn[i]));
        }

        kiss_fftri(m_fplani, m_fpacked, m_fbuf);

        for (int i = 0; i < m_size; ++i) {
            realOut[i] = m_fbuf[i];
        }
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {

        const int hs = m_size/2;

        for (int i = 0; i <= hs; ++i) {
            m_fpacked[i].r = float(log(magIn[i] + 0.000001));
            m_fpacked[i].i = 0.0f;
        }

        kiss_fftri(m_fplani, m_fpacked, m_fbuf);

        for (int i = 0; i < m_size; ++i) {
            cepOut[i] = m_fbuf[i];
        }
    }
    
    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {

        packFloat(realIn, imagIn);
        kiss_fftri(m_fplani, m_fpacked, realOut);
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {

        v_copy((float *)m_fpacked, complexIn, m_size + 2);
        kiss_fftri(m_fplani, m_fpacked, realOut);
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {

        const int hs = m_size/2;

        for (int i = 0; i <= hs; ++i) {
            m_fpacked[i].r = magIn[i] * cosf(phaseIn[i]);
            m_fpacked[i].i = magIn[i] * sinf(phaseIn[i]);
        }

        kiss_fftri(m_fplani, m_fpacked, realOut);
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {

        const int hs = m_size/2;

        for (int i = 0; i <= hs; ++i) {
            m_fpacked[i].r = logf(magIn[i] + 0.000001f);
            m_fpacked[i].i = 0.0f;
        }

        kiss_fftri(m_fplani, m_fpacked, cepOut);
    }

private:
    const int m_size;
    kiss_fftr_cfg m_fplanf;
    kiss_fftr_cfg m_fplani;
    kiss_fft_scalar *m_fbuf;
    kiss_fft_cpx *m_fpacked;
};

#endif /* USE_KISSFFT */

#ifdef USE_BUILTIN_FFT

class D_Cross : public FFTImpl
{
public:
    D_Cross(int size) : m_size(size), m_table(0) {
        
        m_a = new double[size];
        m_b = new double[size];
        m_c = new double[size];
        m_d = new double[size];

        m_table = new int[m_size];
    
        int bits;
        int i, j, k, m;

        for (i = 0; ; ++i) {
            if (m_size & (1 << i)) {
                bits = i;
                break;
            }
        }
        
        for (i = 0; i < m_size; ++i) {
            
            m = i;
            
            for (j = k = 0; j < bits; ++j) {
                k = (k << 1) | (m & 1);
                m >>= 1;
            }
            
            m_table[i] = k;
        }
    }

    ~D_Cross() {
        delete[] m_table;
        delete[] m_a;
        delete[] m_b;
        delete[] m_c;
        delete[] m_d;
    }

    FFT::Precisions
    getSupportedPrecisions() const {
        return FFT::DoublePrecision;
    }

    void initFloat() { }
    void initDouble() { }

    void forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut) {
        basefft(false, realIn, 0, m_c, m_d);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) realOut[i] = m_c[i];
        if (imagOut) {
            for (int i = 0; i <= hs; ++i) imagOut[i] = m_d[i];
        }
    }

    void forwardInterleaved(const double *R__ realIn, double *R__ complexOut) {
        basefft(false, realIn, 0, m_c, m_d);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) complexOut[i*2] = m_c[i];
        for (int i = 0; i <= hs; ++i) complexOut[i*2+1] = m_d[i];
    }

    void forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut) {
        basefft(false, realIn, 0, m_c, m_d);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(m_c[i] * m_c[i] + m_d[i] * m_d[i]);
            phaseOut[i] = atan2(m_d[i], m_c[i]) ;
        }
    }

    void forwardMagnitude(const double *R__ realIn, double *R__ magOut) {
        basefft(false, realIn, 0, m_c, m_d);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(m_c[i] * m_c[i] + m_d[i] * m_d[i]);
        }
    }

    void forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut) {
        for (int i = 0; i < m_size; ++i) m_a[i] = realIn[i];
        basefft(false, m_a, 0, m_c, m_d);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) realOut[i] = m_c[i];
        if (imagOut) {
            for (int i = 0; i <= hs; ++i) imagOut[i] = m_d[i];
        }
    }

    void forwardInterleaved(const float *R__ realIn, float *R__ complexOut) {
        for (int i = 0; i < m_size; ++i) m_a[i] = realIn[i];
        basefft(false, m_a, 0, m_c, m_d);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) complexOut[i*2] = m_c[i];
        for (int i = 0; i <= hs; ++i) complexOut[i*2+1] = m_d[i];
    }

    void forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut) {
        for (int i = 0; i < m_size; ++i) m_a[i] = realIn[i];
        basefft(false, m_a, 0, m_c, m_d);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(m_c[i] * m_c[i] + m_d[i] * m_d[i]);
            phaseOut[i] = atan2(m_d[i], m_c[i]) ;
        }
    }

    void forwardMagnitude(const float *R__ realIn, float *R__ magOut) {
        for (int i = 0; i < m_size; ++i) m_a[i] = realIn[i];
        basefft(false, m_a, 0, m_c, m_d);
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            magOut[i] = sqrt(m_c[i] * m_c[i] + m_d[i] * m_d[i]);
        }
    }

    void inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            double real = realIn[i];
            double imag = imagIn[i];
            m_a[i] = real;
            m_b[i] = imag;
            if (i > 0) {
                m_a[m_size-i] = real;
                m_b[m_size-i] = -imag;
            }
        }
        basefft(true, m_a, m_b, realOut, m_d);
    }

    void inverseInterleaved(const double *R__ complexIn, double *R__ realOut) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            double real = complexIn[i*2];
            double imag = complexIn[i*2+1];
            m_a[i] = real;
            m_b[i] = imag;
            if (i > 0) {
                m_a[m_size-i] = real;
                m_b[m_size-i] = -imag;
            }
        }
        basefft(true, m_a, m_b, realOut, m_d);
    }

    void inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            double real = magIn[i] * cos(phaseIn[i]);
            double imag = magIn[i] * sin(phaseIn[i]);
            m_a[i] = real;
            m_b[i] = imag;
            if (i > 0) {
                m_a[m_size-i] = real;
                m_b[m_size-i] = -imag;
            }
        }
        basefft(true, m_a, m_b, realOut, m_d);
    }

    void inverseCepstral(const double *R__ magIn, double *R__ cepOut) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            double real = log(magIn[i] + 0.000001);
            m_a[i] = real;
            m_b[i] = 0.0;
            if (i > 0) {
                m_a[m_size-i] = real;
                m_b[m_size-i] = 0.0;
            }
        }
        basefft(true, m_a, m_b, cepOut, m_d);
    }

    void inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            float real = realIn[i];
            float imag = imagIn[i];
            m_a[i] = real;
            m_b[i] = imag;
            if (i > 0) {
                m_a[m_size-i] = real;
                m_b[m_size-i] = -imag;
            }
        }
        basefft(true, m_a, m_b, m_c, m_d);
        for (int i = 0; i < m_size; ++i) realOut[i] = m_c[i];
    }

    void inverseInterleaved(const float *R__ complexIn, float *R__ realOut) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            float real = complexIn[i*2];
            float imag = complexIn[i*2+1];
            m_a[i] = real;
            m_b[i] = imag;
            if (i > 0) {
                m_a[m_size-i] = real;
                m_b[m_size-i] = -imag;
            }
        }
        basefft(true, m_a, m_b, m_c, m_d);
        for (int i = 0; i < m_size; ++i) realOut[i] = m_c[i];
    }

    void inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            float real = magIn[i] * cosf(phaseIn[i]);
            float imag = magIn[i] * sinf(phaseIn[i]);
            m_a[i] = real;
            m_b[i] = imag;
            if (i > 0) {
                m_a[m_size-i] = real;
                m_b[m_size-i] = -imag;
            }
        }
        basefft(true, m_a, m_b, m_c, m_d);
        for (int i = 0; i < m_size; ++i) realOut[i] = m_c[i];
    }

    void inverseCepstral(const float *R__ magIn, float *R__ cepOut) {
        const int hs = m_size/2;
        for (int i = 0; i <= hs; ++i) {
            float real = logf(magIn[i] + 0.000001);
            m_a[i] = real;
            m_b[i] = 0.0;
            if (i > 0) {
                m_a[m_size-i] = real;
                m_b[m_size-i] = 0.0;
            }
        }
        basefft(true, m_a, m_b, m_c, m_d);
        for (int i = 0; i < m_size; ++i) cepOut[i] = m_c[i];
    }

private:
    const int m_size;
    int *m_table;
    double *m_a;
    double *m_b;
    double *m_c;
    double *m_d;
    void basefft(bool inverse, const double *R__ ri, const double *R__ ii, double *R__ ro, double *R__ io);
};

void
D_Cross::basefft(bool inverse, const double *R__ ri, const double *R__ ii, double *R__ ro, double *R__ io)
{
    if (!ri || !ro || !io) return;

    int i, j, k, m;
    int blockSize, blockEnd;

    double tr, ti;

    double angle = 2.0 * M_PI;
    if (inverse) angle = -angle;

    const int n = m_size;

    if (ii) {
	for (i = 0; i < n; ++i) {
	    ro[m_table[i]] = ri[i];
        }
	for (i = 0; i < n; ++i) {
	    io[m_table[i]] = ii[i];
	}
    } else {
	for (i = 0; i < n; ++i) {
	    ro[m_table[i]] = ri[i];
        }
	for (i = 0; i < n; ++i) {
	    io[m_table[i]] = 0.0;
	}
    }

    blockEnd = 1;

    for (blockSize = 2; blockSize <= n; blockSize <<= 1) {

	double delta = angle / (double)blockSize;
	double sm2 = -sin(-2 * delta);
	double sm1 = -sin(-delta);
	double cm2 = cos(-2 * delta);
	double cm1 = cos(-delta);
	double w = 2 * cm1;
	double ar[3], ai[3];

	for (i = 0; i < n; i += blockSize) {

	    ar[2] = cm2;
	    ar[1] = cm1;

	    ai[2] = sm2;
	    ai[1] = sm1;

	    for (j = i, m = 0; m < blockEnd; j++, m++) {

		ar[0] = w * ar[1] - ar[2];
		ar[2] = ar[1];
		ar[1] = ar[0];

		ai[0] = w * ai[1] - ai[2];
		ai[2] = ai[1];
		ai[1] = ai[0];

		k = j + blockEnd;
		tr = ar[0] * ro[k] - ai[0] * io[k];
		ti = ar[0] * io[k] + ai[0] * ro[k];

		ro[k] = ro[j] - tr;
		io[k] = io[j] - ti;

		ro[j] += tr;
		io[j] += ti;
	    }
	}

	blockEnd = blockSize;
    }

/* fftw doesn't rescale, so nor will we

    if (inverse) {

	double denom = (double)n;

	for (i = 0; i < n; i++) {
	    ro[i] /= denom;
	    io[i] /= denom;
	}
    }
*/
}

#endif /* USE_BUILTIN_FFT */

} /* end namespace FFTs */

std::string
FFT::m_implementation;

std::set<std::string>
FFT::getImplementations()
{
    std::set<std::string> impls;
#ifdef HAVE_IPP
    impls.insert("ipp");
#endif
#ifdef HAVE_FFTW3
    impls.insert("fftw");
#endif
#ifdef USE_KISSFFT
    impls.insert("kissfft");
#endif
#ifdef HAVE_VDSP
    impls.insert("vdsp");
#endif
#ifdef HAVE_MEDIALIB
    impls.insert("medialib");
#endif
#ifdef HAVE_OPENMAX
    impls.insert("openmax");
#endif
#ifdef HAVE_SFFT
    impls.insert("sfft");
#endif
#ifdef USE_BUILTIN_FFT
    impls.insert("cross");
#endif
    return impls;
}

void
FFT::pickDefaultImplementation()
{
    if (m_implementation != "") return;

    std::set<std::string> impls = getImplementations();

    std::string best = "cross";
    if (impls.find("kissfft") != impls.end()) best = "kissfft";
    if (impls.find("medialib") != impls.end()) best = "medialib";
    if (impls.find("openmax") != impls.end()) best = "openmax";
    if (impls.find("sfft") != impls.end()) best = "sfft";
    if (impls.find("fftw") != impls.end()) best = "fftw";
    if (impls.find("vdsp") != impls.end()) best = "vdsp";
    if (impls.find("ipp") != impls.end()) best = "ipp";
    
    m_implementation = best;
}

std::string
FFT::getDefaultImplementation()
{
    return m_implementation;
}

void
FFT::setDefaultImplementation(std::string i)
{
    m_implementation = i;
}

FFT::FFT(int size, int debugLevel) :
    d(0)
{
    if ((size < 2) ||
        (size & (size-1))) {
        std::cerr << "FFT::FFT(" << size << "): power-of-two sizes only supported, minimum size 2" << std::endl;
#ifndef NO_EXCEPTIONS
        throw InvalidSize;
#else
        abort();
#endif
    }

    if (m_implementation == "") pickDefaultImplementation();
    std::string impl = m_implementation;

    if (debugLevel > 0) {
        std::cerr << "FFT::FFT(" << size << "): using implementation: "
                  << impl << std::endl;
    }

    if (impl == "ipp") {
#ifdef HAVE_IPP
        d = new FFTs::D_IPP(size);
#endif
    } else if (impl == "fftw") {
#ifdef HAVE_FFTW3
        d = new FFTs::D_FFTW(size);
#endif
    } else if (impl == "kissfft") {        
#ifdef USE_KISSFFT
        d = new FFTs::D_KISSFFT(size);
#endif
    } else if (impl == "vdsp") {
#ifdef HAVE_VDSP
        d = new FFTs::D_VDSP(size);
#endif
    } else if (impl == "medialib") {
#ifdef HAVE_MEDIALIB
        d = new FFTs::D_MEDIALIB(size);
#endif
    } else if (impl == "openmax") {
#ifdef HAVE_OPENMAX
        d = new FFTs::D_OPENMAX(size);
#endif
    } else if (impl == "sfft") {
#ifdef HAVE_SFFT
        d = new FFTs::D_SFFT(size);
#endif
    } else if (impl == "cross") {
#ifdef USE_BUILTIN_FFT
        d = new FFTs::D_Cross(size);
#endif
    }

    if (!d) {
        std::cerr << "FFT::FFT(" << size << "): ERROR: implementation "
                  << impl << " is not compiled in" << std::endl;
#ifndef NO_EXCEPTIONS
        throw InvalidImplementation;
#else
        abort();
#endif
    }
}

FFT::~FFT()
{
    delete d;
}

#ifndef NO_EXCEPTIONS
#define CHECK_NOT_NULL(x) \
    if (!(x)) { \
        std::cerr << "FFT: ERROR: Null argument " #x << std::endl;  \
        throw NullArgument; \
    }
#else
#define CHECK_NOT_NULL(x) \
    if (!(x)) { \
        std::cerr << "FFT: ERROR: Null argument " #x << std::endl;  \
        std::cerr << "FFT: Would be throwing NullArgument here, if exceptions were not disabled" << std::endl;  \
        return; \
    }
#endif

void
FFT::forward(const double *R__ realIn, double *R__ realOut, double *R__ imagOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(realOut);
    CHECK_NOT_NULL(imagOut);
    d->forward(realIn, realOut, imagOut);
}

void
FFT::forwardInterleaved(const double *R__ realIn, double *R__ complexOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(complexOut);
    d->forwardInterleaved(realIn, complexOut);
}

void
FFT::forwardPolar(const double *R__ realIn, double *R__ magOut, double *R__ phaseOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(magOut);
    CHECK_NOT_NULL(phaseOut);
    d->forwardPolar(realIn, magOut, phaseOut);
}

void
FFT::forwardMagnitude(const double *R__ realIn, double *R__ magOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(magOut);
    d->forwardMagnitude(realIn, magOut);
}

void
FFT::forward(const float *R__ realIn, float *R__ realOut, float *R__ imagOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(realOut);
    CHECK_NOT_NULL(imagOut);
    d->forward(realIn, realOut, imagOut);
}

void
FFT::forwardInterleaved(const float *R__ realIn, float *R__ complexOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(complexOut);
    d->forwardInterleaved(realIn, complexOut);
}

void
FFT::forwardPolar(const float *R__ realIn, float *R__ magOut, float *R__ phaseOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(magOut);
    CHECK_NOT_NULL(phaseOut);
    d->forwardPolar(realIn, magOut, phaseOut);
}

void
FFT::forwardMagnitude(const float *R__ realIn, float *R__ magOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(magOut);
    d->forwardMagnitude(realIn, magOut);
}

void
FFT::inverse(const double *R__ realIn, const double *R__ imagIn, double *R__ realOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(imagIn);
    CHECK_NOT_NULL(realOut);
    d->inverse(realIn, imagIn, realOut);
}

void
FFT::inverseInterleaved(const double *R__ complexIn, double *R__ realOut)
{
    CHECK_NOT_NULL(complexIn);
    CHECK_NOT_NULL(realOut);
    d->inverseInterleaved(complexIn, realOut);
}

void
FFT::inversePolar(const double *R__ magIn, const double *R__ phaseIn, double *R__ realOut)
{
    CHECK_NOT_NULL(magIn);
    CHECK_NOT_NULL(phaseIn);
    CHECK_NOT_NULL(realOut);
    d->inversePolar(magIn, phaseIn, realOut);
}

void
FFT::inverseCepstral(const double *R__ magIn, double *R__ cepOut)
{
    CHECK_NOT_NULL(magIn);
    CHECK_NOT_NULL(cepOut);
    d->inverseCepstral(magIn, cepOut);
}

void
FFT::inverse(const float *R__ realIn, const float *R__ imagIn, float *R__ realOut)
{
    CHECK_NOT_NULL(realIn);
    CHECK_NOT_NULL(imagIn);
    CHECK_NOT_NULL(realOut);
    d->inverse(realIn, imagIn, realOut);
}

void
FFT::inverseInterleaved(const float *R__ complexIn, float *R__ realOut)
{
    CHECK_NOT_NULL(complexIn);
    CHECK_NOT_NULL(realOut);
    d->inverseInterleaved(complexIn, realOut);
}

void
FFT::inversePolar(const float *R__ magIn, const float *R__ phaseIn, float *R__ realOut)
{
    CHECK_NOT_NULL(magIn);
    CHECK_NOT_NULL(phaseIn);
    CHECK_NOT_NULL(realOut);
    d->inversePolar(magIn, phaseIn, realOut);
}

void
FFT::inverseCepstral(const float *R__ magIn, float *R__ cepOut)
{
    CHECK_NOT_NULL(magIn);
    CHECK_NOT_NULL(cepOut);
    d->inverseCepstral(magIn, cepOut);
}

void
FFT::initFloat() 
{
    d->initFloat();
}

void
FFT::initDouble() 
{
    d->initDouble();
}

FFT::Precisions
FFT::getSupportedPrecisions() const
{
    return d->getSupportedPrecisions();
}

#ifdef FFT_MEASUREMENT

std::string
FFT::tune()
{
    std::ostringstream os;
    os << "FFT::tune()..." << std::endl;

    std::vector<int> sizes;
    std::map<FFTImpl *, int> candidates;
    std::map<int, int> wins;

    sizes.push_back(512);
    sizes.push_back(1024);
    sizes.push_back(4096);
    
    for (unsigned int si = 0; si < sizes.size(); ++si) {

        int size = sizes[si];

        while (!candidates.empty()) {
            delete candidates.begin()->first;
            candidates.erase(candidates.begin());
        }

        FFTImpl *d;
        
#ifdef HAVE_IPP
        std::cerr << "Constructing new IPP FFT object for size " << size << "..." << std::endl;
        d = new FFTs::D_IPP(size);
        d->initFloat();
        d->initDouble();
        candidates[d] = 0;
#endif
        
#ifdef HAVE_FFTW3
        os << "Constructing new FFTW3 FFT object for size " << size << "..." << std::endl;
        d = new FFTs::D_FFTW(size);
        d->initFloat();
        d->initDouble();
        candidates[d] = 1;
#endif

#ifdef USE_KISSFFT
        os << "Constructing new KISSFFT object for size " << size << "..." << std::endl;
        d = new FFTs::D_KISSFFT(size);
        d->initFloat();
        d->initDouble();
        candidates[d] = 2;
#endif        

#ifdef USE_BUILTIN_FFT
        os << "Constructing new Cross FFT object for size " << size << "..." << std::endl;
        d = new FFTs::D_Cross(size);
        d->initFloat();
        d->initDouble();
        candidates[d] = 3;
#endif
        
#ifdef HAVE_VDSP
        os << "Constructing new vDSP FFT object for size " << size << "..." << std::endl;
        d = new FFTs::D_VDSP(size);
        d->initFloat();
        d->initDouble();
        candidates[d] = 4;
#endif
        
#ifdef HAVE_MEDIALIB
        std::cerr << "Constructing new MediaLib FFT object for size " << size << "..." << std::endl;
        d = new FFTs::D_MEDIALIB(size);
        d->initFloat();
        d->initDouble();
        candidates[d] = 5;
#endif
        
#ifdef HAVE_OPENMAX
        os << "Constructing new OpenMAX FFT object for size " << size << "..." << std::endl;
        d = new FFTs::D_OPENMAX(size);
        d->initFloat();
        d->initDouble();
        candidates[d] = 6;
#endif
        
#ifdef HAVE_SFFT
        os << "Constructing new SFFT FFT object for size " << size << "..." << std::endl;
        d = new FFTs::D_SFFT(size);
//        d->initFloat();
        d->initDouble();
        candidates[d] = 6;
#endif

        os << "CLOCKS_PER_SEC = " << CLOCKS_PER_SEC << std::endl;
        float divisor = float(CLOCKS_PER_SEC) / 1000.f;
        
        os << "Timing order is: ";
        for (std::map<FFTImpl *, int>::iterator ci = candidates.begin();
             ci != candidates.end(); ++ci) {
            os << ci->second << " ";
        }
        os << std::endl;

        int iterations = 500;
        os << "Iterations: " << iterations << std::endl;

        double *da = new double[size];
        double *db = new double[size];
        double *dc = new double[size];
        double *dd = new double[size];
        double *di = new double[size + 2];
        double *dj = new double[size + 2];

        float *fa = new float[size];
        float *fb = new float[size];
        float *fc = new float[size];
        float *fd = new float[size];
        float *fi = new float[size + 2];
        float *fj = new float[size + 2];

        for (int type = 0; type < 16; ++type) {
    
            //!!!
            if ((type > 3 && type < 8) ||
                (type > 11)) {
                continue;
            }

            if (type > 7) {
                // inverse transform: bigger inputs, to simulate the
                // fact that the forward transform is unscaled
                for (int i = 0; i < size; ++i) {
                    da[i] = drand48() * size;
                    fa[i] = da[i];
                    db[i] = drand48() * size;
                    fb[i] = db[i];
                }
            } else {    
                for (int i = 0; i < size; ++i) {
                    da[i] = drand48();
                    fa[i] = da[i];
                    db[i] = drand48();
                    fb[i] = db[i];
                }
            }
                
            for (int i = 0; i < size + 2; ++i) {
                di[i] = drand48();
                fi[i] = di[i];
            }

            int low = -1;
            int lowscore = 0;

            const char *names[] = {

                "Forward Cartesian Double",
                "Forward Interleaved Double",
                "Forward Polar Double",
                "Forward Magnitude Double",
                "Forward Cartesian Float",
                "Forward Interleaved Float",
                "Forward Polar Float",
                "Forward Magnitude Float",

                "Inverse Cartesian Double",
                "Inverse Interleaved Double",
                "Inverse Polar Double",
                "Inverse Cepstral Double",
                "Inverse Cartesian Float",
                "Inverse Interleaved Float",
                "Inverse Polar Float",
                "Inverse Cepstral Float"
            };
            os << names[type] << " :: ";

            for (std::map<FFTImpl *, int>::iterator ci = candidates.begin();
                 ci != candidates.end(); ++ci) {

                FFTImpl *d = ci->first;

                double mean = 0;

                clock_t start = clock();
                
                for (int i = 0; i < iterations; ++i) {

                    if (i == 0) {
                        for (int j = 0; j < size; ++j) {
                            dc[j] = 0;
                            dd[j] = 0;
                            fc[j] = 0;
                            fd[j] = 0;
                            fj[j] = 0;
                            dj[j] = 0;
                        }
                    }

                    switch (type) {
                    case 0: d->forward(da, dc, dd); break;
                    case 1: d->forwardInterleaved(da, dj); break;
                    case 2: d->forwardPolar(da, dc, dd); break;
                    case 3: d->forwardMagnitude(da, dc); break;
                    case 4: d->forward(fa, fc, fd); break;
                    case 5: d->forwardInterleaved(fa, fj); break;
                    case 6: d->forwardPolar(fa, fc, fd); break;
                    case 7: d->forwardMagnitude(fa, fc); break;
                    case 8: d->inverse(da, db, dc); break;
                    case 9: d->inverseInterleaved(di, dc); break;
                    case 10: d->inversePolar(da, db, dc); break;
                    case 11: d->inverseCepstral(da, dc); break;
                    case 12: d->inverse(fa, fb, fc); break;
                    case 13: d->inverseInterleaved(fi, fc); break;
                    case 14: d->inversePolar(fa, fb, fc); break;
                    case 15: d->inverseCepstral(fa, fc); break;
                    }

                    if (i == 0) {
                        mean = 0;
                        for (int j = 0; j < size; ++j) {
                            mean += dc[j];
                            mean += dd[j];
                            mean += fc[j];
                            mean += fd[j];
                            mean += fj[j];
                            mean += dj[j];
                        }
                        mean /= size * 6;
                    }
                }

                clock_t end = clock();

                os << float(end - start)/divisor << " (" << mean << ") ";

                if (low == -1 || (end - start) < lowscore) {
                    low = ci->second;
                    lowscore = end - start;
                }
            }

            os << std::endl;

            os << "  size " << size << ", type " << type << ": fastest is " << low << " (time " << float(lowscore)/divisor << ")" << std::endl;

            wins[low]++;
        }
        
        delete[] fa;
        delete[] fb;
        delete[] fc;
        delete[] fd;
        delete[] da;
        delete[] db;
        delete[] dc;
        delete[] dd;
    }

    while (!candidates.empty()) {
        delete candidates.begin()->first;
        candidates.erase(candidates.begin());
    }

    int bestscore = 0;
    int best = -1;

    for (std::map<int, int>::iterator wi = wins.begin(); wi != wins.end(); ++wi) {
        if (best == -1 || wi->second > bestscore) {
            best = wi->first;
            bestscore = wi->second;
        }
    }

    os << "overall winner is " << best << " with " << bestscore << " wins" << std::endl;

    return os.str();
}

#endif

}
