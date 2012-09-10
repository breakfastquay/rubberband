/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */
/* Copyright Chris Cannam - All Rights Reserved */

#include "rubberband/RubberBandStretcher.h"

#include "system/Allocators.h"

#include <jni.h>

using namespace RubberBand;

extern "C" {

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    dispose
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_dispose
  (JNIEnv *, jobject);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    reset
 * Signature: ()V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_reset
  (JNIEnv *, jobject);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    setTimeRatio
 * Signature: (D)V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_setTimeRatio
  (JNIEnv *, jobject, jdouble);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    setPitchScale
 * Signature: (D)V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_setPitchScale
  (JNIEnv *, jobject, jdouble);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    getChannelCount
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_getChannelCount
  (JNIEnv *, jobject);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    getTimeRatio
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_getTimeRatio
  (JNIEnv *, jobject);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    getPitchScale
 * Signature: ()D
 */
JNIEXPORT jdouble JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_getPitchScale
  (JNIEnv *, jobject);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    getLatency
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_getLatency
  (JNIEnv *, jobject);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    setTransientsOption
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_setTransientsOption
  (JNIEnv *, jobject, jint);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    setDetectorOption
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_setDetectorOption
  (JNIEnv *, jobject, jint);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    setPhaseOption
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_setPhaseOption
  (JNIEnv *, jobject, jint);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    setFormantOption
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_setFormantOption
  (JNIEnv *, jobject, jint);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    setPitchOption
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_setPitchOption
  (JNIEnv *, jobject, jint);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    setExpectedInputDuration
 * Signature: (J)V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_setExpectedInputDuration
  (JNIEnv *, jobject, jlong);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    setMaxProcessSize
 * Signature: (I)V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_setMaxProcessSize
  (JNIEnv *, jobject, jint);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    getSamplesRequired
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_getSamplesRequired
  (JNIEnv *, jobject);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    study
 * Signature: ([[FZ)V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_study
  (JNIEnv *, jobject, jobjectArray, jint, jint, jboolean);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    process
 * Signature: ([[FZ)V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_process
  (JNIEnv *, jobject, jobjectArray, jint, jint, jboolean);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    available
 * Signature: ()I
 */
JNIEXPORT jint JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_available
  (JNIEnv *, jobject);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    retrieve
 * Signature: (I)[[F
 */
JNIEXPORT jint JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_retrieve
  (JNIEnv *, jobject, jobjectArray, jint, jint);

/*
 * Class:     com_breakfastquay_rubberband_RubberBandStretcher
 * Method:    initialise
 * Signature: (IIIDD)V
 */
JNIEXPORT void JNICALL Java_com_breakfastquay_rubberband_RubberBandStretcher_initialise
  (JNIEnv *, jobject, jint, jint, jint, jdouble, jdouble);

}

RubberBandStretcher *
getStretcher(JNIEnv *env, jobject obj)
{
    jclass c = env->GetObjectClass(obj);
    jfieldID fid = env->GetFieldID(c, "handle", "J");
    jlong handle = env->GetLongField(obj, fid);
    return (RubberBandStretcher *)handle;
}

void
setStretcher(JNIEnv *env, jobject obj, RubberBandStretcher *stretcher)
{
    jclass c = env->GetObjectClass(obj);
    jfieldID fid = env->GetFieldID(c, "handle", "J");
    jlong handle = (jlong)stretcher;
    env->SetLongField(obj, fid, handle);
}

void
Java_com_breakfastquay_rubberband_RubberBandStretcher_initialise(JNIEnv *env, jobject obj, jint sampleRate, jint channels, jint options, jdouble initialTimeRatio, jdouble initialPitchScale)
{
    setStretcher(env, obj, new RubberBandStretcher
                 (sampleRate, channels, options, initialTimeRatio, initialPitchScale));
}

void
Java_com_breakfastquay_rubberband_RubberBandStretcher_dispose(JNIEnv *env, jobject obj)
{
    delete getStretcher(env, obj);
    setStretcher(env, obj, 0);
}

void
Java_com_breakfastquay_rubberband_RubberBandStretcher_reset(JNIEnv *env, jobject obj)
{
    getStretcher(env, obj)->reset();
}

void
Java_com_breakfastquay_rubberband_RubberBandStretcher_setTimeRatio(JNIEnv *env, jobject obj, jdouble ratio)
{
    getStretcher(env, obj)->setTimeRatio(ratio);
}

void
Java_com_breakfastquay_rubberband_RubberBandStretcher_setPitchScale(JNIEnv *env, jobject obj, jdouble scale)
{
    getStretcher(env, obj)->setPitchScale(scale);
}

jint
Java_com_breakfastquay_rubberband_RubberBandStretcher_getChannelCount(JNIEnv *env, jobject obj)
{
    return getStretcher(env, obj)->getChannelCount();
}

jdouble
Java_com_breakfastquay_rubberband_RubberBandStretcher_getTimeRatio(JNIEnv *env, jobject obj)
{
    return getStretcher(env, obj)->getTimeRatio();
}

jdouble
Java_com_breakfastquay_rubberband_RubberBandStretcher_getPitchScale(JNIEnv *env, jobject obj)
{
    return getStretcher(env, obj)->getPitchScale();
}

jint
Java_com_breakfastquay_rubberband_RubberBandStretcher_getLatency(JNIEnv *env, jobject obj)
{
    return getStretcher(env, obj)->getLatency();
}

void
Java_com_breakfastquay_rubberband_RubberBandStretcher_setTransientsOption(JNIEnv *env, jobject obj, jint options)
{
    getStretcher(env, obj)->setTransientsOption(options);
}

void
Java_com_breakfastquay_rubberband_RubberBandStretcher_setDetectorOption(JNIEnv *env, jobject obj, jint options)
{
    getStretcher(env, obj)->setDetectorOption(options);
}

void
Java_com_breakfastquay_rubberband_RubberBandStretcher_setPhaseOption(JNIEnv *env, jobject obj, jint options)
{
    getStretcher(env, obj)->setPhaseOption(options);
}

void
Java_com_breakfastquay_rubberband_RubberBandStretcher_setFormantOption(JNIEnv *env, jobject obj, jint options)
{
    getStretcher(env, obj)->setFormantOption(options);
}

void
Java_com_breakfastquay_rubberband_RubberBandStretcher_setPitchOption(JNIEnv *env, jobject obj, jint options)
{
    getStretcher(env, obj)->setPitchOption(options);
}

void
Java_com_breakfastquay_rubberband_RubberBandStretcher_setExpectedInputDuration(JNIEnv *env, jobject obj, jlong duration)
{
    getStretcher(env, obj)->setExpectedInputDuration(duration);
}

jint
Java_com_breakfastquay_rubberband_RubberBandStretcher_getSamplesRequired(JNIEnv *env, jobject obj)
{
    return getStretcher(env, obj)->getSamplesRequired();
}

void
Java_com_breakfastquay_rubberband_RubberBandStretcher_study(JNIEnv *env, jobject obj, jobjectArray data, jint offset, jint n, jboolean final)
{
    int channels = env->GetArrayLength(data);
    float **arr = allocate<float *>(channels);
    float **input = allocate<float *>(channels);
    int samples = 0;
    for (int c = 0; c < channels; ++c) {
        jfloatArray cdata = (jfloatArray)env->GetObjectArrayElement(data, c);
        arr[c] = env->GetFloatArrayElements(cdata, 0);
        input[c] = arr[c] + offset;
    }

    getStretcher(env, obj)->study(input, n, final);

    for (int c = 0; c < channels; ++c) {
        jfloatArray cdata = (jfloatArray)env->GetObjectArrayElement(data, c);
        env->ReleaseFloatArrayElements(cdata, arr[c], 0);
    }
}

void
Java_com_breakfastquay_rubberband_RubberBandStretcher_process(JNIEnv *env, jobject obj, jobjectArray data, jint offset, jint n, jboolean final)
{
    int channels = env->GetArrayLength(data);
    float **arr = allocate<float *>(channels);
    float **input = allocate<float *>(channels);
    int samples = 0;
    for (int c = 0; c < channels; ++c) {
        jfloatArray cdata = (jfloatArray)env->GetObjectArrayElement(data, c);
        arr[c] = env->GetFloatArrayElements(cdata, 0);
        input[c] = arr[c] + offset;
    }

    getStretcher(env, obj)->process(input, n, final);

    for (int c = 0; c < channels; ++c) {
        jfloatArray cdata = (jfloatArray)env->GetObjectArrayElement(data, c);
        env->ReleaseFloatArrayElements(cdata, arr[c], 0);
    }

    deallocate(input);
    deallocate(arr);
}

jint
Java_com_breakfastquay_rubberband_RubberBandStretcher_available(JNIEnv *env, jobject obj)
{
    return getStretcher(env, obj)->available();
}

jint
Java_com_breakfastquay_rubberband_RubberBandStretcher_retrieve(JNIEnv *env, jobject obj, jobjectArray output, jint offset, jint n)
{
    RubberBandStretcher *stretcher = getStretcher(env, obj);
    size_t channels = stretcher->getChannelCount();
    
    float **outbuf = allocate_channels<float>(channels, n);
    size_t retrieved = stretcher->retrieve(outbuf, n);

    for (int c = 0; c < channels; ++c) {
        jfloatArray cdata = (jfloatArray)env->GetObjectArrayElement(output, c);
        env->SetFloatArrayRegion(cdata, offset, retrieved, outbuf[c]);
    }
    
    deallocate_channels(outbuf, channels);
    return retrieved;
}

