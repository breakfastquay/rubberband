cmake_minimum_required(VERSION 3.4.1)

set (RUBBERBAND_PATH .)
set (RUBBERBAND_SRC_PATH ${RUBBERBAND_PATH}/src)

set (RUBBERBAND_JNI_FILES
        ${RUBBERBAND_SRC_PATH}/jni/RubberBandStretcherJNI.cpp
)

set (RUBBERBAND_SRC_FILES
        ${RUBBERBAND_SRC_PATH}/base/Profiler.cpp
        ${RUBBERBAND_SRC_PATH}/system/Thread.cpp
        ${RUBBERBAND_SRC_PATH}/system/Allocators.cpp
        ${RUBBERBAND_SRC_PATH}/system/sysutils.cpp
        ${RUBBERBAND_SRC_PATH}/system/VectorOpsComplex.cpp
        ${RUBBERBAND_SRC_PATH}/StretcherChannelData.cpp
        ${RUBBERBAND_SRC_PATH}/dsp/AudioCurveCalculator.cpp
        ${RUBBERBAND_SRC_PATH}/dsp/FFT.cpp
        ${RUBBERBAND_SRC_PATH}/dsp/Resampler.cpp
        ${RUBBERBAND_SRC_PATH}/audiocurves/SilentAudioCurve.cpp
        ${RUBBERBAND_SRC_PATH}/audiocurves/CompoundAudioCurve.cpp
        ${RUBBERBAND_SRC_PATH}/audiocurves/HighFrequencyAudioCurve.cpp
        ${RUBBERBAND_SRC_PATH}/audiocurves/SpectralDifferenceAudioCurve.cpp
        ${RUBBERBAND_SRC_PATH}/audiocurves/ConstantAudioCurve.cpp
        ${RUBBERBAND_SRC_PATH}/audiocurves/PercussiveAudioCurve.cpp
        ${RUBBERBAND_SRC_PATH}/StretcherImpl.cpp
        ${RUBBERBAND_SRC_PATH}/StretcherProcess.cpp
        ${RUBBERBAND_SRC_PATH}/StretchCalculator.cpp
        ${RUBBERBAND_SRC_PATH}/RubberBandStretcher.cpp
        ${RUBBERBAND_SRC_PATH}/rubberband-c.cpp
        ${RUBBERBAND_SRC_PATH}/speex/resample.c
)

set (LOCAL_SRC_FILES
        ${LOCAL_SRC_FILES}
            ${RUBBERBAND_JNI_FILES}
                ${RUBBERBAND_SRC_FILES}
)

set (LOCAL_SRC_FILES
        ${LOCAL_SRC_FILES}
            ${RUBBERBAND_SRC_PATH}/kissfft/kiss_fft.c
                ${RUBBERBAND_SRC_PATH}/kissfft/kiss_fftr.c
)

# Creates and names a library, sets it as either STATIC
# or SHARED, and provides the relative paths to its source code.
# You can define multiple libraries, and CMake builds them for you.
# Gradle automatically packages shared libraries with your APK.

add_library( # Sets the name of the library.
        rubberband-jni

        # Sets the library as a shared library.
        SHARED

        # Provides a relative path to your source file(s).
        ${LOCAL_SRC_FILES})

target_compile_options (rubberband-jni
        PRIVATE

        "$<$<CONFIG:DEBUG>:-g>"
        "$<$<CONFIG:DEBUG>:-mfloat-abi=softfp>"
        "$<$<CONFIG:DEBUG>:-DWANT_TIMING>"
        "$<$<CONFIG:DEBUG>:-DFFT_MEASUREMENT>"

        "$<$<CONFIG:RELEASE>:-O3>"
        "$<$<CONFIG:RELEASE>:-mfpu=neon>"
        "$<$<CONFIG:RELEASE>:-mfloat-abi=softfp>"
        "$<$<CONFIG:RELEASE>:-ffast-math>"
        "$<$<CONFIG:RELEASE>:-ftree-vectorize>"
        "$<$<CONFIG:RELEASE>:-freciprocal-math>"
        "$<$<CONFIG:RELEASE>:-fsingle-precision-constant>"
        "$<$<CONFIG:RELEASE>:-D__ARM_ARCH_7__>"
        "$<$<CONFIG:RELEASE>:-DUSE_POMMIER_MATHFUN>"
        "$<$<CONFIG:RELEASE>:-DNO_TIMING>"
        "$<$<CONFIG:RELEASE>:-DNO_TIMING_COMPLETE_NOOP>"
        -Wall
        -DUSE_SPEEX
        -DUSE_KISSFFT
        -DPROCESS_SAMPLE_TYPE=float
        -DLACK_POSIX_MEMALIGN
        -DUSE_OWN_ALIGNED_MALLOC
        -DLACK_SINCOS
        -DNO_EXCEPTIONS
        -DNO_THREADING
        -DNO_THREAD_CHECKS
)

set (TARGET_ARCH_ABI armeabi-v7a)
set (LOCAL_ARM_MODE arm)
set (LOCAL_ARM_NEON true)

include_directories(${RUBBERBAND_PATH} ${RUBBERBAND_SRC_PATH})

target_link_libraries( # Specifies the target library.
        rubberband-jni

        # Links the target library to the log library
        # included in the NDK.
        log)
