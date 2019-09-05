cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project(rubberband VERSION 2.1.1)

option(RUBBERBAND_EXAMPLE "Build the Rubberband example app" OFF)

### Dependencies ###
find_package(Threads REQUIRED)

# FFTW3
find_path(FFTW3_INCLUDEDIR fftw3.h)
find_library(FFTW3_LIBRARY fftw3)

# libsamplerate
find_path(SAMPLERATE_INCLUDEDIR samplerate.h)
find_library(SAMPLERATE_LIBRARY samplerate)

### Create our target ###
add_library(rubberband)

### Handle general flags, include paths, etc ###
target_compile_options(rubberband
  PRIVATE
    $<$<BOOL:${APPLE}>:-mmacosx-version-min=10.7>
    $<$<AND:$<NOT:$<CXX_COMPILER_ID:MSVC>>,$<CONFIG:Release>>:-ffast-math -march=native -O3 -ftree-vectorize>
)

target_compile_definitions(rubberband
  PRIVATE
    $<$<BOOL:${UNIX}>:USE_PTHREADS>
    $<$<BOOL:${APPLE}>:HAVE_VDSP>
    MALLOC_IS_ALIGNED
    NO_THREAD_CHECKS
    NO_TIMING
)

target_include_directories(rubberband
  PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/rubberband>
      $<INSTALL_INTERFACE:rubberband>
  PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

target_link_libraries(rubberband
  PRIVATE
    Threads::Threads
)

# Handle FFT library
if(FFTW3_INCLUDEDIR AND FFTW3_LIBRARY)
  target_include_directories(rubberband
    PRIVATE
      ${FFTW3_INCLUDEDIR}
  )
  target_compile_definitions(rubberband
    PRIVATE
      HAVE_FFTW3
      FFTW_DOUBLE_ONLY
  )
  target_link_libraries(rubberband
    PUBLIC
      ${FFTW3_LIBRARY}
  )
else()
  target_compile_definitions(rubberband
    PRIVATE
      USE_KISSFFT
  )
  target_sources(rubberband
    PRIVATE
      src/kissfft/kiss_fft.c
      src/kissfft/kiss_fftr.c
  )
endif()

# Handle resampling library
if(SAMPLERATE_INCLUDEDIR AND SAMPLERATE_LIBRARY)
  target_include_directories(rubberband
    PRIVATE
      ${SAMPLERATE_INCLUDEDIR}
  )
  target_compile_definitions(rubberband
    PRIVATE
      HAVE_LIBSAMPLERATE
  )
  target_link_libraries(rubberband
    PUBLIC
      ${SAMPLERATE_LIBRARY}
  )
else()
  target_compile_definitions(rubberband
    PRIVATE
      USE_SPEEX
  )
  target_sources(rubberband
    PRIVATE
      src/speex/resample.c
  )
endif()

### Source files ###
target_sources(rubberband
  PRIVATE
    src/RubberBandStretcher.cpp
    src/StretchCalculator.cpp
    src/StretcherChannelData.cpp
    src/StretcherImpl.cpp
    src/StretcherProcess.cpp
    src/audiocurves/CompoundAudioCurve.cpp
    src/audiocurves/ConstantAudioCurve.cpp
    src/audiocurves/HighFrequencyAudioCurve.cpp
    src/audiocurves/PercussiveAudioCurve.cpp
    src/audiocurves/SilentAudioCurve.cpp
    src/audiocurves/SpectralDifferenceAudioCurve.cpp
    src/base/Profiler.cpp
    src/dsp/AudioCurveCalculator.cpp
    src/dsp/FFT.cpp
    src/dsp/Resampler.cpp
    src/rubberband-c.cpp
    src/system/Allocators.cpp
    src/system/Thread.cpp
    src/system/VectorOpsComplex.cpp
    src/system/sysutils.cpp
)

### Example ###
if(RUBBERBAND_EXAMPLE)
  add_executable(rubberband_test main/main.cpp)
  target_include_directories(rubberband_test PRIVATE src Threads::Threads)
  target_link_libraries(rubberband_test PRIVATE rubberband -lsndfile)
endif()

### Installation ###

install(
  TARGETS rubberband
  EXPORT rubberband-exports
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib)

install(EXPORT rubberband-exports
        DESTINATION lib/cmake/rubberband)
export(EXPORT rubberband-exports)
