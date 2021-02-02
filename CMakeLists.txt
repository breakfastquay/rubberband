cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(rubberband VERSION 2.1.1)

option(BUILD_EXECUTABLE "Build the Rubberband executable" OFF)
option(BUILD_JNI "Build the Java Native Interface (JNI)" OFF)
option(OPTIMIZE "Build Rubberband optimized for the build machine's CPU" OFF)
option(USE_FFTW "Use FFTW instead of kissFFT" OFF)
option(USE_LIBSAMPLERATE "Use libsamplerate instead of Speex" ON)

add_library(rubberband
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

if(BUILD_JNI)
  add_library(rubberband-jni src/jni/RubberBandStretcherJNI.cpp)
  target_include_directories(rubberband-jni PRIVATE src)

  find_package(JNI REQUIRED)
  target_include_directories(rubberband-jni PRIVATE ${JNI_INCLUDE_DIRS})
  target_link_libraries(rubberband-jni PRIVATE ${JNI_LIBRARIES})
endif()


#
# Dependencies
#

# Threads
find_package(Threads REQUIRED)
target_link_libraries(rubberband PRIVATE Threads::Threads)

# FFTW / KissFFT
if(USE_FFTW)
  message(STATUS "Using FFT library: fftw")
  find_path(FFTW3_INCLUDEDIR fftw3.h REQUIRED)
  find_library(FFTW3_LIBRARY fftw3 REQUIRED)
  target_include_directories(rubberband PRIVATE ${FFTW3_INCLUDEDIR})
  target_link_libraries(rubberband PUBLIC ${FFTW3_LIBRARY})
  target_compile_definitions(rubberband PRIVATE HAVE_FFTW3 FFTW_DOUBLE_ONLY)
else()
  message(STATUS "Using resampling library: kissfft")
  target_compile_definitions(rubberband PRIVATE USE_KISSFFT)
  target_sources(rubberband PRIVATE src/kissfft/kiss_fft.c src/kissfft/kiss_fftr.c)
endif()

# libsamplerate / Speex
if(USE_LIBSAMPLERATE)
  message(STATUS "Using resampling library: libsamplerate")
  find_path(SAMPLERATE_INCLUDEDIR samplerate.h REQUIRED)
  find_library(SAMPLERATE_LIBRARY samplerate REQUIRED)
  target_include_directories(rubberband PRIVATE ${SAMPLERATE_INCLUDEDIR})
  target_link_libraries(rubberband PUBLIC ${SAMPLERATE_LIBRARY})
  target_compile_definitions(rubberband PRIVATE HAVE_LIBSAMPLERATE)
else()
  message(STATUS "Using resampling library: speex")
  target_compile_definitions(rubberband PRIVATE USE_SPEEX)
  target_sources(rubberband PRIVATE src/speex/resample.c)
endif()

#
# Handle general flags, include paths, etc.
#
include(CheckCXXCompilerFlag)

if(OPTIMIZE)
  check_cxx_compiler_flag("-march=native" MARCH_NATIVE_SUPPORTED)
  if(MARCH_NATIVE_SUPPORTED AND NOT MSVC)
    target_compile_options(rubberband PRIVATE -march=native)
  endif()
endif()

include(GNUInstallDirs)
target_include_directories(rubberband
  PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/rubberband>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/rubberband>
  PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

target_compile_definitions(rubberband
  PRIVATE
    $<$<BOOL:${CMAKE_USE_PTHREADS_INIT}>:USE_PTHREADS>
    MALLOC_IS_ALIGNED
    NO_THREAD_CHECKS
    NO_TIMING
)

if(APPLE)
  if(CMAKE_OSX_DEPLOYMENT_TARGET VERSION_LESS 10.11)
    message(FATAL_ERROR "macOS deployment target must be >= 10.11")
  endif()
  target_compile_definitions(rubberband PRIVATE HAVE_VDSP)

  find_library(ACCELERATE_FRAMEWORK NAMES Accelerate REQUIRED)
  target_link_libraries(rubberband PUBLIC ${ACCELERATE_FRAMEWORK})
endif()

if(WIN32)
    target_compile_definitions(rubberband PRIVATE _WINDOWS)
    if(CMAKE_SIZEOF_VOID_P EQUAL 4)
        target_compile_definitions(rubberband PRIVATE WIN32)
    endif()

    if(MSVC)
      target_compile_definitions(rubberband PRIVATE
          __MSVC__
          _LIB
          NOMINMAX
          _USE_MATH_DEFINES
          $<$<CONFIG:Release>:NDEBUG>
          $<$<CONFIG:Debug>:_DEBUG NO_TIMING NO_THREAD_CHECKS>
      )
      set_target_properties(rubberband PROPERTIES WINDOWS_EXPORT_ALL_SYMBOLS ON)
    endif()
endif()

#
# Command-Line Program
#
if(BUILD_EXECUTABLE)
  find_library(SNDFILE_LIBRARY sndfile REQUIRED)
  add_executable(rubberband-program main/main.cpp)
  target_include_directories(rubberband-program PRIVATE src Threads::Threads)
  target_link_libraries(rubberband-program PRIVATE rubberband ${SNDFILE_LIBRARY})
endif()

#
# Installation
#
install(
  TARGETS rubberband
  EXPORT rubberband-targets
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
)

# Header files
install(
    DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/rubberband"
    TYPE INCLUDE
)

# pkg-config file
if(UNIX)
    set(PREFIX "${CMAKE_INSTALL_PREFIX}")
    configure_file("${CMAKE_CURRENT_SOURCE_DIR}/rubberband.pc.in" "${CMAKE_CURRENT_BINARY_DIR}/rubberband.pc" @ONLY)
    install(
        FILES "${CMAKE_CURRENT_BINARY_DIR}/rubberband.pc"
        DESTINATION "${CMAKE_INSTALL_LIBDIR}/pkgconfig"
    )
endif()

# CMake config
include(CMakePackageConfigHelpers)
set(RUBBERBAND_INSTALL_CMAKEDIR "lib/cmake/rubberband")
install(
  EXPORT rubberband-targets
  FILE rubberband-targets.cmake
  NAMESPACE rubberband::
  DESTINATION "${RUBBERBAND_INSTALL_CMAKEDIR}"
)
configure_package_config_file(rubberband-config.cmake.in
  "${CMAKE_CURRENT_BINARY_DIR}/rubberband-config.cmake"
  INSTALL_DESTINATION "${RUBBERBAND_INSTALL_CMAKEDIR}"
)
write_basic_package_version_file("${CMAKE_CURRENT_BINARY_DIR}/rubberband-config-version.cmake"
  VERSION "${CMAKE_PROJECT_VERSION}"
  COMPATIBILITY SameMajorVersion
)
install(
  FILES
      "${CMAKE_CURRENT_BINARY_DIR}/rubberband-config.cmake"
      "${CMAKE_CURRENT_BINARY_DIR}/rubberband-config-version.cmake"
  DESTINATION "${RUBBERBAND_INSTALL_CMAKEDIR}"
)
