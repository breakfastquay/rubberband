cmake_minimum_required(VERSION 3.9 FATAL_ERROR)
project(rubberband VERSION 2.1.1)

option(RUBBERBAND_EXAMPLE "Build the Rubberband example app" OFF)
option(RUBBERBAND_OPTIMIZED "Build Rubberband optimized for the build machine's CPU" OFF)
### Dependencies ###

# ARM cross-compilers do not support march=native, only mtune=native
set(MARCH_NATIVE_FLAG)

include(CheckCXXCompilerFlag)

if(RUBBERBAND_OPTIMIZED)
  check_cxx_compiler_flag("-march=native" MARCH_NATIVE_SUPPORTED)
  if(MARCH_NATIVE_SUPPORTED)
    set(MARCH_NATIVE_FLAG "-march=native")
  else()
  endif()
endif()

find_package(Threads REQUIRED)

# FFTW3
find_path(FFTW3_INCLUDEDIR fftw3.h)
find_library(FFTW3_LIBRARY fftw3)

# libsamplerate
if(NOT TARGET samplerate)
  find_path(SAMPLERATE_INCLUDEDIR samplerate.h)
  find_library(SAMPLERATE_LIBRARY samplerate)
endif()

# Accelerate framework
if(APPLE)
  find_library(ACCELERATE_FRAMEWORK NAMES Accelerate)
endif()

### Create our target ###
add_library(rubberband)

### Handle general flags, include paths, etc ###
target_compile_options(rubberband
  PRIVATE
    $<$<BOOL:${APPLE}>:-mmacosx-version-min=10.11>
    $<$<AND:$<NOT:$<CXX_COMPILER_ID:MSVC>>,$<CONFIG:Release>>:-ffast-math ${MARCH_NATIVE_FLAG} -O3 -ftree-vectorize>
)

target_compile_definitions(rubberband
  PRIVATE
    $<$<BOOL:${UNIX}>:USE_PTHREADS>
    $<$<BOOL:${APPLE}>:HAVE_VDSP>
    $<$<BOOL:$<CXX_COMPILER_ID:MSVC>>:__MSVC__ WIN32 _LIB NOMINMAX _USE_MATH_DEFINES USE_KISSFFT>
    $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Release>>:NDEBUG>
    $<$<AND:$<CXX_COMPILER_ID:MSVC>,$<CONFIG:Debug>>:_DEBUG NO_TIMING NO_THREAD_CHECKS>
    MALLOC_IS_ALIGNED
    NO_THREAD_CHECKS
    NO_TIMING
)

target_include_directories(rubberband
  PUBLIC
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/rubberband>
      $<INSTALL_INTERFACE:rubberband>
  PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src>
)

if(WIN32)
    target_compile_definitions(rubberband PRIVATE _WINDOWS)
    if(CMAKE_SIZEOF_VOID_P EQUAL 4)
        target_compile_definitions(rubberband PRIVATE WIN32)
    endif()

    if(MSVC AND BUILD_SHARED_LIBS)
        target_sources(rubberband PRIVATE
            rubberband-dll/dllmain.cpp
            rubberband-dll/rubberband-dll.cpp
            rubberband-dll/stdafx.cpp
        )
        target_include_directories(rubberband PRIVATE rubberband-dll)
        target_compile_definitions(rubberband PRIVATE RUBBERBANDDLL_EXPORTS _USRDLL)
    endif()
    target_include_directories(rubberband PRIVATE $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/rubberband-dll>)
endif()

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

# Apple framework
if(APPLE)
  target_link_libraries(rubberband
    PRIVATE
      "${ACCELERATE_FRAMEWORK}"
  )
endif()

# Handle resampling library
if(TARGET samplerate)
  target_compile_definitions(rubberband
    PRIVATE
      HAVE_LIBSAMPLERATE
  )

  target_link_libraries(rubberband
    PRIVATE
      samplerate
  )
elseif(SAMPLERATE_INCLUDEDIR AND SAMPLERATE_LIBRARY)
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

include(GNUInstallDirs)
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
