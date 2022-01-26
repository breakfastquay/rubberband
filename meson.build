
project(
  'Rubber Band Library',
  'c', 'cpp',
  version: '2.0.2',
  license: 'GPL-2.0-or-later',
  default_options: [
    'cpp_std=c++14',
    'buildtype=release',
    'default_library=both',
    'b_ndebug=if-release',
    'b_lundef=true',
  ],
  meson_version: '>= 0.53.0'
)

rubberband_dynamic_library_version = '2.1.7'

system = host_machine.system()
architecture = host_machine.cpu_family()
cpp = meson.get_compiler('cpp')

pkg = import('pkgconfig')


# Define the project source sets

public_headers = [
  'rubberband/rubberband-c.h',
  'rubberband/RubberBandStretcher.h',
]

library_sources = [
  'src/rubberband-c.cpp',
  'src/RubberBandStretcher.cpp',
  'src/StretcherProcess.cpp',
  'src/StretchCalculator.cpp',
  'src/base/Profiler.cpp',
  'src/dsp/AudioCurveCalculator.cpp',
  'src/audiocurves/CompoundAudioCurve.cpp',
  'src/audiocurves/SpectralDifferenceAudioCurve.cpp',
  'src/audiocurves/HighFrequencyAudioCurve.cpp',
  'src/audiocurves/SilentAudioCurve.cpp',
  'src/audiocurves/ConstantAudioCurve.cpp',
  'src/audiocurves/PercussiveAudioCurve.cpp',
  'src/dsp/Resampler.cpp',
  'src/dsp/FFT.cpp',
  'src/system/Allocators.cpp',
  'src/system/sysutils.cpp',
  'src/system/Thread.cpp',
  'src/StretcherChannelData.cpp',
  'src/StretcherImpl.cpp',
]

jni_sources = [
  'src/jni/RubberBandStretcherJNI.cpp',
]

java_sources = [
  'com/breakfastquay/rubberband/RubberBandStretcher.java',
]

program_sources = [
  'main/main.cpp',
]

if system == 'windows'
  program_sources += [
    'src/getopt/getopt.c',
    'src/getopt/getopt_long.c'
  ]
endif

vamp_sources = [
  'vamp/RubberBandVampPlugin.cpp',
  'vamp/libmain.cpp',
]

ladspa_sources = [
  'ladspa-lv2/libmain-ladspa.cpp',
]

lv2_sources = [
  'ladspa-lv2/libmain-lv2.cpp',
]

general_include_dirs = [
  'rubberband',
  'src',
]


# Scan for any dependencies we may use later; all are optional

extra_include_args = []
foreach d: get_option('extra_include_dirs')
  extra_include_args += [ '-I' + d ]
endforeach

fftw3_dep = dependency('fftw3', version: '>= 3.0.0', required: false)
samplerate_dep = dependency('samplerate', version: '>= 0.1.8', required: false)
sndfile_dep = dependency('sndfile', version: '>= 1.0.16', required: false)
vamp_dep = dependency('vamp-sdk', version: '>= 2.9', required: false)
thread_dep = dependency('threads')
have_ladspa = cpp.has_header('ladspa.h', args: extra_include_args)
have_lv2 = cpp.has_header('lv2.h', args: extra_include_args)
have_sincos = cpp.has_function('sincos',
                               prefix: '#define _GNU_SOURCE\n#include <math.h>',
                               args: '-lm')
have_jni = cpp.has_header('jni.h', args: extra_include_args)

javac = find_program('javac', required: false)
jar = find_program('jar', required: false)


# Check FFT and resampler options and set up dependencies and paths

feature_dependencies = []
feature_defines = []
feature_libraries = []
feature_sources = []
pkgconfig_requirements = []
pkgconfig_libraries = []
arch_flags = []

config_summary = {}
target_summary = {}

resampler = get_option('resampler')
fft = get_option('fft')
ipp_path = get_option('ipp_path')
ipp_needed = false

if fft == 'auto'
  if system == 'darwin'
    fft = 'vdsp'
  else
    fft = 'builtin'
  endif
endif

if resampler == 'auto'
  resampler = 'builtin'
endif

if fft == 'builtin'
  config_summary += { 'FFT': 'Built-in' }
  message('For FFT: using built-in implementation')
  if fftw3_dep.found()
    message('(to use FFTW instead, reconfigure with -Dfft=fftw)')
  endif
  feature_defines += ['-DUSE_BUILTIN_FFT']

elif fft == 'kissfft'
  config_summary += { 'FFT': 'KissFFT' }
  message('For FFT: using KissFFT')
  if fftw3_dep.found()
    message('(to use FFTW instead, reconfigure with -Dfft=fftw)')
  endif
  feature_sources += ['src/kissfft/kiss_fft.c', 'src/kissfft/kiss_fftr.c']
  feature_defines += ['-DHAVE_KISSFFT']
  general_include_dirs += 'src/kissfft'

elif fft == 'fftw'
  if fftw3_dep.found()
    config_summary += { 'FFT': 'FFTW' }
    message('For FFT: using FFTW')
    pkgconfig_requirements += fftw3_dep
  else 
    fftw_dep = cpp.find_library('fftw3',
                                dirs: get_option('extra_lib_dirs'),
                                has_headers: ['fftw3.h'],
                                header_args: extra_include_args,
                                required: true)
  endif
  feature_dependencies += fftw3_dep
  feature_defines += ['-DHAVE_FFTW3', '-DFFTW_DOUBLE_ONLY']

elif fft == 'vdsp'
  config_summary += { 'FFT': 'vDSP' }
  message('For FFT: using vDSP')
  feature_defines += ['-DHAVE_VDSP']
  feature_libraries += ['-framework', 'Accelerate']
  pkgconfig_libraries += ['-framework', 'Accelerate']
  
elif fft == 'ipp'
  if ipp_path != ''
    config_summary += { 'FFT': 'Intel IPP' }
    message('For FFT: using IPP')
    message('IPP path defined as ' + ipp_path)
  else 
    error('For FFT: IPP selected, but ipp_path not specified')
  endif
  ipp_needed = true

else
  error('Unknown or unsupported FFT option: ' + fft)

endif # fft

if resampler == 'builtin'
  config_summary += { 'Resampler': 'Built-in' }
  message('For resampler: using built-in implementation')
  if samplerate_dep.found()
    message('(to use libsamplerate instead, reconfigure with -Dresampler=libsamplerate)')
  endif
  library_sources += 'src/dsp/BQResampler.cpp'
  feature_defines += ['-DUSE_BQRESAMPLER']

elif resampler == 'libsamplerate'
  if samplerate_dep.found()
    config_summary += { 'Resampler': 'libsamplerate' }
    message('For resampler: using libsamplerate')
    pkgconfig_requirements += samplerate_dep
  else 
    samplerate_dep = cpp.find_library('samplerate',
                                      dirs: get_option('extra_lib_dirs'),
                                      has_headers: ['samplerate.h'],
                                      header_args: extra_include_args,
                                      required: true)
  endif
  feature_dependencies += samplerate_dep
  feature_defines += ['-DHAVE_LIBSAMPLERATE']
  
elif resampler == 'speex'
  config_summary += { 'Resampler': 'Speex' }
  message('For resampler: using Speex')
  message('(consider libsamplerate if time-varying pitch shift is required)')
  feature_sources += ['src/speex/resample.c']
  feature_defines += ['-DUSE_SPEEX']
  
elif resampler == 'ipp'
  if ipp_path != ''
    config_summary += { 'Resampler': 'Intel IPP' }
    message('For resampler: using IPP')
    message('(consider libsamplerate if time-varying pitch shift is required)')
    message('IPP path defined as ' + ipp_path)
  else 
    error('For resampler: IPP selected, but ipp_path not specified')
  endif
  ipp_needed = true

else
  error('Unknown or unsupported resampler option: ' + resampler)

endif # resampler

if not have_sincos
  feature_defines += [ '-DLACK_SINCOS' ]
endif

if ipp_needed
  feature_defines += [
    '-DHAVE_IPP',
    '-DUSE_IPP_STATIC',
    '-I' + ipp_path / 'include'
  ]
  if architecture == 'x86'
    feature_libraries += [
      '-L' + ipp_path / 'lib/ia32',
    ]
  elif architecture == 'x86_64'
    feature_libraries += [
      '-L' + ipp_path / 'lib/intel64',
    ]
  else
    error('IPP is not supported for this architecture')
  endif
  if system == 'windows'
    feature_libraries += [
      '-lippsmt', '-lippvmmt', '-lippcoremt',
    ]
  elif system == 'linux'
    feature_libraries += [
      '-Wl,-Bstatic', '-lipps', '-lippvm', '-lippcore', '-Wl,-Bdynamic',
    ]
  else   
    feature_libraries += [
      '-lipps', '-lippvm', '-lippcore',
    ]
  endif
endif # ipp_needed

if not vamp_dep.found()
  vamp_dep = cpp.find_library('VampPluginSDK',
                              dirs: get_option('extra_lib_dirs'),
                              has_headers: ['vamp-sdk.h'],
                              header_args: extra_include_args,
                              required: false)
  if not vamp_dep.found()
    vamp_dep = cpp.find_library('vamp-sdk',
                                dirs: get_option('extra_lib_dirs'),
                                has_headers: ['vamp-sdk.h'],
                                header_args: extra_include_args,
                                required: false)
  endif
endif
have_vamp = vamp_dep.found()

if not sndfile_dep.found()
  sndfile_dep = cpp.find_library('sndfile',
                                 dirs: get_option('extra_lib_dirs'),
                                 has_headers: ['sndfile.h'],
                                 header_args: extra_include_args,
                                 required: false)
  if not sndfile_dep.found()
    sndfile_dep = cpp.find_library('sndfile-1',
                                   dirs: get_option('extra_lib_dirs'),
                                   has_headers: ['sndfile.h'],
                                   header_args: extra_include_args,
                                   required: false)
  endif
endif
have_sndfile = sndfile_dep.found()


# General platform and compiler expectations

ladspa_symbol_args = []
lv2_symbol_args = []
vamp_symbol_args = []

if get_option('buildtype').startswith('release')
  config_summary += { 'Build type': 'Release' }
  feature_defines += ['-DNO_THREAD_CHECKS', '-DNO_TIMING', '-DNDEBUG']
else
  config_summary += { 'Build type': 'Debug' }
endif

if system == 'darwin'
  feature_defines += ['-DUSE_PTHREADS', '-DMALLOC_IS_ALIGNED']
  ladspa_symbol_args += [
    '-exported_symbols_list', meson.current_source_dir() / 'ladspa-lv2/ladspa-plugin.list'
  ]
  lv2_symbol_args += [
    '-exported_symbols_list', meson.current_source_dir() / 'ladspa-lv2/lv2-plugin.list'
  ]
  vamp_symbol_args += [
    '-exported_symbols_list', meson.current_source_dir() / 'vamp/vamp-plugin.list'
  ]

  if architecture == 'aarch64'
    arch_flags += [
      '-arch', 'arm64',
    ]
  elif architecture == 'x86_64'
    arch_flags += [
      '-arch', 'x86_64',
    ]
  else # begin architecture != 'aarch64' or 'x86_64'
    error('Build for architecture ' + architecture + ' is not supported on this platform')
  endif # end architecture

  have_version_min = false
  foreach arg: get_option('cpp_args')
    if arg.contains('version-min')
      have_version_min = true
      bits = arg.split('=')
      if bits.length() > 1
	config_summary += { 'Target OS': bits[1] + '+' }
      else
	config_summary += { 'Target OS': '(unknown)' }
      endif
    endif
  endforeach

  if not have_version_min
    message('Using default minimum target OS version')
    message('(consider specifying this in cross-file if earlier target is desired)')
    if architecture == 'aarch64'
      arch_flags += [ '-mmacosx-version-min=11' ]
      config_summary += { 'Target OS': '11+' }
    else 
      arch_flags += [ '-mmacosx-version-min=10.13' ]
      config_summary += { 'Target OS': '10.13+' }
    endif
  endif

elif system == 'windows'
  feature_defines += ['-D_WIN32', '-DNOMINMAX', '-D_USE_MATH_DEFINES', '-DGETOPT_API=']
  if cpp.get_id() == 'msvc'
    ladspa_symbol_args += ['-EXPORT:ladspa_descriptor']
    lv2_symbol_args += ['-EXPORT:lv2_descriptor']
    vamp_symbol_args += ['-EXPORT:vampGetPluginDescriptor']
  endif

else # system not darwin or windows
  feature_defines += ['-DUSE_PTHREADS', '-DHAVE_POSIX_MEMALIGN']
  ladspa_symbol_args += [
    '-Wl,--version-script=' + meson.current_source_dir() / 'ladspa-lv2/ladspa-plugin.map'
  ]
  lv2_symbol_args += [
    '-Wl,--version-script=' + meson.current_source_dir() / 'ladspa-lv2/lv2-plugin.map'
  ]
  vamp_symbol_args += [
    '-Wl,--version-script=' + meson.current_source_dir() / 'vamp/vamp-plugin.map'
  ]
endif # system


general_include_dirs += get_option('extra_include_dirs')
general_compile_args = [ arch_flags, feature_defines ]
general_dependencies = [ feature_dependencies, thread_dep ]

rubberband_additional_static_lib = ''

if cpp.get_id() == 'msvc'
  #
  # In the MSVC world we have a quandary, partly as a result of
  # wanting to use naming compatible with our previous/other build
  # systems.
  #
  # Meson would like to use rubberband.dll for the dynamic library,
  # rubberband.lib for the import library, and librubberband.a for the
  # static library. This is kind of ok, even though lib*.a is not a
  # very familiar naming style here - except that previously we called
  # the static library rubberband-static.lib. (It would be usual to
  # expect some .lib file to be produced as a static library,
  # especially if default_library=static is set.)
  #
  # Our "solution" is to leave alone if default_library=shared (when
  # the Meson and MSVC ways are the same), but emit an additional
  # static .lib called rubberband-static.lib otherwise.
  #
  if get_option('default_library') != 'shared'
    rubberband_additional_static_lib = 'rubberband-static'
  endif    
  rubberband_library_name = 'rubberband'
  rubberband_program_name = 'rubberband-program'
  rubberband_ladspa_name = 'ladspa-rubberband'
  rubberband_lv2_name = 'lv2-rubberband'
  rubberband_vamp_name = 'vamp-rubberband'
  rubberband_jni_name = 'rubberband-jni'
else
  rubberband_library_name = 'rubberband'
  rubberband_dynamic_name = 'rubberband'
  rubberband_program_name = 'rubberband'
  rubberband_ladspa_name = 'ladspa-rubberband'
  rubberband_lv2_name = 'lv2-rubberband'
  rubberband_vamp_name = 'vamp-rubberband'
  rubberband_jni_name = 'rubberband-jni'
endif  

rubberband_objlib = static_library(
  'rubberband_objlib',
  library_sources,
  feature_sources,
  include_directories: general_include_dirs,
  cpp_args: general_compile_args,
  c_args: general_compile_args,
  dependencies: general_dependencies,
  pic: true,
  install: false,
)

rubberband_objlib_dep = declare_dependency(
  link_with: rubberband_objlib,
)


# And the build targets: Static and dynamic libraries, command-line
# utility, LADSPA and LV2 plugins, Vamp plugin, JNI library

if get_option('default_library') == 'shared'
  message('Not building Rubber Band Library static library: default_library option is set to shared')
  target_summary += { 'Static library': false }
else
  message('Will build Rubber Band Library static library')
  if rubberband_additional_static_lib != ''
    target_summary += { 'Static library': [ true, 'Name: ' + rubberband_additional_static_lib ] }
  else 
    target_summary += { 'Static library': [ true, 'Name: ' + rubberband_library_name ] }
  endif
endif

if get_option('default_library') == 'static'
  message('Not building Rubber Band Library dynamic library: default_library option is set to static')
  target_summary += { 'Dynamic library': false }
else
  message('Will build Rubber Band Library dynamic library')
  target_summary += { 'Dynamic library': [ true, 'Name: ' + rubberband_library_name ] }
endif
  
rubberband_library = library(
  rubberband_library_name,
  # We would like to write "link_with: rubberband_objlib",
  # but that didn't work with MSVC when I tried it: no
  # DLL entry point found
  objects: rubberband_objlib.extract_all_objects(recursive: true),
  link_args: [
    arch_flags,
    feature_libraries,
  ],
  dependencies: general_dependencies,
  version: rubberband_dynamic_library_version,
  install: true,
)

# This dependency is not used in this build file, but is declared
# for use when including this project as a subproject using Wrap
#
rubberband_dep = declare_dependency(
  link_with: rubberband_library,
  include_directories: '.',
)

if get_option('default_library') != 'shared' and rubberband_additional_static_lib != ''
  rubberband_additional_library = static_library(
    rubberband_additional_static_lib,
    link_with: rubberband_objlib,
    name_prefix: '',
    name_suffix: 'lib',
    install: true
  )
endif

if have_jni and javac.found() and jar.found()
  target_summary += { 'JNI library': [ true, 'Name: ' + rubberband_jni_name ] }
  message('Will build Java Native Interface')
  rubberband_jni = shared_library(
    rubberband_jni_name,
    jni_sources,
    include_directories: general_include_dirs,
    cpp_args: general_compile_args,
    c_args: general_compile_args,
    link_args: [
      arch_flags,
      feature_libraries,
    ],
    dependencies: [
      rubberband_objlib_dep,
      general_dependencies,
    ],
    # NB the JNI library is not versioned
    install: true,
  )
  rubberband_class = custom_target(
    'rubberband_class',
    input: 'com/breakfastquay/rubberband/RubberBandStretcher.java',
    output: 'RubberBandStretcher.class',
    command: [ javac, '@INPUT@', '-d', '@OUTDIR@' ],
  )
  rubberband_jar = custom_target(
    'rubberband_jar',
    input: rubberband_class,
    output: 'rubberband.jar',
    command: [ jar, 'cvf', '@OUTPUT@', 'com/breakfastquay/rubberband/@INPUT@' ],
    build_by_default: true,
  )
else
  target_summary += { 'JNI library': false }
  if not have_jni
    message('Not building Java Native Interface: jni.h header not found')
  else 
    message('Not building Java Native Interface: Java compiler not found')
  endif
endif

install_headers(
  [ 'rubberband/RubberBandStretcher.h',
    'rubberband/rubberband-c.h'
  ],
  subdir: 'rubberband'
)

if have_ladspa
  target_summary += { 'LADSPA plugin': [ true, 'Name: ' + rubberband_ladspa_name ] }
  message('Will build LADSPA plugin')
  rubberband_ladspa = shared_library(
    rubberband_ladspa_name,
    ladspa_sources,
    include_directories: general_include_dirs,
    cpp_args: general_compile_args,
    c_args: general_compile_args,
    link_args: [
      arch_flags,
      feature_libraries,
      ladspa_symbol_args,
    ],
    dependencies: [
      rubberband_objlib_dep,
      general_dependencies,
    ],
    name_prefix: '',
    install: true,
    install_dir: get_option('libdir') / 'ladspa',
  )
  install_data(
    'ladspa-lv2/ladspa-rubberband.cat',
    install_dir: get_option('libdir') / 'ladspa',
  )
  install_data(
    'ladspa-lv2/ladspa-rubberband.rdf',
    install_dir: get_option('datadir') / 'ladspa/rdf',
  )
else
  target_summary += { 'LADSPA plugin': false }
  message('Not building LADSPA plugin: ladspa.h header not found')
endif

if have_lv2
  target_summary += { 'LV2 plugin': [ true, 'Name: ' + rubberband_lv2_name ] }
  message('Will build LV2 plugin')
  rubberband_lv2 = shared_library(
    rubberband_lv2_name,
    lv2_sources,
    include_directories: general_include_dirs,
    cpp_args: general_compile_args,
    c_args: general_compile_args,
    link_args: [
      arch_flags,
      feature_libraries,
      lv2_symbol_args,
    ],
    dependencies: [
      rubberband_objlib_dep,
      general_dependencies,
    ],
    name_prefix: '',
    install: true,
    install_dir: get_option('libdir') / 'lv2/rubberband.lv2',
  )
  install_data(
    'ladspa-lv2/rubberband.lv2/manifest.ttl',
    'ladspa-lv2/rubberband.lv2/lv2-rubberband.ttl',
    install_dir: get_option('libdir') / 'lv2/rubberband.lv2',
  )
else
  target_summary += { 'LV2 plugin': false }
  message('Not building LV2 plugin: lv2.h header not found')
endif

if have_vamp
  target_summary += { 'Vamp plugin': [ true, 'Name: ' + rubberband_vamp_name ] }
  message('Will build Vamp plugin')
  rubberband_vamp = shared_library(
    rubberband_vamp_name,
    vamp_sources,
    include_directories: general_include_dirs,
    cpp_args: general_compile_args,
    c_args: general_compile_args,
    link_args: [
      arch_flags,
      feature_libraries,
      vamp_symbol_args,
    ],
    dependencies: [
      rubberband_objlib_dep,
      general_dependencies,
      vamp_dep,
    ],
    name_prefix: '',
    install: true,
    install_dir: get_option('libdir') / 'vamp',
  )
  install_data(
    'vamp/vamp-rubberband.cat',
    install_dir: get_option('libdir') / 'vamp',
  )
else
  target_summary += { 'Vamp plugin': false }
  message('Not building Vamp plugin: Vamp dependency not found')
endif

if have_sndfile
  target_summary += { 'Command-line utility': [ true, 'Name: ' + rubberband_program_name ] }
  message('Will build command-line utility')
  rubberband_program = executable(
    rubberband_program_name,
    program_sources,
    include_directories: general_include_dirs,
    cpp_args: general_compile_args,
    c_args: general_compile_args,
    link_args: [
      arch_flags,
      feature_libraries,
    ],
    dependencies: [
      rubberband_objlib_dep,
      general_dependencies,
      sndfile_dep,
    ],
    install: true,
  )
else 
  target_summary += { 'Command-line utility': false }
  message('Not building command-line utility: libsndfile dependency not found')
endif

pkg.generate(
  name: 'rubberband',
  description: 'Audio time-stretching and pitch-shifting library',
  url: 'https://breakfastquay.com/rubberband/',
  version: meson.project_version(),
  requires: pkgconfig_requirements,
  libraries: ['-L${libdir} -lrubberband'] + pkgconfig_libraries,
  extra_cflags: '-I${includedir}',
)

summary({'prefix': get_option('prefix'),
         'bindir': get_option('bindir'),
         'libdir': get_option('libdir'),
         'datadir': get_option('datadir'),
        }, section: 'Directories')

summary(config_summary + { 'Architecture': architecture },
        section: 'Configuration', bool_yn: true)
summary(target_summary, section: 'Build targets', bool_yn: true)

if system == 'darwin'
  foreach arg: get_option('cpp_args')
    if arg.contains('iPhone')
      summary({'Please note': 'You cannot legally distribute the Rubber Band Library\n             in an iOS app on the App Store, unless you have first obtained a\n             commercial licence.'}, section: '***')
      break
    endif
  endforeach
endif
