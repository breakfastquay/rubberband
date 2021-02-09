
project(
  'Rubber Band Library',
  'c', 'cpp',
  version: '1.9.0',
  license: 'GPL-2.0',
  default_options: [
    'cpp_std=c++11',
    'buildtype=release',
    'b_ndebug=if-release',
    'b_lundef=true',
  ],
  meson_version: '>= 0.53.0'
)

rubberband_dynamic_library_version = '2.1.2'

system = build_machine.system()
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
  'ladspa/RubberBandPitchShifter.cpp',
  'ladspa/libmain.cpp',
]

general_include_dirs = [
  'rubberband',
  'src',
]


# Scan for any dependencies we may use later; all are optional
# We collect these first so our summary messages about selected build
# features all appear together at the end

fftw3_dep = dependency('fftw3', version: '>= 3.0.0', required: false)
samplerate_dep = dependency('samplerate', version: '>= 0.1.8', required: false)
sndfile_dep = dependency('sndfile', version: '>= 1.0.16', required: false)
vamp_dep = dependency('vamp-sdk', version: '>= 2.9', required: false)
thread_dep = dependency('threads')
have_ladspa = cpp.has_header('ladspa.h')


# Check FFT and resampler options and set up dependencies and paths

feature_dependencies = []
feature_defines = []
feature_libraries = []
feature_sources = []
pkgconfig_requirements = []
arch_flags = []

resampler = get_option('resampler')
fft = get_option('fft')
ipp_path = get_option('ipp_path')
ipp_needed = false

if fft == 'auto'
  if system == 'darwin'
    fft = 'vdsp'
  else
    fft = 'kissfft'
  endif
endif

if resampler == 'auto'
  if samplerate_dep.found()
    resampler = 'libsamplerate'
  else
    resampler = 'speex'
  endif
endif

extra_include_args = []
foreach d: get_option('extra_include_dirs')
  extra_include_args += [ '-I' + d ]
endforeach

if fft == 'kissfft'
  message('For FFT: using KissFFT')
  if fftw3_dep.found()
    message('(to use FFTW instead, reconfigure with -Dfft=fftw)')
  endif
  feature_sources += ['src/kissfft/kiss_fft.c', 'src/kissfft/kiss_fftr.c']
  feature_defines += ['-DUSE_KISSFFT']

elif fft == 'fftw'
  if fftw3_dep.found()
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
  message('For FFT: using vDSP')
  feature_defines += ['-DHAVE_VDSP']
  feature_libraries += ['-framework', 'Accelerate']
  
elif fft == 'ipp'
  if ipp_path != ''
    message('For FFT: using IPP')
    message('IPP path defined as ' + ipp_path)
  else 
    error('For FFT: IPP selected, but ipp_path not specified')
  endif
  ipp_needed = true

else
  error('Unknown or unsupported FFT option: ' + fft)

endif # fft

if resampler == 'libsamplerate'
  if samplerate_dep.found()
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
  message('For resampler: using Speex')
  message('(consider libsamplerate if time-varying pitch shift is required)')
  feature_sources += ['src/speex/resample.c']
  feature_defines += ['-DUSE_SPEEX']
  
elif resampler == 'ipp'
  if ipp_path != ''
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
endif
have_sndfile = sndfile_dep.found()


# General platform and compiler expectations

ladspa_symbol_args = []
vamp_symbol_args = []

if get_option('buildtype').startswith('release')
  feature_defines += ['-DNO_THREAD_CHECKS', '-DNO_TIMING', '-DNDEBUG']
endif

if system == 'darwin'
  feature_defines += ['-DUSE_PTHREADS', '-DMALLOC_IS_ALIGNED']
  ladspa_symbol_args += [
    '-exported_symbols_list', meson.source_root() / 'ladspa/ladspa-plugin.list'
  ]
  vamp_symbol_args += [
    '-exported_symbols_list', meson.source_root() / 'vamp/vamp-plugin.list'
  ]
  if architecture == 'aarch64'
    arch_flags = [
      '-mmacosx-version-min=11',
      '-arch', 'arm64'
    ]
  elif architecture == 'x86_64'
    arch_flags = [
      '-mmacosx-version-min=10.11',
      '-arch', 'x86_64',
    ]
  endif # architecture

elif system == 'windows'
  if cpp.get_id() == 'msvc'
    feature_defines += ['-DWIN32', '-D__MSVC__', '-DNOMINMAX', '-D_USE_MATH_DEFINES', '-DGETOPT_API=']
    ladspa_symbol_args += ['-EXPORT:ladspa_descriptor']
    vamp_symbol_args += ['-EXPORT:vampGetPluginDescriptor']
  endif

else # system not darwin or windows
  feature_defines += ['-DUSE_PTHREADS', '-DHAVE_POSIX_MEMALIGN']
  ladspa_symbol_args += [
    '-Wl,--version-script=' + meson.source_root() / 'ladspa/ladspa-plugin.map'
  ]
  vamp_symbol_args += [
    '-Wl,--version-script=' + meson.source_root() / 'vamp/vamp-plugin.map'
  ]
endif # system


general_include_dirs += get_option('extra_include_dirs')
general_compile_args = [ arch_flags, feature_defines ]
general_dependencies = [ feature_dependencies, thread_dep ]

if system == 'windows'
  if get_option('no_shared')
    rubberband_static_name = 'rubberband'
  else
    rubberband_static_name = 'rubberband-static'
  endif    
  rubberband_dynamic_name = 'rubberband'
  rubberband_program_name = 'rubberband-program'
  # Meson likes libxxx.a even on Windows, and so might we for a
  # new library, but not here
  platform_static_name_prefix = ''
  platform_static_name_suffix = 'lib'
else
  rubberband_static_name = 'rubberband'
  rubberband_dynamic_name = 'rubberband'
  rubberband_program_name = 'rubberband'
  platform_static_name_prefix = 'lib'
  platform_static_name_suffix = 'a'
endif  


# And the build targets: Static and dynamic libraries, command-line
# utility, LADSPA plugin, Vamp plugin

message('Will build Rubber Band Library static library')
rubberband_static = static_library(
  rubberband_static_name,
  library_sources,
  feature_sources,
  include_directories: general_include_dirs,
  cpp_args: general_compile_args,
  c_args: general_compile_args,
  dependencies: general_dependencies,
  name_prefix: platform_static_name_prefix,
  name_suffix: platform_static_name_suffix,
  pic: true,
  install: true,
)

rubberband_static_dep = declare_dependency(
  link_with: rubberband_static,
)

if not get_option('no_shared')
  message('Will build Rubber Band Library dynamic library')
  rubberband_dynamic = shared_library(
    rubberband_dynamic_name,
    objects: rubberband_static.extract_all_objects(),
    link_args: [
      arch_flags,
      feature_libraries,
    ],
    dependencies: general_dependencies,
    version: rubberband_dynamic_library_version,
    install: true,
  )
else
  message('Not building Rubber Band Library dynamic library: no_shared option set')
endif

if have_ladspa
  message('Will build LADSPA plugin')
  rubberband_ladspa = shared_library(
    'ladspa-rubberband',
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
      rubberband_static_dep,
      general_dependencies,
    ],
    name_prefix: '',
    install: true,
    install_dir: get_option('libdir') / 'ladspa',
  )
  install_data(
    'ladspa/ladspa-rubberband.cat',
    install_dir: get_option('libdir') / 'ladspa',
  )
  install_data(
    'ladspa/ladspa-rubberband.rdf',
    install_dir: get_option('datadir') / 'ladspa/rdf',
  )
else
  message('Not building LADSPA plugin: ladspa.h header not found')
endif

if have_vamp
  message('Will build Vamp plugin')
  rubberband_vamp = shared_library(
    'vamp-rubberband',
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
      rubberband_static_dep,
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
  message('Not building Vamp plugin: Vamp dependency not found')
endif

if have_sndfile
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
      rubberband_static_dep,
      general_dependencies,
      sndfile_dep,
    ],
    install: true,
  )
else 
  message('Not building command-line utility: libsndfile dependency not found')
endif

pkg.generate(
  name: 'rubberband',
  description: 'Audio time-stretching and pitch-shifting library',
  url: 'https://breakfastquay.com/rubberband/',
  version: meson.project_version(),
  requires: pkgconfig_requirements,
  libraries: '-L${libdir} -lrubberband',
  extra_cflags: '-I${includedir}',
)
