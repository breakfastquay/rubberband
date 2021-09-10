
# Rubber Band

An audio time-stretching and pitch-shifting library and utility program.

Written by Chris Cannam, chris.cannam@breakfastquay.com.
Published by Particular Programs Ltd t/a Breakfast Quay.
Copyright 2007-2021 Particular Programs Ltd.

Rubber Band is a library and utility program that permits changing the
tempo and pitch of an audio recording independently of one another.

* About Rubber Band: https://breakfastquay.com/rubberband/
* Code repository: https://hg.sr.ht/~breakfastquay/rubberband
* Issue tracker: https://todo.sr.ht/~breakfastquay/rubberband
* Github mirror: https://github.com/breakfastquay/rubberband

CI builds:

* [![Build status](https://builds.sr.ht/~breakfastquay/rubberband.svg)](https://builds.sr.ht/~breakfastquay/rubberband?) (Linux)
* [![Build Status](https://github.com/breakfastquay/rubberband/workflows/macOS%20and%20iOS%20CI/badge.svg)](https://github.com/breakfastquay/rubberband/actions?query=workflow%3A%22macOS+and+iOS+CI%22) (macOS, iOS)
* [![Build Status](https://ci.appveyor.com/api/projects/status/hhhhpf718jwhpyf6?svg=true)](https://ci.appveyor.com/project/breakfastquay/rubberband) (Windows)


## Licence

Rubber Band Library is distributed under the GNU General Public
License (GPL). You can redistribute it and/or modify it under the
terms of the GPL; either version 2 of the License, or (at your option)
any later version. See the file COPYING for more information.

If you wish to distribute code using the Rubber Band Library under
terms other than those of the GNU General Public License, you must
obtain a commercial licence from us before doing so. In particular,
you may not legally distribute through any Apple App Store unless you
have a commercial licence.  See https://breakfastquay.com/rubberband/
for licence terms.

If you have obtained a valid commercial licence, your licence
supersedes this README and the enclosed COPYING file and you may
redistribute and/or modify Rubber Band under the terms described in
that licence. Please refer to your licence agreement for more details.

Rubber Band includes a .NET interface generously contributed by
Jonathan Gilbert under a BSD-like licence. The files in the
dotnet/rubberband-dll and dotnet/rubberband-sharp directories fall
under this licence. If you make use of this interface, please ensure
you comply with the terms of its licence.

Rubber Band may link with other libraries or be compiled with other
source code, depending on its build configuration. It is your
responsibility to ensure that you redistribute these only in
accordance with their own licence terms, regardless of the conditions
under which you are redistributing the Rubber Band code itself. The
licences for some relevant library code are as follows, to the best of
our knowledge. See also the end of this README for detailed terms.

 * FFTW3 - GPL; proprietary licence needed for redistribution
 * Intel IPP - Proprietary; licence needed for redistribution
 * KissFFT - BSD-like
 * libsamplerate - BSD-like from version 0.1.9 onwards
 * libresample - LGPL
 * Speex - BSD-like
 * Pommier math functions - BSD-like
 

## Contents of this README

1. Code components
2. Using the Rubber Band command-line tool
3. Using the Rubber Band Library
4. Compiling Rubber Band
    a. Building on Linux
    b. Building on macOS
    c. Building for iOS
    d. Building on Windows
    e. Building for Android and Java integration
    f. FFT and resampler selection
    g. Other supported #defines
5. Copyright notes for bundled libraries


## 1. Code components

Rubber Band consists of:

 * The Rubber Band Library code.  This is the code that will normally
   be used by your applications.  The headers for this are in the
   rubberband/ directory, and the source code is in src/.
   The Rubber Band Library may also depend upon external resampler
   and FFT code; see section 3a below for details.

 * The Rubber Band command-line tool.  This is in main/main.cpp.
   This program uses the Rubber Band Library and also requires libsndfile
   (http://www.mega-nerd.com/libsndfile/, licensed under the GNU Lesser
   General Public License) for audio file loading.

 * A pitch-shifter LADSPA audio effects plugin.  This is in ladspa/.
   It requires the LADSPA SDK header ladspa.h (not included).

 * A Vamp audio analysis plugin which may be used to inspect the
   dynamic stretch ratios and other decisions taken by the Rubber Band
   Library when in use.  This is in vamp/.  It requires the Vamp
   plugin SDK (https://www.vamp-plugins.org/develop.html) (not included).


## 2. Using the Rubber Band command-line tool

The Rubber Band command-line tool builds as bin/rubberband.  The basic
incantation is

```
  $ rubberband -t <timeratio> -p <pitchratio> <infile.wav> <outfile.wav>
```

For example,

```
  $ rubberband -t 1.5 -p 2.0 test.wav output.wav
```

stretches the file test.wav to 50% longer than its original duration,
shifts it up in pitch by one octave, and writes the output to output.wav.

Several further options are available: run "rubberband -h" for help.
In particular, different types of music may benefit from different
"crispness" options (-c flag with a numerical argument from 0 to 6).


## 3. Using the Rubber Band Library

The Rubber Band Library has a public API that consists of one C++
class, called RubberBandStretcher in the RubberBand namespace.  You
should `#include <rubberband/RubberBandStretcher.h>` to use this
class.  There is extensive documentation in the class header.

A header with C language bindings is also provided in
`<rubberband/rubberband-c.h>`.  This is a wrapper around the C++
implementation, and as the implementation is the same, it also
requires linkage against the C++ standard libraries.  It is not yet
documented separately from the C++ header.  You should include only
one of the two headers, not both.

A .NET interface is also included, contributed by Jonathan Gilbert;
see the files in the `dotnet` directory for details.

The source code for the command-line utility (`main/main.cpp`)
provides a good example of how to use Rubber Band in offline mode; the
LADSPA pitch shifter plugin (`ladspa/RubberBandPitchShifter.cpp`) may
be used as an example of Rubber Band in real-time mode.

IMPORTANT: Please ensure you have read and understood the licensing
terms for Rubber Band before using it in your application.  This
library is provided under the GNU General Public License, which means
that any application that uses it must also be published under the GPL
or a compatible licence (i.e. with its full source code also available
for modification and redistribution) unless you have separately
acquired a commercial licence from the author.


## 4. Compiling the Rubber Band Library

The primary supported build system for the Rubber Band Library on all
platforms is Meson (https://mesonbuild.com). The Meson build system
can be used to build all targets (static and dynamic library,
command-line utility, and plugins) and to cross-compile.

If you only need a static library and don't wish to use Meson, some
alternative build files (Makefiles and Visual C++ projects) are
included in the `otherbuilds` directory. See the platform-specific
build sections below for more details.

To build with Meson, ensure Meson and Ninja are installed and run:

```
$ meson build && ninja -C build
```

This checks for necessary dependencies, reports what it finds, and if
all is well, builds the code into a subdirectory called `build`. It
will build everything it can find the requisite dependencies for:
static and dynamic libraries, LADSPA and Vamp plugins, and
command-line utility.

Some configuration options are provided, described in the
`meson_options.txt` file. To set one of these, add a `-D` option to
Meson:

```
$ meson build -Dipp_path=/opt/intel/ipp
```

The options are documented in the library- and platform-specific
sections below.

The Rubber Band Library is written entirely in C++ to the C++98
standard. It is unlikely to make any difference (performance or
otherwise) which C++ standard your compiler uses - as long as it's no
older than C++98!

If you are building this software using either of the Speex or KissFFT
library options, please be sure to review the terms for those
libraries in `src/speex/COPYING` and `src/kissfft/COPYING` as
applicable.


### 4a. Building on Linux

For best results, and to ensure the command-line tool and plugins are
built, first install libsamplerate, libsndfile, and the LADSPA and
Vamp plugin headers so they can be found using `pkg-config`. Then

```
$ meson build && ninja -C build
```

See "FFT and resampler selection" below for further build options.

Alternatively, if you only need the static library and prefer a
Makefile, try

```
$ make -f otherbuilds/Makefile.linux
```


### 4b. Building on macOS

Ensure the Xcode command-line tools are installed, and ideally also
install libsamplerate and libsndfile.

To build for the default architecture:

```
$ meson build && ninja -C build
```

Which architecture is the default may depend on the version of Meson
and/or the current shell. To force a particular architecture you can
use a Meson cross-file, as follows.

To build for Apple Silicon (arm64):

```
$ meson build --cross-file cross/macos-arm64.txt && ninja -C build
```

To build for Intel (x86_64):

```
$ meson build --cross-file cross/macos-x86_64.txt && ninja -C build
```

You can build a universal binary library for both architectures like
this:

```
$ meson build --cross-file cross/macos-universal.txt && ninja -C build
```

Note that the universal cross file also sets the minimum OS version to
the earliest supported macOS versions for both architectures. (Note
that actual compatibility will also depend on how any dependent
libraries have been compiled.)  You can edit this in the
`cross/macos-universal.txt` file if you want a specific target.

See "FFT and resampler selection" below for further build options.

Note that you cannot legally distribute applications using Rubber Band
in the Mac App Store, unless you have first obtained a commercial
licence for the Rubber Band Library.  GPL code is not permitted in the
app store.  See https://breakfastquay.com/technology/license.html for
commercial terms.


### 4c. Building for iOS

Ensure the Xcode command-line tools are installed, and

```
$ meson build_ios --cross-file cross/ios.txt && ninja -C build_ios
```

The output files will be found in the `build_ios` directory.

To build for the simulator,

```
$ meson build_sim --cross-file cross/ios-simulator.txt && ninja -C build_sim
```

The output files will be found in the `build_sim` directory.

See "FFT and resampler selection" below for further build options.

Note that you cannot legally distribute applications using Rubber Band
in the iOS App Store, unless you have a first obtained a commercial
licence for the Rubber Band Library. GPL code is not permitted in the
app store.  See https://breakfastquay.com/technology/license.html for
commercial terms.


### 4d. Building on Windows

If you only need to build the static library for integration into your
project, and you prefer a Visual Studio project file, you can find a
simple one in `otherbuilds\rubberband-library.vcxproj`.

The rest of this section describes the "full" build system, which uses
Meson just as on the other platforms. So to build this way, start by
ensuring Meson and Ninja are installed and available. Then, in a
terminal window with the compiler tools available in the path (e.g. a
Visual Studio command-line prompt for the relevant build architecture)
run

```
> meson build
> ninja -C build
```

The output files will be found in the `build` directory.

The Rubber Band code is compatible with both the traditional Visual
C++ compiler (`cl`) and the Clang front-end (`clang`), and the build
system will use whichever appears (first) in your path.

To build against a specific Visual C++ runtime, use the built-in Meson
option `b_vscrt`:

```
> meson build -Db_vscrt=mt
```

See "FFT and resampler selection" below for further build options.


### 4e. Building for Android and Java integration

Currently only a very old Android NDK build file is provided, as
`otherbuilds/Android.mk`. This includes compile definitions for a
shared library built for ARM architectures which can be loaded from a
Java application using the Java native interface (i.e. the Android
NDK).

The Java side of the interface can be found in
`com/breakfastquay/rubberband/RubberBandStretcher.java`.

See
https://hg.sr.ht/~breakfastquay/rubberband-android-simple-sample
for a very trivial example of integration with Android Java code.

The supplied `.mk` file uses KissFFT and the Speex resampler.


### 4f. FFT and resampler selection

Rubber Band requires the selection of library code for FFT calculation
and resampling.  Several libraries are supported.  The selection is
controlled (in Meson) using `-D` options and (in the code itself)
using preprocessor flags set by the build system. These options and
flags are detailed in the tables below.

At least one resampler implementation and one FFT implementation must
be enabled. It is technically possible to enable more than one, but
it's confusing and not often useful.

If you are building this software using the bundled Speex or KissFFT
library code, please be sure to review the terms for those libraries
in `src/speex/COPYING` and `src/kissfft/COPYING` as applicable.

If you are proposing to package Rubber Band for a Linux distribution,
please select either the built-in FFT (simpler for you) or FFTW (a bit
faster) and use libsamplerate.

#### FFT libraries supported

```
Library     Build option    CPP define     Notes
----        ------------    ----------     -----

Built-in    -Dfft=builtin   -DUSE_BUILTIN_FFT
                                           Default except on macOS/iOS.
                                           Can be distributed with either
                                           the Rubber Band GPL or
                                           commercial licence.

KissFFT     -Dfft=kissfft   -DHAVE_KISSFFT
                                           Single precision.
                                           Only indicated for use with
                                           single-precision sample type
                                           (see below).
                                           Bundled, can be distributed with
                                           either the Rubber Band GPL or
                                           commercial licence.

Accelerate  -Dfft=vdsp      -DHAVE_VDSP    Default on macOS/iOS.
                                           Best option on these platforms.

FFTW3       -Dfft=fftw      -DHAVE_FFTW3   GPL.

Intel IPP   -Dfft=ipp       -DHAVE_IPP     Proprietary, can only be used with
                                           Rubber Band commercial licence.
```

#### Resampler libraries supported

```
Library     Build option    CPP define     Notes
----        ------------    ----------     -----

libsamplerate               -DHAVE_LIBSAMPLERATE
            -Dresampler=libsamplerate      Best choice in most cases.

Speex                       -DUSE_SPEEX
            -Dresampler=speex              Bundled, can be distributed with
	    				   either the Rubber Band GPL or
					   commercial licence.
```

### 4g. Other supported #defines

Other known preprocessor symbols are as follows. (Usually the supplied
build files will handle these for you.)

    -DLACK_BAD_ALLOC
    Define on systems lacking std::bad_alloc in the C++ library.

    -DLACK_POSIX_MEMALIGN
    Define on systems lacking posix_memalign.

    -DUSE_OWN_ALIGNED_MALLOC
    Define on systems lacking any aligned malloc implementation.

    -DLACK_SINCOS
    Define on systems lacking sincos().
    
    -DNO_EXCEPTIONS
    Build without use of C++ exceptions.

    -DNO_THREADING
    Build without any multithread support.

    -DUSE_PTHREADS
    Use the pthreads library (required unless NO_THREADING or on Windows)

    -DPROCESS_SAMPLE_TYPE=float
    Select single precision for internal calculations. The default is
    double precision. Consider in conjunction with single-precision
    KissFFT for mobile architectures with slower double-precision
    support.

    -DUSE_POMMIER_MATHFUN
    Select the Julien Pommier implementations of trig functions for ARM
    NEON or x86 SSE architectures. These are usually faster but may be
    of lower precision than system implementations. Consider using this
    for mobile architectures.


## 5. Copyright notes for bundled libraries

### 5a. Speex

```
[files in src/speex]

Copyright 2002-2007     Xiph.org Foundation
Copyright 2002-2007     Jean-Marc Valin
Copyright 2005-2007     Analog Devices Inc.
Copyright 2005-2007     Commonwealth Scientific and Industrial Research 
                        Organisation (CSIRO)
Copyright 1993, 2002, 2006 David Rowe
Copyright 2003          EpicGames
Copyright 1992-1994     Jutta Degener, Carsten Bormann

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

- Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in the
documentation and/or other materials provided with the distribution.

- Neither the name of the Xiph.org Foundation nor the names of its
contributors may be used to endorse or promote products derived from
this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

### 5b. KissFFT

```
[files in src/kissfft]

Copyright (c) 2003-2004 Mark Borgerding

All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.
    * Neither the author nor the names of any contributors may be used
      to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
```

### 5c. Pommier math functions

```
[files in src/pommier]

Copyright (C) 2011  Julien Pommier

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.
```

### 5d. float_cast

```
[files in src/float_cast]

Copyright (C) 2001 Erik de Castro Lopo <erikd AT mega-nerd DOT com>

Permission to use, copy, modify, distribute, and sell this file for any 
purpose is hereby granted without fee, provided that the above copyright 
and this permission notice appear in all copies.  No representations are
made about the suitability of this software for any purpose.  It is 
provided "as is" without express or implied warranty.
```

### 5e. getopt

```
[files in src/getopt, used by command-line tool on some platforms]

Copyright (c) 2000 The NetBSD Foundation, Inc.
All rights reserved.

This code is derived from software contributed to The NetBSD Foundation
by Dieter Baron and Thomas Klausner.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:
1. Redistributions of source code must retain the above copyright
   notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
   notice, this list of conditions and the following disclaimer in the
   documentation and/or other materials provided with the distribution.
3. All advertising materials mentioning features or use of this software
   must display the following acknowledgement:
       This product includes software developed by the NetBSD
       Foundation, Inc. and its contributors.
4. Neither the name of The NetBSD Foundation nor the names of its
   contributors may be used to endorse or promote products derived
   from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE NETBSD FOUNDATION, INC. AND CONTRIBUTORS
``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE FOUNDATION OR CONTRIBUTORS
BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
POSSIBILITY OF SUCH DAMAGE.
```

### 5f. rubberband-sharp

```
[files in rubberband-dll and rubberband-sharp]

Copyright 2018-2019 Jonathan Gilbert

Permission is hereby granted, free of charge, to any person
obtaining a copy of this software and associated documentation
files (the "Software"), to deal in the Software without
restriction, including without limitation the rights to use, copy,
modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be
included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS BE LIABLE FOR
ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Except as contained in this notice, the name of Jonathan Gilbert
shall not be used in advertising or otherwise to promote the sale,
use or other dealings in this Software without prior written
authorization.
```
