
Rubber Band
===========

An audio time-stretching and pitch-shifting library and utility program.

Written by Chris Cannam, chris.cannam@breakfastquay.com.
Copyright 2007-2018 Particular Programs Ltd.

Rubber Band is a library and utility program that permits changing the
tempo and pitch of an audio recording independently of one another.

See http://breakfastquay.com/rubberband/ for more information.


Licence
=======

Rubber Band is distributed under the GNU General Public License. See
the file COPYING for more information.

If you wish to distribute code using the Rubber Band Library under
terms other than those of the GNU General Public License, you must
obtain a commercial licence from us before doing so. In particular,
you may not legally distribute through any Apple App Store unless you
have a commercial licence.  See http://breakfastquay.com/rubberband/
for licence terms.

If you have obtained a valid commercial licence, your licence
supersedes this README and the enclosed COPYING file and you may
redistribute and/or modify Rubber Band under the terms described in
that licence. Please refer to your licence agreement for more details.

Note that Rubber Band may link with other GPL libraries or with
proprietary libraries, depending on its build configuration. See the
section "FFT and resampler selection" below. It is your responsibility
to ensure that you redistribute only in accordance with the licence
terms of any other libraries you may build with.


Contents of this README
-----------------------

1. Code components
2. Using the Rubber Band command-line tool
3. Using the Rubber Band Library
4. Compiling Rubber Band
   a. FFT and resampler selection
   b. Other supported #defines
   c. Building on Linux
   d. Building on macOS
   e. Building for iOS
   f. Building on Windows with Visual C++
   g. Building for Android and Java integration
5. Copyright notes for bundled libraries


1. Code components
------------------

Rubber Band consists of:

 * The Rubber Band library code.  This is the code that will normally
   be used by your applications.  The headers for this are in the
   rubberband/ directory, and the source code is in src/.
   The Rubber Band library depends upon resampler and FFT code; see
   section 3a below for details.

 * The Rubber Band command-line tool.  This is in main/main.cpp.
   This program uses the Rubber Band library and also requires libsndfile
   (http://www.mega-nerd.com/libsndfile/, licensed under the GNU Lesser
   General Public License) for audio file loading.

 * A pitch-shifter LADSPA audio effects plugin.  This is in ladspa/.
   It requires the LADSPA SDK header ladspa.h (not included).

 * A Vamp audio analysis plugin which may be used to inspect the
   dynamic stretch ratios and other decisions taken by the Rubber Band
   library when in use.  This is in vamp/.  It requires the Vamp
   plugin SDK (http://www.vamp-plugins.org/develop.html) (not included).


2. Using the Rubber Band command-line tool
------------------------------------------

The Rubber Band command-line tool builds as bin/rubberband.  The basic
incantation is

  $ rubberband -t <timeratio> -p <pitchratio> <infile.wav> <outfile.wav>

For example,

  $ rubberband -t 1.5 -p 2.0 test.wav output.wav

stretches the file test.wav to 50% longer than its original duration,
shifts it up in pitch by one octave, and writes the output to output.wav.

Several further options are available: run "rubberband -h" for help.
In particular, different types of music may benefit from different
"crispness" options (-c <n> where <n> is from 0 to 6).


3. Using the Rubber Band library
--------------------------------

The Rubber Band library has a public API that consists of one C++
class, called RubberBandStretcher in the RubberBand namespace.  You
should #include <rubberband/RubberBandStretcher.h> to use this class.
There is extensive documentation in the class header.

A header with C language bindings is also provided in
<rubberband/rubberband-c.h>.  This is a wrapper around the C++
implementation, and as the implementation is the same, it also
requires linkage against the C++ standard libraries.  It is not yet
documented separately from the C++ header.  You should include only
one of the two headers, not both.

The source code for the command-line utility (main/main.cpp) provides
a good example of how to use Rubber Band in offline mode; the LADSPA
pitch shifter plugin (ladspa/RubberBandPitchShifter.cpp) may be used
as an example of Rubber Band in real-time mode.

IMPORTANT: Please ensure you have read and understood the licensing
terms for Rubber Band before using it in your application.  This
library is provided under the GNU General Public License, which means
that any application that uses it must also be published under the GPL
or a compatible licence (i.e. with its full source code also available
for modification and redistribution) unless you have separately
acquired a commercial licence from the author.


4. Compiling Rubber Band
------------------------

4a. FFT and resampler selection
-------------------------------

Rubber Band requires additional library code for FFT calculation and
resampling.  Several libraries are supported.  The selection is
controlled using preprocessor flags at compile time, as detailed in
the tables below.

Flags that declare that you want to use an external library begin with
HAVE_; flags that select from the bundled options begin with USE_.

You must enable one resampler implementation and one FFT
implementation.  Do not enable more than one of either unless you know
what you're doing.

If you are building this software using one of the bundled library
options (Speex or KissFFT), please be sure to review the terms for
those libraries in src/speex/COPYING and src/kissfft/COPYING as
applicable.

FFT libraries supported
-----------------------

Name           Flags required        Notes
----           --------------        -----   

FFTW3	       -DHAVE_FFTW3	     GPL.

Accelerate     -DHAVE_VDSP	     Platform library on OS/X and iOS.

Intel IPP      -DHAVE_IPP            Proprietary library, can only be used with
      	    			     Rubber Band commercial licence. Define
				     USE_IPP_STATIC as well to build with static
				     IPP libraries.

KissFFT        -DUSE_KISSFFT	     Bundled, can be used with GPL or commercial
	    			     licence.  Single-precision. Slower than the
				     above options.

Resampler libraries supported
-----------------------------

Name           Flags required        Notes
----           --------------        -----   

libsamplerate  -DHAVE_LIBSAMPLERATE  GPL until v0.1.8, BSD for v0.1.9 and later.

libresample    -DHAVE_LIBRESAMPLE    LGPL.

Speex	       -DUSE_SPEEX	     Bundled, can be used with GPL or commercial
	       			     licence.


4b. Other supported #defines
----------------------------

Other symbols you may define at compile time are as follows. (Usually
the supplied build files will handle these for you.)

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
   double precision. Consider using for mobile architectures with
   slower double-precision support.

   -DUSE_POMMIER_MATHFUN
   Select the Julien Pommier implementations of trig functions for ARM
   NEON or x86 SSE architectures. These are usually faster but may be
   of lower precision than system implementations. Consider using this
   for mobile architectures.


4c. Building on Linux
---------------------

A GNU-style configure script is included for use on Linux and similar
systems.

Run ./configure, then adjust the generated Makefile according to your
preference for FFT and resampler implementations.  The default is to
use FFTW3 and libsamplerate.

The following Makefile targets are available:

  static  -- build static libraries only
  dynamic -- build dynamic libraries only
  library -- build static and dynamic libraries only
  program -- build the command-line tool
  vamp    -- build Vamp plugin
  ladspa  -- build LADSPA plugin
  all     -- build everything.

The default target is "all".


4d. Building on macOS
---------------------

A Makefile for macOS is provided as Makefile.osx.

Adjust the Makefile according to your preference for compiler and
platform SDK, FFT and resampler implementations.  The default is to
use the Accelerate framework and the Speex resampler.  Then run
e.g. "make -f Makefile.osx library" in a terminal window to build.
You will need the Xcode command-line tools installed.

(You probably don't want to use the configure script on macOS -- just
use Makefile.osx directly.)

The following Makefile targets are available:

  static  -- build static libraries only
  dynamic -- build dynamic libraries only
  library -- build static and dynamic libraries only
  program -- build the command-line tool
  vamp    -- build Vamp plugin
  ladspa  -- build LADSPA plugin
  all     -- build everything.

The default target is to build the static and dynamic libraries and
the command line tool.  The sndfile library is required for the
command line tool.

If you prefer to add the Rubber Band library files to an existing
build project instead of using the Makefile, the files in src/ (except
for RubberBandStretcherJNI.cpp) and the API headers in rubberband/
should be all you need.

Note that you cannot legally distribute applications using Rubber Band
in the Mac App Store, unless you have first obtained a commercial
licence for the Rubber Band Library.  GPL code is not permitted in the
app store.  See http://breakfastquay.com/technology/license.html for
commercial terms.


4e. Building for iOS
--------------------

A Makefile for iOS is provided as Makefile.ios.  It produces a single
static library containing both simulator and device binaries, in both
32- and 64-bit architectures.

Run e.g. "make -f Makefile.ios" in a terminal window to build.  You
will need the Xcode command-line tools installed.

If you prefer to add the Rubber Band library files to an existing
build project instead of using the Makefile, the files in src/ (except
for RubberBandStretcherJNI.cpp) and the API headers in rubberband/
should be all you need.

Note that you cannot legally distribute applications using Rubber Band
in the iOS App Store, unless you have a first obtained a commercial
licence for the Rubber Band Library. GPL code is not permitted in the
app store.  See http://breakfastquay.com/technology/license.html for
commercial terms.


4f. Building on Windows with Visual C++
---------------------------------------

A Visual Studio solution, targeted to VC 2015, with two projects is
supplied. The rubberband-library project builds the Rubber Band static
libraries only. The rubberband-program project builds the Rubber Band
command-line tool (which requires the Rubber Band libraries, and
libsndfile).

You will need to adjust the project settings so as to set the compile
flags according to your preference for FFT and resampler
implementation, and set the include path and library path
appropriately.  The default is to use the bundled KissFFT and the
Speex resampler.

If you prefer to add the Rubber Band library files to an existing
build project instead of using the supplied one, the files in src/
(except for RubberBandStretcherJNI.cpp) and the API headers in
rubberband/ should be all you need.


4g. Building for Android and Java integration
---------------------------------------------

An Android NDK build file is provided as Android.mk. This includes
compile definitions for a shared library built for ARM architectures
which can be loaded from a Java application using the Java native
interface (i.e. the Android NDK).

The Java side of the interface can be found in
com/breakfastquay/rubberband/RubberBandStretcher.java.

See
https://bitbucket.org/breakfastquay/rubberband-android-simple-sample
for a very trivial example of integration with Android Java code.

The supplied .mk file uses KissFFT and the Speex resampler.


5. Copyright notes for bundled libraries
========================================

5a. Speex
---------

[files in src/speex]

Copyright 2002-2007 	Xiph.org Foundation
Copyright 2002-2007 	Jean-Marc Valin
Copyright 2005-2007	Analog Devices Inc.
Copyright 2005-2007	Commonwealth Scientific and Industrial Research 
                        Organisation (CSIRO)
Copyright 1993, 2002, 2006 David Rowe
Copyright 2003 		EpicGames
Copyright 1992-1994	Jutta Degener, Carsten Bormann

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


5b. KissFFT
-----------

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


5c. Pommier math functions
--------------------------

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


5d. float_cast
--------------

[files in src/float_cast]

Copyright (C) 2001 Erik de Castro Lopo <erikd AT mega-nerd DOT com>

Permission to use, copy, modify, distribute, and sell this file for any 
purpose is hereby granted without fee, provided that the above copyright 
and this permission notice appear in all copies.  No representations are
made about the suitability of this software for any purpose.  It is 
provided "as is" without express or implied warranty.


5d. getopt
----------

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
