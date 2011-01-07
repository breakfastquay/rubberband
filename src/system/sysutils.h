/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band
    An audio time-stretching and pitch-shifting library.
    Copyright 2007-2011 Chris Cannam.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#ifndef _RUBBERBAND_SYSUTILS_H_
#define _RUBBERBAND_SYSUTILS_H_


#ifdef __GNUC__
#define R__ __restrict__
#endif

#ifndef R__
#define R__
#endif

#ifdef __MINGW32__
#include <malloc.h>
#else
#include <alloca.h>
#endif


#include <stdint.h>

#include <math.h>

namespace RubberBand {

extern const char *system_get_platform_tag();
extern bool system_is_multiprocessor();
extern void system_specific_initialise();
extern void system_specific_application_initialise();

#ifdef _WIN32

struct timeval { long tv_sec; long tv_usec; };
void gettimeofday(struct timeval *p, void *tz);

#endif


enum ProcessStatus { ProcessRunning, ProcessNotRunning, UnknownProcessStatus };
extern ProcessStatus GetProcessStatus(int pid);

inline double mod(double x, double y) { return x - (y * floor(x / y)); }
inline float modf(float x, float y) { return x - (y * float(floor(x / y))); }

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

inline double princarg(double a) { return mod(a + M_PI, -2.0 * M_PI) + M_PI; }
inline float princargf(float a) { return modf(a + (float)M_PI, -2.f * (float)M_PI) + (float)M_PI; }

} // end namespace

// The following should be functions in the RubberBand namespace, really

#ifdef _WIN32

#define MLOCK(a,b)   1
#define MUNLOCK(a,b) 1
#define MUNLOCK_SAMPLEBLOCK(a) 1

#define DLOPEN(a,b)  LoadLibrary((a).toStdWString().c_str())
#define DLSYM(a,b)   GetProcAddress((HINSTANCE)(a),(b))
#define DLCLOSE(a)   FreeLibrary((HINSTANCE)(a))
#define DLERROR()    ""

#else

#include <sys/mman.h>
#include <dlfcn.h>
#include <stdio.h>

#define MLOCK(a,b)   ::mlock((char *)(a),(b))
#define MUNLOCK(a,b) (::munlock((char *)(a),(b)) ? (::perror("munlock failed"), 0) : 0)
#define MUNLOCK_SAMPLEBLOCK(a) do { if (!(a).empty()) { const float &b = *(a).begin(); MUNLOCK(&b, (a).capacity() * sizeof(float)); } } while(0);

#define DLOPEN(a,b)  dlopen((a).toStdString().c_str(),(b))
#define DLSYM(a,b)   dlsym((a),(b))
#define DLCLOSE(a)   dlclose((a))
#define DLERROR()    dlerror()

#endif

#endif
