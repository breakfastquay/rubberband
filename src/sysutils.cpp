/* -*- c-basic-offset: 4 indent-tabs-mode: nil -*-  vi:set ts=8 sts=4 sw=4: */

/*
    Rubber Band
    An audio time-stretching and pitch-shifting library.
    Copyright 2007 Chris Cannam.
    
    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 2 of the
    License, or (at your option) any later version.  See the file
    COPYING included with this distribution for more information.
*/

#include "sysutils.h"

#ifdef _WIN32
#else /* !_WIN32 */
#ifdef __APPLE__
#include <sys/sysctl.h>
#else /* !__APPLE__, !_WIN32 */
#include <stdio.h>
#include <string.h>
#endif /* !__APPLE__, !_WIN32 */
#endif /* !_WIN32 */

namespace RubberBand {

bool
system_is_multiprocessor()
{
    static bool tested = false, mp = false;

    if (tested) return mp;
    int count = 0;

#ifdef _WIN32

    //...

#else /* !_WIN32 */
#ifdef __APPLE__
    
    size_t sz = sizeof(count);
    if (sysctlbyname("hw.ncpu", &count, &sz, NULL, 0)) {
        mp = false;
    } else {
        mp = (count > 1);
    }

#else /* !__APPLE__, !_WIN32 */
    
    //...

    FILE *cpuinfo = fopen("/proc/cpuinfo", "r");
    if (!cpuinfo) return false;

    char buf[256];
    while (!feof(cpuinfo)) {
        fgets(buf, 256, cpuinfo);
        if (!strncmp(buf, "processor", 9)) {
            ++count;
        }
        if (count > 1) break;
    }

    fclose(cpuinfo);

#endif /* !__APPLE__, !_WIN32 */
#endif /* !_WIN32 */

    mp = (count > 1);
    tested = true;
    return mp;
}

}



