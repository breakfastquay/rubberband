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

#include "sysutils.h"

#ifdef _WIN32
#include <windows.h>
#include <fcntl.h>
#include <io.h>
#else /* !_WIN32 */
#include <signal.h>
#include <unistd.h>
#ifdef __APPLE__
#include <sys/sysctl.h>
#include <mach/mach.h>
#include <mach/mach_time.h>
#else /* !__APPLE__, !_WIN32 */
#include <stdio.h>
#include <string.h>
#endif /* !__APPLE__, !_WIN32 */
#endif /* !_WIN32 */

#ifdef __sun
#include <sys/processor.h>
#endif

#include <cstdlib>
#include <iostream>



#ifdef _WIN32
#include <fstream>
#endif


namespace RubberBand {

const char *
system_get_platform_tag()
{
#ifdef _WIN32
    return "win32";
#else /* !_WIN32 */
#ifdef __APPLE__
    return "osx";
#else 
#ifdef __LINUX__
    if (sizeof(long) == 8) {
        return "linux64";
    } else {
        return "linux";
    }
#else 
    return "posix";
#endif 
#endif 
#endif /* !_WIN32 */
}

bool
system_is_multiprocessor()
{
    static bool tested = false, mp = false;

    if (tested) return mp;
    int count = 0;

#ifdef _WIN32

    SYSTEM_INFO sysinfo;
    GetSystemInfo(&sysinfo);
    count = sysinfo.dwNumberOfProcessors;

#else /* !_WIN32 */
#ifdef __APPLE__
    
    size_t sz = sizeof(count);
    if (sysctlbyname("hw.ncpu", &count, &sz, NULL, 0)) {
        count = 0;
        mp = false;
    } else {
        mp = (count > 1);
    }

#else /* !__APPLE__, !_WIN32 */
#ifdef __sun

    processorid_t i, n;
    n = sysconf(_SC_CPUID_MAX);
    for (i = 0; i <= n; ++i) {
        int status = p_online(i, P_STATUS);
        if (status == P_ONLINE) {
            ++count;
        }
        if (count > 1) break;
    }

#else /* !__sun, !__APPLE__, !_WIN32 */

    //...

    FILE *cpuinfo = fopen("/proc/cpuinfo", "r");
    if (!cpuinfo) return false;

    char buf[256];
    while (!feof(cpuinfo)) {
        if (!fgets(buf, 256, cpuinfo)) break;
        if (!strncmp(buf, "processor", 9)) {
            ++count;
        }
        if (count > 1) break;
    }

    fclose(cpuinfo);

#endif /* !__sun, !__APPLE__, !_WIN32 */
#endif /* !__APPLE__, !_WIN32 */
#endif /* !_WIN32 */

    mp = (count > 1);
    tested = true;
    return mp;
}

#ifdef _WIN32

void gettimeofday(struct timeval *tv, void *tz)
{
    union { 
	long long ns100;  
	FILETIME ft; 
    } now; 
    
    ::GetSystemTimeAsFileTime(&now.ft); 
    tv->tv_usec = (long)((now.ns100 / 10LL) % 1000000LL); 
    tv->tv_sec = (long)((now.ns100 - 116444736000000000LL) / 10000000LL); 
}

void clock_gettime(int, struct timespec *ts)
{
    static LARGE_INTEGER cps;
    static bool haveCps = false;
    
    if (!haveCps) {
        QueryPerformanceFrequency(&cps);
        haveCps = true;
    }

    LARGE_INTEGER counter;
    QueryPerformanceCounter(&counter);

    //!!! check this
    ts->tv_sec = counter.QuadPart / cps.QuadPart;
    double sub = counter.QuadPart % cps.QuadPart;
    sub = sub / cps.QuadPart;
    sub = sub * 1000000000.;
    ts->tv_nsec = long(sub) ;
}

void usleep(unsigned long usec)
{
    ::Sleep(usec == 0 ? 0 : usec < 1000 ? 1 : usec / 1000);
}

#endif

#ifdef __APPLE__

void clock_gettime(int, struct timespec *ts)
{
    uint64_t t = mach_absolute_time();
    static mach_timebase_info_data_t sTimebaseInfo;
    if (sTimebaseInfo.denom == 0) (void)mach_timebase_info(&sTimebaseInfo);
    uint64_t n = t * sTimebaseInfo.numer / sTimebaseInfo.denom;
    ts->tv_sec = n / 1000000000;
    ts->tv_nsec = n % 1000000000;
}

#endif

void system_specific_initialise()
{
}

void system_specific_application_initialise()
{
}


ProcessStatus
system_get_process_status(int pid)
{
#ifdef _WIN32
    HANDLE handle = OpenProcess(PROCESS_QUERY_INFORMATION, FALSE, pid);
    if (!handle) {
        return ProcessNotRunning;
    } else {
        CloseHandle(handle);
        return ProcessRunning;
    }
#else
    if (kill(getpid(), 0) == 0) {
        if (kill(pid, 0) == 0) {
            return ProcessRunning;
        } else {
            return ProcessNotRunning;
        }
    } else {
        return UnknownProcessStatus;
    }
#endif
}

#ifdef _WIN32
void system_memorybarrier()
{
    LONG Barrier = 0;
    __asm__ __volatile__("xchgl %%eax,%0 "
                         : "=r" (Barrier));
}
#else /* !_WIN32 */
#if (__GNUC__ > 4) || (__GNUC__ == 4 && __GNUC_MINOR__ >= 1)
// Not required
#else
#include <pthread.h>
void system_memorybarrier()
{
    pthread_mutex_t dummy = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_lock(&dummy);
    pthread_mutex_unlock(&dummy);
}
#endif
#endif

}



