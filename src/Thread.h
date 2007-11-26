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

#ifndef _RUBBERBAND_THREAD_H_
#define _RUBBERBAND_THREAD_H_

#include <pthread.h>

namespace RubberBand
{

class Thread
{
public:
    typedef pthread_t Id;

    Thread();
    virtual ~Thread();

    Id id();

    void start();
    void wait();

    static bool threadingAvailable();

protected:
    virtual void run() = 0;

private:
    pthread_t m_id;
    bool m_extant;
    static void *staticRun(void *);
};

class Mutex
{
public:
    Mutex();
    ~Mutex();

    void lock();
    void unlock();
    bool trylock();

private:
    pthread_mutex_t m_mutex;
    bool m_locked;
};

class MutexLocker
{
public:
    MutexLocker(Mutex *);
    ~MutexLocker();

private:
    Mutex *m_mutex;
};

class Condition
{
public:
    Condition();
    ~Condition();
    
    void lock();
    void unlock();
    void wait(int us = 0);
    void signal();
    
private:
    pthread_mutex_t m_mutex;
    bool m_locked;
    pthread_cond_t m_condition;
};

}

#endif
