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

#include "Thread.h"

#include <iostream>

#include <sys/time.h>
#include <time.h>

using std::cerr;
using std::endl;

namespace RubberBand
{


Thread::Thread() :
    m_id(0),
    m_extant(false)
{
}

Thread::~Thread()
{
    if (m_extant) {
        pthread_join(m_id, 0);
    }
}

void
Thread::start()
{
    if (pthread_create(&m_id, 0, staticRun, this)) {
        //!!!
        cerr << "ERROR: thread creation failed" << endl;
        exit(1);
    } else {
        m_extant = true;
    }
}    

void 
Thread::wait()
{
    if (m_extant) {
        pthread_join(m_id, 0);
        m_extant = false;
    }
}

Thread::Id
Thread::id()
{
    return m_id;
}

bool
Thread::threadingAvailable()
{
    return true;
}

void *
Thread::staticRun(void *arg)
{
    Thread *thread = (Thread *)arg;
    thread->run();
    return 0;
}

Mutex::Mutex()
{
    pthread_mutex_init(&m_mutex, 0);
}

Mutex::~Mutex()
{
    pthread_mutex_destroy(&m_mutex);
}

void
Mutex::lock()
{
    pthread_mutex_lock(&m_mutex);
}

void
Mutex::unlock()
{
    pthread_mutex_unlock(&m_mutex);
}

bool
Mutex::trylock()
{
    if (pthread_mutex_trylock(&m_mutex)) {
        return false;
    } else {
        return true;
    }
}

Condition::Condition()
{
    pthread_mutex_init(&m_mutex, 0);
    m_locked = false;
    pthread_cond_init(&m_condition, 0);
}

Condition::~Condition()
{
    if (m_locked) pthread_mutex_unlock(&m_mutex);
    pthread_cond_destroy(&m_condition);
    pthread_mutex_destroy(&m_mutex);
}

void
Condition::lock()
{
    if (m_locked) return;
    pthread_mutex_lock(&m_mutex);
    m_locked = true;
}

void 
Condition::wait(int us)
{
    lock();
    if (us == 0) {

        pthread_cond_wait(&m_condition, &m_mutex);

    } else {

        struct timeval now;
        gettimeofday(&now, 0);

        now.tv_usec += us;
        while (now.tv_usec > 1000000) {
            now.tv_usec -= 1000000;
            ++now.tv_sec;
        }

        struct timespec timeout;
        timeout.tv_sec = now.tv_sec;
        timeout.tv_nsec = now.tv_usec * 1000;
    
        pthread_cond_timedwait(&m_condition, &m_mutex, &timeout);
    }

    pthread_mutex_unlock(&m_mutex);
    m_locked = false;
}

void
Condition::signal()
{
    pthread_cond_signal(&m_condition);
}


MutexLocker::MutexLocker(Mutex *mutex) :
    m_mutex(mutex)
{
    if (m_mutex) {
        m_mutex->lock();
    }
}

MutexLocker::~MutexLocker()
{
    if (m_mutex) {
        m_mutex->unlock();
    }
}

}

