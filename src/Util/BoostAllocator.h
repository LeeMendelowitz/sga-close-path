//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// BoostAllocator - Wrapper around Boost Memory Pool
//
#ifndef BOOSTALLOCATOR_H
#define BOOSTALLOCATOR_H

#include <boost/pool/object_pool.hpp>
#include "Allocator.h"

template<class T>
class BoostAllocator : public Allocator
{
    typedef boost::object_pool<T> Pool;

    public:
        BoostAllocator() {}

        ~BoostAllocator() { }

        void* alloc()
        {
            T * const ptr = m_pool.malloc();
            return static_cast<void *>(ptr);
        }

        void dealloc(void* ptr)
        {
            T * const p = static_cast<T *>(ptr);
            m_pool.free(p);
        }

        void reset() {};

    private:
        Pool m_pool;
};

#endif
