//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// SimpleAllocator - High-level manager of SimplePool
// memory pools. See SimplePool.h for description of allocation
// strategy
//
#ifndef SIMPLEALLOCATOR_H
#define SIMPLEALLOCATOR_H

#include "Allocator.h"
#include <list>
#include "SimplePool.h"

template<class T>
class SimpleAllocator : public Allocator
{
    typedef SimplePool<T> StorageType;
    typedef std::list<StorageType* > StorageList;

    public:
        SimpleAllocator() {}

        ~SimpleAllocator()
        {
            for(typename StorageList::iterator iter = m_pPoolList.begin(); iter != m_pPoolList.end(); ++iter)
            {
                delete *iter;
            }
            m_pPoolList.clear();
        }

        void* alloc()
        {
            if(m_pPoolList.empty() || m_pPoolList.back()->isFull())
            {
                // new storage must be allocated
                m_pPoolList.push_back(new StorageType);
            }

            // allocate from the last pool
            return m_pPoolList.back()->alloc();
        }

        void dealloc(void* ptr)
        {
            return;
            // deallocation not tracked in this strategy
        }

        // Reset the Allocator to hold one empty pool.
        // The pool is reset, without freeing any memory.
        void reset()
        {
            if (m_pPoolList.empty()) return;

            typename StorageList::iterator iter = m_pPoolList.begin();
            const typename StorageList::iterator E = m_pPoolList.end();

            // Reset the first pool
            StorageType * pFirstPool = *iter;
            pFirstPool->reset();
            iter++;

            // Erase all additional pools
            for( ; iter != E; iter++)
                delete *iter;
            m_pPoolList.resize(1);
        }

        // If the poolList is empty, add an empty pool
        void makePool()
        {
            if(m_pPoolList.empty())
                m_pPoolList.push_back(new StorageType);
        }

    private:

        StorageList m_pPoolList;
};

#endif
