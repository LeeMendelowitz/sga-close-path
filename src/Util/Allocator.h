//-----------------------------------------------
// Copyright 2010 Wellcome Trust Sanger Institute
// Written by Jared Simpson (js18@sanger.ac.uk)
// Released under the GPL license
//-----------------------------------------------
//
// Allocator interface

#ifndef ALLOCATOR_H
#define ALLOCATOR_H


class Allocator
{
    public:
        Allocator() {}
        virtual ~Allocator() {}

        // Allocate some memory
        virtual void* alloc() = 0;

        // Deallocate memory
        virtual void dealloc(void *) {};

        // Reset the allocator
        virtual void reset() = 0;
};

#endif
