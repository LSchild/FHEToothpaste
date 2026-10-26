//
// Created by leonard on 3/23/26.
//

#ifndef INTERFACES_H
#define INTERFACES_H

#include <vector>
#include <cstdint>
#include "hexl/hexl.hpp"
#include "static/container_types.h"

static const size_t ALIGN_TO = 64;

template <class T, std::size_t Align = ALIGN_TO>
struct AlignedAllocator {
    using value_type = T;
    AlignedAllocator() noexcept = default;
    template <class U>
    AlignedAllocator(const AlignedAllocator<U, Align>&) noexcept {}
    template <class U>
    struct rebind {
        using other = AlignedAllocator<U, Align>;
    };
    T* allocate(std::size_t n) {
        return static_cast<T*>(::operator new(n * sizeof(T), std::align_val_t(Align)));
    }
    void deallocate(T* p, std::size_t) noexcept { ::operator delete(p, std::align_val_t(Align)); }
    friend bool operator==(const AlignedAllocator&, const AlignedAllocator&) { return true; }
    friend bool operator!=(const AlignedAllocator&, const AlignedAllocator&) { return false; }
};

using AlignedBuffer = std::vector<uint64_t, AlignedAllocator<uint64_t, ALIGN_TO>>;


#endif //INTERFACES_H
