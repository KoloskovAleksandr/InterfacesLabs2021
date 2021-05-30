#pragma once
#include "../include/IVector.h"
#include "../include/IMultiIndex.h"
#include "Interfacedllexport.h"
#include <cstddef>

class LIB_EXPORT ICompactControlBlock {
public:
    /*
    * Control block moves iterator forward in grid of compact. When reached last node of an axis, next axis position will be incremented. Next axis defined by bypassOrder
    *
    * @param [in] currentIndex Multi-index of current iterator position
    *
    * @param [in] bypassOrder Multi-index, that defining bypass orred of axis
    */
    virtual RC get(IMultiIndex* const& currentIndex, IMultiIndex const* const& bypassOrder) const = 0;
    /*
    * Control block calculates vector corresponding to multi-index
    *
    * @param [in] currentIndex Multi-index of current iterator position
    *
    * @param [in] val Buffer for vector data
    */
    virtual RC get(IMultiIndex const* const& currentIndex, IVector* const& val) const = 0;

    virtual ~ICompactControlBlock() = 0;

private:
    ICompactControlBlock(ICompactControlBlock const&);
    ICompactControlBlock& operator=(ICompactControlBlock const&);

protected:
    ICompactControlBlock() = default;
};
