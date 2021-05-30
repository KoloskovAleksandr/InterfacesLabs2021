#pragma once
#include "../include/ISetControlBlock.h"

class Set;

class ControlBlock : public ISetControlBlock {
public:
    ControlBlock(Set* set);
    RC getNext(IVector* const& vec, size_t& index, size_t indexInc = 1) const override;
    RC getPrevious(IVector* const& vec, size_t& index, size_t indexInc = 1) const override;
    RC getBegin(IVector* const& vec, size_t& index) const override;
    RC getEnd(IVector* const& vec, size_t& index) const override;

private:
    Set* set_;
};