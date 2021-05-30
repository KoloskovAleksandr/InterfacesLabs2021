#include "KoloskovControlBlock.h"
#include "KoloskovSet.h"

ControlBlock::ControlBlock(Set* set) {	set_ = set;  }
RC ControlBlock::getNext(IVector* const& vec, size_t& index, size_t indexInc) const
{
	return set_->MoveBySteps(vec, index, Set::Course::FORWARD, indexInc);
}
RC ControlBlock::getPrevious(IVector* const& vec, size_t& index, size_t indexInc) const
{
	return set_->MoveBySteps(vec, index, Set::Course::BACK, indexInc);	
}
RC ControlBlock::getBegin(IVector* const& vec, size_t& index) const
{
	return set_->MoveBegin(vec, index);	
}
RC ControlBlock::getEnd(IVector* const& vec, size_t& index) const
{
	return set_->MoveEnd(vec, index);
}
ISetControlBlock::~ISetControlBlock() {}