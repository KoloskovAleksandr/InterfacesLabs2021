#include "KoloskovSet.h"
#include <limits>

ILogger* Set::Iterator::logger_ = nullptr;

Set::Iterator::Iterator(const std::shared_ptr<ControlBlock>& controlBlock, size_t index, IVector* vector)
{
    controlBlock_ = controlBlock;
    index_ = index;
    vector_ = vector;
}

RC Set::Iterator::setLogger(ILogger* const logger) {
    logger_ = logger;
    return RC::SUCCESS;
}

RC ISet::IIterator::setLogger(ILogger* const logger) {
    Set::Iterator::setLogger(logger);
    return RC::SUCCESS;
}

Set::Iterator::~Iterator() {
    delete vector_;
}
ISet::IIterator::~IIterator(){}

RC Set::Iterator::makeBegin() {
    return controlBlock_->getBegin(vector_, index_);
}

RC Set::Iterator::makeEnd() {
    return controlBlock_->getEnd(vector_, index_);
}

ISet::IIterator* Set::Iterator::getNext(size_t indexInc) const {
    IIterator* copy = clone();
    if (copy == nullptr) {
        if (Iterator::logger_) SendWarning(Iterator::logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    copy->next(indexInc);
    return copy;
}

ISet::IIterator* Set::Iterator::getPrevious(size_t indexInc) const {
    IIterator* copy = clone();
    if (copy == nullptr) {
        if (Iterator::logger_) SendWarning(Iterator::logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    copy->previous(indexInc);
    return copy;
}

ISet::IIterator* Set::Iterator::clone() const {
    IVector* vector = vector_->clone();
    if (vector == nullptr) {
        if (Iterator::logger_) SendWarning(Iterator::logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    Set::Iterator* iterator = new(std::nothrow) Set::Iterator(controlBlock_, index_, vector);
    if (iterator == nullptr) {
        if (Iterator::logger_) SendWarning(Iterator::logger_, RC::ALLOCATION_ERROR);
        delete vector;
        return nullptr;
    }
    return (ISet::IIterator*)iterator;
}

RC Set::Iterator::next(size_t indexInc) {
    return controlBlock_->getNext(vector_, index_, indexInc);
}

RC Set::Iterator::previous(size_t indexInc) {
    return controlBlock_->getPrevious(vector_, index_, indexInc);
}

bool Set::Iterator::isValid() const
{
    return index_ != std::numeric_limits<size_t>::max();
}


RC Set::Iterator::getVectorCopy(IVector*& val) const {
    IVector* copy = vector_->clone();
    if (copy == nullptr) {
        if (Iterator::logger_) SendWarning(Iterator::logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    val = copy;
    return RC::SUCCESS;
}

RC Set::Iterator::getVectorCoords(IVector* const& val) const {
    return val->setData(vector_->getDim(), vector_->getData());
}

size_t Set::Iterator::getIndex() const {
    return index_;
}

ISet::IIterator* Set::getIterator(size_t index) const {
    if (containedQuantity_ == 0 || index >= containedQuantity_) {
        if (Set::logger_) SendWarning(Set::logger_, RC::INDEX_OUT_OF_BOUND);
        return  nullptr;
    }

    IVector* vector = nullptr;
    RC flag = getCopy(index, vector);
    if (flag != RC::SUCCESS) {
        if (Set::logger_) SendWarning(Set::logger_, flag);
        return  nullptr;
    }

    Set::Iterator* iterator = new(std::nothrow) Set::Iterator(controlBlock_, map_[index], vector);
    if (iterator == nullptr) {
        if (Set::logger_) SendWarning(Set::logger_, RC::ALLOCATION_ERROR);
        delete vector;
        return nullptr;
    }
    Set::Iterator::setLogger(logger_);
    return (ISet::IIterator*)iterator;
}

ISet::IIterator* Set::getBegin() const {
    return getIterator(0);
}

ISet::IIterator* Set::getEnd() const {
    return getIterator(containedQuantity_ - 1);
}

RC Set::findKey(size_t const key, size_t& index, Course course, size_t steps) const {
    size_t left = 0;
    size_t right = containedQuantity_ - 1;
    while (right - left > 1) {
        if (map_[left + (right - left) / 2] < key)
            left = left + (right - left) / 2;
        else
            right = left + (right - left) / 2;
    }

    size_t start = 0;
    if (map_[left] == key || map_[right] == key) {
        if (map_[right] == key)  start = right;
        if (map_[left] == key)  start = left;
    }
    else {
        if (course == Course::BACK)  start = right;
        else  start = left;
    }

    if(course == Course::FORWARD)   index = start + steps;
    else   index = start - steps;
    return RC::SUCCESS;
}

RC Set::MoveBySteps(IVector* const& vec, size_t& key, Course course, size_t steps) const {
    if (key == std::numeric_limits<size_t>::max()) {
        if (Set::logger_) SendWarning(Set::logger_, RC::INDEX_OUT_OF_BOUND);
        return RC::INDEX_OUT_OF_BOUND;
    }
    size_t index;
    findKey(key, index, course, steps);
    if (index >= containedQuantity_) {
        key = std::numeric_limits<size_t>::max();
        return RC::INDEX_OUT_OF_BOUND;
    }
    key = map_[index];
    return getCoords(index, vec);
}

RC Set::MoveBegin(IVector* const& vec, size_t& key) const {
    if (containedQuantity_ == 0) {
        key = std::numeric_limits<size_t>::max();
        if (Set::logger_) SendWarning(Set::logger_, RC::INDEX_OUT_OF_BOUND);
        return RC::INDEX_OUT_OF_BOUND;
    }
    key = map_[0];
    return getCoords(0, vec);
}

RC Set::MoveEnd(IVector* const& vec, size_t& key) const {
    if (containedQuantity_ == 0) {
        key = std::numeric_limits<size_t>::max();
        if (Set::logger_) SendWarning(Set::logger_, RC::INDEX_OUT_OF_BOUND);
        return RC::INDEX_OUT_OF_BOUND;
    }
    key = map_[containedQuantity_ - 1];
    return getCoords(containedQuantity_ - 1, vec);
}
