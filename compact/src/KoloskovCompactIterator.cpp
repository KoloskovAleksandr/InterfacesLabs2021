#include "KoloskovCompact.h"

ILogger* Compact::Iterator::logger_ = nullptr;

RC Compact::Iterator::setLogger(ILogger* const pLogger) {
    logger_ = pLogger;
    return RC::SUCCESS;
}

RC ICompact::IIterator::setLogger(ILogger* const pLogger) {
    return Compact::Iterator::setLogger(pLogger);
}

CompactControlBlock::CompactControlBlock(Compact* compact) {
    compact_ = compact;
}

RC CompactControlBlock::get(IMultiIndex* const& currentIndex, IMultiIndex const* const& bypassOrder) const
{
    return compact_->CalculateShift(currentIndex, bypassOrder);
}

RC CompactControlBlock::get(IMultiIndex const* const& currentIndex, IVector* const& val) const
{
    return compact_->getVectorCoords(currentIndex, val);
}

Compact::Iterator::Iterator(std::shared_ptr<CompactControlBlock> const& cb, IMultiIndex* startPos, IVector* vector, IMultiIndex* order, bool isAtEnd)
{
    isAtEnd_ = isAtEnd;
    vector_ = vector;
    coords_ = startPos;
    order_ = order;
    controlBlock_ = cb;
}

ICompact::IIterator* Compact::Iterator::getMovedIterator()
{
    ICompact::IIterator* copy = clone();
    if (copy == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    copy->moveForwardIterator();
    return copy;
}

bool Compact::Iterator::isAtEnd() const {
    return isAtEnd_;
}

ICompact::IIterator* Compact::Iterator::clone() const {
    IVector* vector = vector_->clone();
    IMultiIndex* coords = coords_->clone();
    IMultiIndex* order = order_->clone();
    if (vector == nullptr || coords == nullptr || order == nullptr) {
        delete vector;
        delete coords;
        delete order;
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    return new(std::nothrow) Compact::Iterator(controlBlock_, coords, vector, order, isAtEnd_);
}

RC Compact::Iterator::moveForwardIterator() {
    size_t zeroCoord;
    size_t increment;

    RC flag = order_->getAxisIndex(0, zeroCoord);
    if (flag != RC::SUCCESS) {
        if (logger_) SendWarning(logger_, flag);
        return flag;
    }
    flag = coords_->getAxisIndex(zeroCoord, increment);
    if (flag != RC::SUCCESS) {
        if (logger_) SendWarning(logger_, flag);
        return flag;
    }

    coords_->setAxisIndex(zeroCoord, increment + 1);
    flag = controlBlock_->get(coords_, order_);
    if (flag != RC::SUCCESS && flag != RC::INDEX_OUT_OF_BOUND) {
        if (logger_) SendWarning(logger_, flag);
        return flag;
    }
    if (flag == RC::INDEX_OUT_OF_BOUND) {
        isAtEnd_ = true;
        return RC::INDEX_OUT_OF_BOUND;
    }

    return controlBlock_->get(coords_, vector_);
}

RC Compact::Iterator::getVectorCopy(IVector*& val) const {
    val = vector_->clone();
    if (val == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    return RC::SUCCESS;
}

RC Compact::Iterator::getVectorCoords(IVector* const& val) const {
    if (val->getDim() != vector_->getDim()) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return RC::MISMATCHING_DIMENSIONS;
    }
    return val->setData(vector_->getDim(), vector_->getData());
}

ICompact::IIterator* Compact::getIterator(IMultiIndex const* const& index, IMultiIndex const* const& bypassOrder) const
{
    if (index->getDim() != dim_ || bypassOrder->getDim() != dim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return nullptr;
    }
    if (isOrderValid(bypassOrder) != RC::SUCCESS)
        return nullptr;

    IMultiIndex* coords = index->clone();
    IMultiIndex* order = bypassOrder->clone();
    IVector* vector = nullptr;
    getVectorCopy(index, vector);
    if (coords == nullptr || order == nullptr || vector == nullptr) {
        delete coords;
        delete order;
        delete vector;
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }

    RC flag = CalculateShift(coords, order);
    if(flag == RC::SUCCESS)
        return (ICompact::IIterator*)new(std::nothrow) Compact::Iterator(controlBlock_, coords, vector, order, false);
    if(flag == RC::INDEX_OUT_OF_BOUND)
        return (ICompact::IIterator*)new(std::nothrow) Compact::Iterator(controlBlock_, coords, vector, order, true);
    if (logger_) SendWarning(logger_, flag);
    return nullptr;
}

ICompact::IIterator* Compact::getBegin(IMultiIndex const* const& bypassOrder) const
{
    if (bypassOrder->getDim() != dim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return nullptr;
    }
    if (isOrderValid(bypassOrder) != RC::SUCCESS)
        return nullptr;

    size_t* index = (size_t*)calloc(dim_, sizeof(size_t));
    if (index == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    IMultiIndex* order = bypassOrder->clone();
    IVector* vector = nullptr;
    getLeftBoundary(vector);
    IMultiIndex* coords = IMultiIndex::createMultiIndex(dim_, index);
    if (coords == nullptr || order == nullptr || vector == nullptr) {
///IAA: delete index?
        free(index);
        delete coords;
        delete order;
        delete vector;
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    free(index);
    return (ICompact::IIterator*)new(std::nothrow) Compact::Iterator(controlBlock_, coords, vector, order, false);
}

ICompact::IIterator* Compact::getEnd(IMultiIndex const* const& bypassOrder) const
{
    if (bypassOrder->getDim() != dim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return nullptr;
    }
    if (isOrderValid(bypassOrder) != RC::SUCCESS)
        return nullptr;

    size_t* index = (size_t*)calloc(dim_, sizeof(size_t));
    if (index == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    const size_t* gridData = grid_->getData();
    for (size_t i = 0; i < dim_; i++) {
        index[i] = gridData[i] - 1;
    }

    IMultiIndex* order = bypassOrder->clone();
    IVector* vector = nullptr;
    getRightBoundary(vector);
    IMultiIndex* coords = IMultiIndex::createMultiIndex(dim_, index);
    if (coords == nullptr || order == nullptr || vector == nullptr) {
        delete coords;
        delete order;
        delete vector;
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    free(index);
    return (ICompact::IIterator*)new(std::nothrow) Compact::Iterator(controlBlock_, coords, vector, order, true);
}

RC Compact::CalculateShift(IMultiIndex* const& currentIndex, IMultiIndex const* const& bypassOrder) const
{
    if (isIndexValid(currentIndex)) return RC::SUCCESS;
    if (dim_ == 1) {
        size_t gridData;
        grid_->getAxisIndex(0, gridData);
        currentIndex->setAxisIndex(0, gridData - 1);
        return RC::SUCCESS;
    }

    const size_t* order = bypassOrder->getData();
    const size_t* curCoords = currentIndex->getData();
    const size_t* gridData = grid_->getData();
    bool isIndexValid = false;
    size_t coordIndex = 0;
    size_t i, j;
    while (isIndexValid == false) {
        i = order[coordIndex];
        if (curCoords[i] >= gridData[i]) {
            j = order[coordIndex + 1];
            currentIndex->setAxisIndex(i, 0);
            currentIndex->incAxisIndex(j, 1);
            if (coordIndex + 1 == dim_ - 1 && curCoords[j] >= gridData[j]) {
                for (size_t k = 0; k < dim_; k++)
                    currentIndex->setAxisIndex(k, gridData[k] - 1);
                return RC::INDEX_OUT_OF_BOUND;
            }
        }
        else {
            isIndexValid = true;
        }
        coordIndex++;
    }
    return RC::SUCCESS;
}

bool Compact::isIndexValid(IMultiIndex const* const& index) const
{
    size_t const* indexData = index->getData();
    size_t const* gridData = grid_->getData();
    for (size_t i = 0; i < dim_; i++)
        if (indexData[i] >= gridData[i])
            return false;
    return true;
}

RC Compact::isOrderValid(IMultiIndex const* const& bypassOrder) const
{
    const size_t* order = bypassOrder->getData();
    for (size_t i = 0; i < dim_; i++) {
        bool isValid = false;
        for (size_t j = 0; j < dim_; j++) {
            if (order[j] == i)
                isValid = true;
        }
        if (isValid == false) {
            if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
            return RC::INVALID_ARGUMENT;
        }
    }
    return RC::SUCCESS;
}

Compact::Iterator::~Iterator() {
    delete coords_;
    delete vector_;
    delete order_;
}
ICompact::IIterator::~IIterator() {}
CompactControlBlock::~CompactControlBlock() {}
ICompactControlBlock::~ICompactControlBlock() {}
