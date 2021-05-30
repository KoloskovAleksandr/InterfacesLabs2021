#include "KoloskovMultiIndex.h"
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cstddef>
#include <functional>

ILogger* MultiIndex::logger_ = nullptr;

size_t* MultiIndex::getMutableData(size_t i) {
    return (size_t*)((uint8_t*)this + sizeof(MultiIndex)) + i;
}

MultiIndex::MultiIndex(size_t dim) {
    dim_ = dim;
}

MultiIndex* MultiIndex::MultiIndexFactory(size_t dim, const size_t* indices)
{
    void* mem = malloc(sizeof(MultiIndex) + dim * sizeof(size_t));
    if (mem == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }

    MultiIndex* multiIndex = new(mem) MultiIndex(dim);

    memcpy(multiIndex->getMutableData(), indices, dim * sizeof(size_t));
    return multiIndex;
}

IMultiIndex* IMultiIndex::createMultiIndex(size_t dim, const size_t* indices)
{
    return (IMultiIndex*)MultiIndex::MultiIndexFactory(dim, indices);
}

RC MultiIndex::setLogger(ILogger* const logger) {
    logger_ = logger;
    return RC::SUCCESS;
}

RC IMultiIndex::setLogger(ILogger* const logger) {
    MultiIndex::setLogger(logger);
    return RC::SUCCESS;
}

IMultiIndex* MultiIndex::clone() const
{
    return IMultiIndex::createMultiIndex(getDim(), getData());
}

size_t MultiIndex::getDim() const
{
    return dim_;
}

const size_t* MultiIndex::getData() const
{
    return (const size_t*)((uint8_t*)this + sizeof(MultiIndex));
}

RC MultiIndex::setData(size_t dim, size_t const* const& ptr_data)
{
    if(ptr_data == nullptr){
        if (logger_) SendWarning(logger_, RC::NULLPTR_ERROR);
        return RC::NULLPTR_ERROR;
    }
    if (getDim() != dim) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return RC::MISMATCHING_DIMENSIONS;
    }

    memcpy(getMutableData(), ptr_data, dim * sizeof(size_t));
    return RC::SUCCESS;
}

RC MultiIndex::getAxisIndex(size_t axisIndex, size_t& val) const
{
    if (axisIndex >= getDim()) {
        if (logger_) SendWarning(logger_, RC::INDEX_OUT_OF_BOUND);
        return RC::INDEX_OUT_OF_BOUND;
    }
    val = *(getData() + axisIndex);
    return RC::SUCCESS;
}

RC MultiIndex::setAxisIndex(size_t axisIndex, size_t val)
{
    if (dim_ <= axisIndex) {
        if (logger_) SendWarning(logger_, RC::INDEX_OUT_OF_BOUND);
        return RC::INDEX_OUT_OF_BOUND;
    }

    *(getMutableData(axisIndex)) = val;
    return RC::SUCCESS;
}

RC MultiIndex::incAxisIndex(size_t axisIndex, size_t val)
{
    if (axisIndex >= dim_) {
        if (logger_) SendWarning(logger_, RC::INDEX_OUT_OF_BOUND);
        return RC::INDEX_OUT_OF_BOUND;
    }

    getMutableData()[axisIndex] += val;
    return RC::SUCCESS;
}

MultiIndex::~MultiIndex(){}
IMultiIndex::~IMultiIndex(){}
