#include "KoloskovVector.h"
#include <cstring>
#include <cmath>
using namespace std;

ILogger* Vector::logger_ = nullptr;

Vector::Vector(size_t dim) { dim_ = dim; }

double* Vector::getMutableData(size_t i) {
    return (double*)((uint8_t*)this + sizeof(Vector)) + i;
}

double const* Vector::getData() const {
    return (double const*)((uint8_t*)this + sizeof(Vector));
}

RC Vector::getCord(size_t index, double& val) const  {
    if (index >= getDim()) {
        if (logger_) SendWarning(logger_, RC::INDEX_OUT_OF_BOUND);
        return RC::INDEX_OUT_OF_BOUND;
    }
    val = *(getData() + index);
    return RC::SUCCESS;
}

size_t Vector::getDim() const { return dim_; }

RC Vector::setData(size_t dim, double const* const& ptr_data) {
    ///IAA: ptr_data == nullptr?
    if(ptr_data == nullptr){
        if (logger_) SendWarning(logger_, RC::NULLPTR_ERROR);
        return RC::NULLPTR_ERROR;
    }
    if (getDim() != dim) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return RC::MISMATCHING_DIMENSIONS;
    }
    for (size_t i = 0; i < dim; i++) {
        if (isnan(ptr_data[i]) || isinf(ptr_data[i])) {
            if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
            return RC::INVALID_ARGUMENT;
        }
    }
    memcpy(getMutableData(), ptr_data, dim * sizeof(double));
    return RC::SUCCESS;
}

RC Vector::setCord(size_t index, double val) {
    if (dim_ <= index) {
        if (logger_) SendWarning(logger_, RC::INDEX_OUT_OF_BOUND);
        return RC::INDEX_OUT_OF_BOUND;
    }
    if (isnan(val) || isinf(val)) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }
    *(getMutableData(index)) = val;
    return RC::SUCCESS;
}

Vector* Vector::VectorFactory(size_t dim, double const* const& ptr_data) {
    void* mem = malloc(sizeof(Vector) + dim * sizeof(double));
    if (mem == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }

    Vector* vector = new(mem) Vector(dim);
    RC isSetData = vector->setData(dim, ptr_data);
    if (isSetData != RC::SUCCESS) {
        if (logger_) SendWarning(logger_, isSetData);
        free(mem);
        return nullptr;
    }
    return vector;
}

IVector* IVector::createVector(size_t dim, double const* const& ptr_data) {
    return (IVector*)Vector::VectorFactory(dim, ptr_data);
}

RC Vector::setLogger(ILogger* const logger) {
    logger_ = logger;
    return RC::SUCCESS;
}

RC IVector::setLogger(ILogger* const logger) {
    Vector::setLogger(logger);
    return RC::SUCCESS;
}

IVector::~IVector(){}

size_t Vector::sizeAllocated() const {
    return sizeof(Vector) + getDim() * sizeof(double);
}

RC IVector::copyInstance(IVector* const dest, IVector const* const& src) {
    size_t destSize = dest->sizeAllocated();
    size_t srcSize = src->sizeAllocated();
    if (destSize != srcSize) {
        if (Vector::logger_) SendWarning(Vector::logger_, RC::MISMATCHING_DIMENSIONS);
        return RC::MISMATCHING_DIMENSIONS;
    }
    if ((dest <= src && (uint8_t*)dest >= (uint8_t*)src - destSize) ||
        (src <= dest && (uint8_t*)src >= (uint8_t*)dest - destSize)) {
        if (Vector::logger_) SendWarning(Vector::logger_, RC::MEMORY_INTERSECTION);
        return RC::MEMORY_INTERSECTION;
    }
    memcpy(dest, src, srcSize);
    return RC::SUCCESS;
}

RC IVector::moveInstance(IVector* const dest, IVector*& src) {
    RC isCopy = copyInstance(dest, src);
    if (isCopy == RC::SUCCESS) {
        if (Vector::logger_) SendWarning(Vector::logger_, isCopy);
        delete src;
        src = nullptr;
    }
    return isCopy;
}

IVector* Vector::clone() const {
    return IVector::createVector(getDim(), getData());
}

RC Vector::scale(double multiplier){
    if (isnan(multiplier) || isinf(multiplier)) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }
    double* data = getMutableData();
    for (size_t i = 0; i < dim_; i++) {
        if (isinf(data[i] * multiplier)) {
            if (logger_) SendWarning(logger_, RC::INFINITY_OVERFLOW);
            ///IAA: and invalid vector!
            if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
            return RC::INFINITY_OVERFLOW;
        }
        data[i] *= multiplier;
    }
    return RC::SUCCESS;
}

RC Vector::inc(IVector const* const& op){
    if (op->getDim() != dim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return RC::MISMATCHING_DIMENSIONS;
    }
    double* data = getMutableData();
    double const* opData = op->getData();
    for (size_t i = 0; i < dim_; i++) {
        data[i] += opData[i];
        if (isnan(data[i]) || isinf(data[i])) {
            if (logger_) SendWarning(logger_, RC::INFINITY_OVERFLOW);
            return RC::INFINITY_OVERFLOW;
        }
    }
    return RC::SUCCESS;
}

RC Vector::dec(IVector const* const& op) {
    if (op->getDim() != dim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return RC::MISMATCHING_DIMENSIONS;
    }
    double* data = getMutableData();
    double const* opData = op->getData();
    for (size_t i = 0; i < dim_; i++) {
        data[i] -= opData[i];
        if (isnan(data[i]) || isinf(data[i])) {
            if (logger_) SendWarning(logger_, RC::INFINITY_OVERFLOW);
            return RC::INFINITY_OVERFLOW;
        }
    }
    return RC::SUCCESS;
}

IVector* IVector::add(IVector const* const& op1, IVector const* const& op2){
    IVector* result = createVector(op1->getDim(), op1->getData());
    if (result == nullptr) {
        if (Vector::logger_) SendWarning(Vector::logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    RC isIncremented = result->inc(op2);
    if (isIncremented != RC::SUCCESS) {
        if (Vector::logger_) SendWarning(Vector::logger_, isIncremented);
        return nullptr;
    }
    return result;
}

IVector* IVector::sub(IVector const* const& op1, IVector const* const& op2){
    IVector* result = createVector(op1->getDim(), op1->getData());
    if (result == nullptr) {
        if (Vector::logger_) SendWarning(Vector::logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    RC isDecremented = result->dec(op2);
    if (isDecremented != RC::SUCCESS) {
        if (Vector::logger_) SendWarning(Vector::logger_, isDecremented);
        return nullptr;
    }
    return result;
}

double Vector::norm(NORM n) const {
    double norm = 0;
    double const* data = getData();
    switch (n) {
    case NORM::CHEBYSHEV:
        for (size_t i = 0; i < dim_; i++)
            norm = fmax(norm, fabs(data[i]));
        break;
    case NORM::FIRST:
        for (size_t i = 0; i < dim_; i++)
            norm += fabs(data[i]);
        break;
    case NORM::SECOND:
        for (size_t i = 0; i < dim_; i++)
            norm += pow(data[i], 2.0f);
        norm = sqrt(norm);
        break;
    default:
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return (double)NAN;
    }

    if (isinf(norm)) {
        if (logger_) SendWarning(logger_, RC::INFINITY_OVERFLOW);
        return (double)NAN;
    }
    return norm;
}

double IVector::dot(IVector const* const& op1, IVector const* const& op2) {
    if (op1->getDim() != op2->getDim()) {
        if (Vector::logger_) SendWarning(Vector::logger_, RC::MISMATCHING_DIMENSIONS);
        return (double)NAN;
    }
    double const* data1 = op1->getData();
    double const* data2 = op2->getData();
    double dot = 0;
    for (size_t i = 0; i < op1->getDim(); i++)
        dot += data1[i] * data2[i];

    if (isinf(dot)) {
        if (Vector::logger_) SendWarning(Vector::logger_, RC::INFINITY_OVERFLOW);
        return (double)NAN;
    }
    return dot;
}

bool IVector::equals(IVector const* const& op1, IVector const* const& op2, NORM n, double tol) {
    IVector* delta = IVector::sub(op1, op2);
    if (delta == nullptr) {
        if (Vector::logger_) SendWarning(Vector::logger_, RC::INVALID_ARGUMENT);
        return false;
    }
    double norm = delta->norm(n);
    delete delta;

    if (isnan(norm)) {
        if (Vector::logger_) SendWarning(Vector::logger_, RC::NOT_NUMBER);
        return false;
    }
    return norm <= tol;
}

RC Vector::applyFunction(const std::function<double(double)>& fun) {
    double* data = getMutableData();
    double result;
    for (size_t i = 0; i < dim_; i++) {
        result = fun(data[i]);
        if (isnan(result) || isinf(result)) {
            if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
            return RC::INVALID_ARGUMENT;
        }
        data[i] = result;
    }
    return RC::SUCCESS;
}

RC Vector::foreach(const std::function<void(double)>& fun) const {
    double const* data = getData();
    for (size_t i = 0; i < dim_; i++)
        fun(data[i]);
    return RC::SUCCESS;
}
