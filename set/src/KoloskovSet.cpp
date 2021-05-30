#include "KoloskovSet.h"
#include <cstring>

ILogger* Set::logger_ = nullptr;

RC Set::setLogger(ILogger* const logger) {
    logger_ = logger;
    return RC::SUCCESS;
}

RC ISet::setLogger(ILogger* const logger) {
    Set::setLogger(logger);
    return RC::SUCCESS;
}

Set::Set():controlBlock_(new ControlBlock(this)) {
    dim_ = 0;
    memAllocated_ = 0;
    containedQuantity_ = 0;
    lastKey_ = 0;
    data_ = nullptr;
    map_ = nullptr;
}

///IAA: Зачем этот метод?
//Set* Set::SetFactory() {


ISet* ISet::createSet() {
    return (ISet*)new(std::nothrow) Set();
}

Set::~Set() {
    free(data_);
    free(map_);
}

ISet::~ISet(){}

size_t Set::getDim() const {
    return dim_;
}

size_t Set::getSize() const {
    return containedQuantity_;
}

RC Set::findIndex(IVector const* const& pat, IVector::NORM n, double tol, size_t& index) const
{
    if (dim_ != pat->getDim()) {
        if (Set::logger_) SendWarning(Set::logger_, RC::MISMATCHING_DIMENSIONS);
        return RC::MISMATCHING_DIMENSIONS;
    }

    IVector* inSet = IVector::createVector(dim_, data_);
    if (inSet == nullptr) {
        if (Set::logger_) SendWarning(Set::logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    for (size_t i = 0; i < containedQuantity_; i++) {
        inSet->setData(dim_, data_ + i * dim_);
        if(IVector::equals(inSet, pat, n, tol)){
            index = i;
            delete inSet;
            return RC::SUCCESS;
        }
    }
    delete inSet;
    return RC::VECTOR_NOT_FOUND;
}

RC Set::insert(IVector const* const& val, IVector::NORM n, double tol) {
    if (dim_ == 0) dim_ = val->getDim();
    if (dim_ != val->getDim()) {
        if (Set::logger_) SendWarning(Set::logger_, RC::MISMATCHING_DIMENSIONS);
        return RC::MISMATCHING_DIMENSIONS;
    }
    if (dim_ * containedQuantity_ == memAllocated_) {
        double* buf = (double*)malloc((dim_ * reallocSize_ + memAllocated_) * sizeof(double));
        size_t* fub = (size_t*)malloc((reallocSize_ + memAllocated_ / dim_) * sizeof(size_t));
        if (buf == nullptr || fub == nullptr) {
            if (Set::logger_) SendWarning(Set::logger_, RC::ALLOCATION_ERROR);
            ///IAA: memory leaks
            free(buf);
            free(fub);
            return RC::ALLOCATION_ERROR;
        }
        if (data_ != nullptr) {
            memmove(buf, data_, dim_ * containedQuantity_ * sizeof(double));
            free(data_);
        }
        if (map_ != nullptr) {
            memmove(fub, map_, containedQuantity_ * sizeof(size_t));
            free(map_);
        }
        data_ = buf;
        map_ = fub;
        memAllocated_ += dim_ * reallocSize_;
    }

    size_t index = 0;
    if (findIndex(val, n, tol, index) == RC::VECTOR_NOT_FOUND) {
        memmove(data_ + containedQuantity_ * dim_, val->getData(), dim_ * sizeof(double));
        map_[containedQuantity_] = lastKey_;
        containedQuantity_++;
        lastKey_++;
    }
    return RC::SUCCESS;
}

RC Set::remove(size_t index) {
    if (containedQuantity_ <= index) {
        if (Set::logger_) SendWarning(Set::logger_, RC::INDEX_OUT_OF_BOUND);
        return RC::INDEX_OUT_OF_BOUND;
    }
    memcpy(data_ + dim_ * index, data_ + dim_ * (index + 1), (containedQuantity_ - index - 1) * dim_ * sizeof(double));
    memcpy(map_ + index, map_ + index + 1, (containedQuantity_ - index - 1) * sizeof(size_t));
    containedQuantity_--;
    return RC::SUCCESS;
}

RC Set::remove(IVector const* const& pat, IVector::NORM n, double tol) {
    RC flag = RC::SUCCESS;
    size_t index;
    while (flag == RC::SUCCESS) {
        flag = findIndex(pat, n, tol, index);
        if (flag == RC::SUCCESS)
            remove(index);
    }
    if (flag == RC::VECTOR_NOT_FOUND)
        return RC::SUCCESS;
    return flag;
}

ISet* Set::clone() const {
    Set* copy = new(std::nothrow) Set();
    if (copy == nullptr) {
        if (Set::logger_) SendWarning(Set::logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    copy->setLogger(logger_);
    copy->dim_ = dim_;
    copy->memAllocated_ = memAllocated_;
    copy->containedQuantity_ = containedQuantity_;
    copy->lastKey_ = lastKey_;

    copy->data_ = (double*)malloc(memAllocated_ * sizeof(double));
    copy->map_ = (size_t*)malloc(memAllocated_ / dim_ * sizeof(size_t));
    if (copy->data_ == nullptr || copy->map_ == nullptr) {
        if (Set::logger_) SendWarning(Set::logger_, RC::ALLOCATION_ERROR);
        free(copy->data_);
        free(copy->map_);
        delete copy;
        return nullptr;
    }
    memcpy(copy->data_, data_, memAllocated_ * sizeof(double));
    memcpy(copy->map_, map_, memAllocated_ / dim_ * sizeof(size_t));
    for (size_t i = 0; i < containedQuantity_; i++)
        printf("%i ", copy->map_[i]);
    printf("\n");
    return copy;
}

RC Set::getCopy(size_t index, IVector*& val) const {
    if (containedQuantity_ <= index) {
        if (Set::logger_) SendWarning(Set::logger_, RC::INDEX_OUT_OF_BOUND);
        return RC::INDEX_OUT_OF_BOUND;
    }
    IVector* vector = IVector::createVector(dim_, data_ + dim_ * index);
    if (vector == nullptr) {
        if (Set::logger_) SendWarning(Set::logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    val = vector;
    return RC::SUCCESS;
}

RC Set::findFirstAndCopy(IVector const* const& pat, IVector::NORM n, double tol, IVector*& val) const {
    size_t index;
    RC isFind = findIndex(pat, n, tol, index);
    if (isFind != RC::SUCCESS) return isFind;
    return getCopy(index, val);
}

RC Set::getCoords(size_t index, IVector* const& val) const {
    if (dim_ != val->getDim()) {
        if (Set::logger_) SendWarning(Set::logger_, RC::MISMATCHING_DIMENSIONS);
        return RC::MISMATCHING_DIMENSIONS;
    }
    if (containedQuantity_ <= index) {
        if (Set::logger_) SendWarning(Set::logger_, RC::INDEX_OUT_OF_BOUND);
        return RC::INDEX_OUT_OF_BOUND;
    }
    val->setData(dim_, data_ + dim_ * index);
    return RC::SUCCESS;
}

RC Set::findFirstAndCopyCoords(IVector const* const& pat, IVector::NORM n, double tol, IVector* const& val) const {
    size_t index;
    RC isFind = findIndex(pat, n, tol, index);
    if (isFind != RC::SUCCESS) return isFind;
    return getCoords(index, val);
}

ISet* ISet::makeIntersection(ISet const* const& op1, ISet const* const& op2, IVector::NORM n, double tol) {
    if (op1->getDim() != op2->getDim())  return nullptr;
    ISet* Intersection = createSet();
    if (Intersection == nullptr) return nullptr;
    if (op1->getSize() == 0 || op2->getSize() == 0)
        return Intersection;

    RC flag = RC::SUCCESS;
    IVector* iterOp1 = nullptr;
    IVector* iterOp2 = nullptr;
    flag = op1->getCopy(0, iterOp1);
    flag = op2->getCopy(0, iterOp2);
    if (flag != RC::SUCCESS) {
        delete iterOp1;
        delete iterOp2;
        delete Intersection;
        return nullptr;
    }

    for (size_t i = 0; i < op1->getSize() && flag == RC::SUCCESS; i++) {
        op1->getCoords(i, iterOp1);
        for (size_t j = 0; j < op2->getSize() && flag == RC::SUCCESS; j++) {
            op2->getCoords(j, iterOp2);
            if (IVector::equals(iterOp1, iterOp2, n, tol))
                flag = Intersection->insert(iterOp2, n, tol);
        }
    }
    delete iterOp1;
    delete iterOp2;
    if (flag != RC::SUCCESS) {
        delete Intersection;
        return nullptr;
    }
    return Intersection;
}

ISet* ISet::makeUnion(ISet const* const& op1, ISet const* const& op2, IVector::NORM n, double tol) {
    if (op1->getDim() != op2->getDim()) return nullptr;
    ISet* Union = createSet();
    if (Union == nullptr) return nullptr;

    if (op1->getSize() == 0) return op2->clone();

    RC flag = RC::SUCCESS;
    IVector* iterator = nullptr;
    flag = op1->getCopy(0, iterator);
    if (flag != RC::SUCCESS) {
        delete Union;
        return nullptr;
    }

    for (size_t i = 0; i < op1->getSize() && flag == RC::SUCCESS; i++) {
        op1->getCoords(i, iterator);
        flag = Union->insert(iterator, n, tol);
    }
    for (size_t i = 0; i < op2->getSize() && flag == RC::SUCCESS; i++) {
        op2->getCoords(i, iterator);
        flag = Union->insert(iterator, n, tol);
    }
    delete iterator;

    if (flag != RC::SUCCESS) {
        delete Union;
        return nullptr;
    }
    return Union;
}

ISet* ISet::sub(ISet const* const& op1, ISet const* const& op2, IVector::NORM n, double tol) {
    if (op1->getDim() != op2->getDim()) return nullptr;
    ISet* Sub = createSet();
    if (Sub == nullptr) return nullptr;
    if (op1->getSize() == 0) {
        return Sub;
    }

    RC flag = RC::SUCCESS;
    IVector* iterator = nullptr;
    flag = op1->getCopy(0, iterator);
    if (flag != RC::SUCCESS) {
        delete Sub;
        return nullptr;
    }

    for (size_t i = 0; i < op1->getSize() && flag == RC::SUCCESS; i++) {
        op1->getCoords(i, iterator);
        flag = Sub->insert(iterator, n, tol);
    }
    for (size_t i = 0; i < op2->getSize() && flag == RC::SUCCESS; i++) {
        op2->getCoords(i, iterator);
        flag = Sub->remove(iterator, n, tol);
    }
    delete iterator;

    if (flag != RC::SUCCESS) {
        delete Sub;
        return nullptr;
    }
    return Sub;
}

ISet* ISet::symSub(ISet const* const& op1, ISet const* const& op2, IVector::NORM n, double tol) {
    ISet* Union = Set::makeUnion(op1, op2, n, tol);
    ISet* Intersection = Set::makeIntersection(op1, op2, n, tol);
    if (Union == nullptr || Union->getSize() == 0 ||
        Intersection == nullptr || Intersection->getSize() == 0){
        ///IAA: memory leaks
        delete Intersection;
        return Union;
    }

    IVector* iterator = nullptr;
    Intersection->getCopy(0, iterator);
    if (iterator == nullptr) {
        delete Union;
        delete Intersection;
        return nullptr;
    }

    RC flag = RC::SUCCESS;
    for (size_t i = 0; i < Intersection->getSize() && flag == RC::SUCCESS; i++) {
        Intersection->getCoords(i, iterator);
        Union->remove(iterator, n, tol);
    }
    delete iterator;
    delete Intersection;

    if (flag != RC::SUCCESS) {
        delete Union;
        return nullptr;
    }
    return Union;
}

bool ISet::subSet(ISet const* const& op1, ISet const* const& op2, IVector::NORM n, double tol) {
    if (op1->getDim() != op2->getDim()) return false;

    RC flag = RC::SUCCESS;
    IVector* iterOp1 = nullptr;
    IVector* iterOp2 = nullptr;
    flag = op1->getCopy(0, iterOp1);
    flag = op2->getCopy(0, iterOp2);
    if (flag != RC::SUCCESS) {
        delete iterOp1;
        delete iterOp2;
        return false;
    }

    bool isSubset = true;
    for (size_t i = 0; i < op1->getSize() && isSubset == true; i++) {
        isSubset = false;
        op1->getCoords(i, iterOp1);
        for (size_t j = 0; j < op2->getSize() && isSubset == false; j++) {
            op2->getCoords(j, iterOp2);
            isSubset = IVector::equals(iterOp1, iterOp2, n, tol);
        }
    }
    delete iterOp1;
    delete iterOp2;
    return isSubset;
}

bool ISet::equals(ISet const* const& op1, ISet const* const& op2, IVector::NORM n, double tol) {
    return subSet(op1, op2, n, tol) && subSet(op1, op2, n, tol);
}
