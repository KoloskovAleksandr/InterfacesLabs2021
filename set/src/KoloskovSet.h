#pragma once
#include <memory>
#include "KoloskovControlBlock.h"
#include "..\include\ISet.h"

#define SendWarning(Logger, Code) Logger->warning((Code), __FILE__, __func__, __LINE__)

class Set : public ISet {
private:
    const size_t reallocSize_ = 10;
    size_t dim_;
    size_t memAllocated_;
    size_t containedQuantity_;
    double* data_;
    size_t* map_;
    size_t lastKey_;
    std::shared_ptr<ControlBlock> controlBlock_;


    RC findIndex(IVector const* const& pat, IVector::NORM n, double tol, size_t& index) const;

public:
    enum class Course {
        FORWARD,
        BACK,
    };
    static ILogger* logger_;
    Set();
    static RC setLogger(ILogger* const logger);
    ~Set();

    size_t getDim() const override;
    size_t getSize() const override;

    RC insert(IVector const* const& val, IVector::NORM n, double tol) override;
    RC remove(size_t index) override;
    RC remove(IVector const* const& pat, IVector::NORM n, double tol) override;

    ISet* clone() const override;
    RC getCopy(size_t index, IVector*& val) const override;
    RC findFirstAndCopy(IVector const* const& pat, IVector::NORM n, double tol, IVector*& val) const override;
    RC getCoords(size_t index, IVector* const& val) const override;
    RC findFirstAndCopyCoords(IVector const* const& pat, IVector::NORM n, double tol, IVector* const& val) const override;


    class Iterator : public IIterator {
    private:
        static ILogger* logger_;
        std::shared_ptr<ControlBlock> controlBlock_;
        IVector* vector_;
        size_t index_;
    public:
        Iterator(const std::shared_ptr<ControlBlock>& controlBlock, size_t index, IVector* vector);
        static RC setLogger(ILogger* logger);
        ~Iterator();

        IIterator* getNext(size_t indexInc = 1) const override;
        IIterator* getPrevious(size_t indexInc = 1) const override;
        IIterator* clone() const override;

        RC next(size_t indexInc = 1) override;
        RC previous(size_t indexInc = 1) override;

        bool isValid() const override;

        RC makeBegin() override;
        RC makeEnd() override;

        RC getVectorCopy(IVector*& val) const override;
        RC getVectorCoords(IVector* const& val) const override;

    protected:
        size_t getIndex() const override;

    };

private:
    RC findKey(size_t const key, size_t& index, Course course, size_t steps) const;

public:
    RC MoveBySteps(IVector* const& vec, size_t& key, Course course, size_t steps) const;
    RC MoveBegin(IVector* const& vec, size_t& key) const;
    RC MoveEnd(IVector* const& vec, size_t& key) const;

    IIterator* getIterator(size_t index) const override;
    IIterator* getBegin() const override;
    IIterator* getEnd() const override;

};
