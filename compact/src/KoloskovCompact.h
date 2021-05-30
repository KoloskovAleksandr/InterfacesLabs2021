#pragma once
#include "..\include\ICompact.h"
#include "..\include\ICompactControlBlock.h"
#include <memory>
#define SendWarning(Logger, Code) Logger->warning((Code), __FILE__, __func__, __LINE__)

class CompactControlBlock;

class Compact : public ICompact {
private:
    size_t dim_;
    IVector* leftMost_;
    IVector* rightMost_;
    IMultiIndex* grid_;
    std::shared_ptr<CompactControlBlock> controlBlock_;

    Compact(size_t dim);
    RC findLeftestAndRightest(IVector const* const& left, IVector const* const& right);
    RC setGrid(IMultiIndex const* nodeQuantities);
public:
    static ILogger* logger_;
    static RC setLogger(ILogger* const logger);
    static Compact* createCompact(IVector const* vec1, IVector const* vec2, IMultiIndex const* nodeQuantities);
    ~Compact();


    ICompact* clone() const override;

    size_t getDim() const override;
    IMultiIndex* getGrid() const override;

    bool isInside(IVector const* const& vec) const override;

    RC getVectorCopy(IMultiIndex const* index, IVector*& val) const override;
    RC getVectorCoords(IMultiIndex const* index, IVector* const& val) const override;

    RC getLeftBoundary(IVector*& val) const override;
    RC getRightBoundary(IVector*& val) const override;

    class Iterator : public Compact::IIterator {
    private:
        bool isAtEnd_;
        IVector* vector_;
        IMultiIndex* coords_;
        IMultiIndex* order_;
        std::shared_ptr<CompactControlBlock> controlBlock_;
    public:
        static ILogger* logger_;
        static RC setLogger(ILogger* const pLogger);
        Iterator(std::shared_ptr<CompactControlBlock> const& block, IMultiIndex* startPos, IVector* vector, IMultiIndex* order, bool isAtEnd);
        IIterator* getMovedIterator() override;
        IIterator* clone() const override;

        RC moveForwardIterator() override;
        bool isAtEnd() const override;

        RC getVectorCopy(IVector*& val) const override;
        RC getVectorCoords(IVector* const& val) const override;

        ~Iterator();
    };

    IIterator* getIterator(IMultiIndex const* const& index, IMultiIndex const* const& bypassOrder) const override;
    IIterator* getBegin(IMultiIndex const* const& bypassOrder) const override;
    IIterator* getEnd(IMultiIndex const* const& bypassOrder) const override;
    RC CalculateShift(IMultiIndex* const& currentIndex, IMultiIndex const* const& bypassOrder) const;
private:
    RC isOrderValid(IMultiIndex const* const& bypassOrder) const;
    bool isIndexValid(IMultiIndex const* const& index) const;
};

class CompactControlBlock : public ICompactControlBlock {
private:
    Compact* compact_;
public:
    CompactControlBlock(Compact* compact);
    RC get(IMultiIndex* const& currentIndex, IMultiIndex const* const& bypassOrder) const override;
    RC get(IMultiIndex const* const& currentIndex, IVector* const& val) const override;
    ~CompactControlBlock();
};
