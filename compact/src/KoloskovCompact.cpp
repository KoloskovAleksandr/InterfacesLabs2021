#include "KoloskovCompact.h"
constexpr size_t defaultGridSize_ = 10;
ILogger* Compact::logger_ = nullptr;

RC Compact::setLogger(ILogger* const logger) {
    logger_ = logger;
    return RC::SUCCESS;
}

RC ICompact::setLogger(ILogger* const logger) {
    Compact::setLogger(logger);
    return RC::SUCCESS;
}

Compact::Compact(size_t dim) : controlBlock_(new(std::nothrow) CompactControlBlock(this)) {
    dim_ = dim;
///IAA: тогда уж leftmost
    leftMost_ = nullptr;
    rightMost_ = nullptr;
    grid_ = nullptr;
}

RC Compact::findLeftestAndRightest(IVector const* const& left, IVector const* const& right)
{
    if (left == nullptr || right == nullptr) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }

    double* leftestCoords = (double*)malloc(dim_ * sizeof(double));
    double* rightestCoords = (double*)malloc(dim_ * sizeof(double));
    if (leftestCoords == NULL || rightestCoords == NULL) {
        if (leftestCoords) free(leftestCoords);
        if (rightestCoords) free(rightestCoords);
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }

    const double* leftCoords = left->getData();
    const double* rightCoords = right->getData();
    for (size_t i = 0; i < dim_; i++) {
        leftestCoords[i] = std::min(leftCoords[i], rightCoords[i]);
        rightestCoords[i] = std::max(leftCoords[i], rightCoords[i]);
    }

    leftMost_ = IVector::createVector(dim_, leftestCoords);
    rightMost_ = IVector::createVector(dim_, rightestCoords);
    if (leftMost_ == nullptr || rightMost_ == nullptr) {
        delete leftMost_;
        delete rightMost_;
        free(leftestCoords);
        free(rightestCoords);
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    free(leftestCoords);
    free(rightestCoords);
    return RC::SUCCESS;
}

RC Compact::setGrid(IMultiIndex const* nodeQuantities)
{
    if (nodeQuantities->getDim() != dim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return RC::MISMATCHING_DIMENSIONS;
    }
    IMultiIndex* grid = nodeQuantities->clone();
    if (grid == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
///IAA: current grid_ == nullptr?
    if (grid_ != nullptr)
        delete grid_;
    grid_ = grid;
    return RC::SUCCESS;
}

Compact* Compact::createCompact(IVector const* vec1, IVector const* vec2, IMultiIndex const* nodeQuantities) {
    if (!vec1 || !vec2 || !nodeQuantities) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return nullptr;
    }
    if (vec1->getDim() != vec2->getDim() || vec1->getDim() != nodeQuantities->getDim()) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return nullptr;
    }

    Compact* compact = new(std::nothrow) Compact(vec1->getDim());
    if (compact == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }

    RC flag = compact->findLeftestAndRightest(vec1, vec2);
    if (flag != RC::SUCCESS) {
        if (logger_) SendWarning(logger_, flag);
        delete compact;
        return nullptr;
    }

    flag = compact->setGrid(nodeQuantities);
    if (flag != RC::SUCCESS) {
        if (logger_) SendWarning(logger_, flag);
        delete compact;
        return nullptr;
    }
    return compact;
}

Compact::~Compact() {
    delete leftMost_;
    delete rightMost_;
    delete grid_;
}
ICompact::~ICompact(){}

ICompact* ICompact::createCompact(IVector const* vec1, IVector const* vec2, IMultiIndex const* nodeQuantities) {
    return (ICompact*)Compact::createCompact(vec1, vec2, nodeQuantities);
}

ICompact* Compact::clone() const
{
	return (ICompact*)Compact::createCompact(leftMost_, rightMost_, grid_);
}

size_t Compact::getDim() const
{
    return dim_;
}

IMultiIndex* Compact::getGrid() const
{
    IMultiIndex* copy = grid_->clone();
    if (copy == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    return copy;
}

bool Compact::isInside(IVector const* const& vec) const
{
    if (vec->getDim() != dim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return false;
    }
    const double* leftCoords = leftMost_->getData();
    const double* rightCoords = rightMost_->getData();
    const double* vecCoords = vec->getData();
    for (size_t i = 0; i < dim_; i++) {
        if (vecCoords[i] < leftCoords[i] || vecCoords[i] > rightCoords[i])
            return false;
    }
    return true;
}

RC Compact::getVectorCopy(IMultiIndex const* index, IVector*& val) const
{
    IVector* result = leftMost_->clone();
    if (result == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    RC flag = getVectorCoords(index, result);
    if (flag != RC::SUCCESS) {
        delete result;
        if (logger_) SendWarning(logger_, flag);
        return flag;
    }
    val = result;
    return RC::SUCCESS;
}

RC Compact::getVectorCoords(IMultiIndex const* index, IVector* const& val) const
{
    if (index->getDim() != dim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return RC::MISMATCHING_DIMENSIONS;
    }
    const size_t* indexCoords = index->getData();
    const size_t* gridLength = grid_->getData();
    for (size_t i = 0; i < dim_; i++) {
        if (indexCoords[i] > gridLength[i]) {
            if (logger_) SendWarning(logger_, RC::INDEX_OUT_OF_BOUND);
            return RC::INDEX_OUT_OF_BOUND;
        }
    }

    const double* leftCoords = leftMost_->getData();
    const double* rightCoords = rightMost_->getData();
    for (size_t i = 0; i < dim_; i++) {
        double lambda = (double)(indexCoords[i]) / (gridLength[i] - 1);
        val->setCord(i, (1.0 - lambda) * leftCoords[i] + lambda * rightCoords[i]);
    }
    return RC::SUCCESS;
}

RC Compact::getLeftBoundary(IVector*& val) const
{
    IVector* result = leftMost_->clone();
    if (result == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    val = result;
    return RC::SUCCESS;
}

RC Compact::getRightBoundary(IVector*& val) const
{
    IVector* result = rightMost_->clone();
    if (result == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    val = result;
    return RC::SUCCESS;
}

//for simplicity, features of my implementation are used
ICompact* ICompact::createIntersection(ICompact const* op1, ICompact const* op2, IMultiIndex const* const grid, double tol) {
    if (op1->getDim() != op2->getDim()) return nullptr;
    IVector* op1Left = nullptr;
    IVector* op1Right = nullptr;
    IVector* op2Left = nullptr;
    IVector* op2Right = nullptr;
    op1->getLeftBoundary(op1Left);
    op2->getLeftBoundary(op2Left);
    op1->getRightBoundary(op1Right);
    op2->getRightBoundary(op2Right);
    if (op1Left == nullptr || op1Right == nullptr || op2Left == nullptr || op2Right == nullptr) {
        delete op1Left;
        delete op2Left;
        delete op1Right;
        delete op2Right;
        return nullptr;
    }
    IVector* intersectionLeft = op1Left->clone();
    IVector* intersectionRight = op1Left->clone();
    IMultiIndex* intersectionGrid = op1->getGrid();
    if (intersectionLeft == nullptr || intersectionRight == nullptr || intersectionGrid == nullptr) {
        delete op1Left;
        delete op2Left;
        delete op1Right;
        delete op2Right;
        delete intersectionLeft;
        delete intersectionRight;
        delete intersectionGrid;
        return nullptr;
    }

    const double* op1LeftCoords = op1Left->getData();
    const double* op1RightCoords = op1Right->getData();
    const double* op2LeftCoords = op2Left->getData();
    const double* op2RightCoords = op2Right->getData();
    bool isIntersect = true;
    for (size_t i = 0; i < op1->getDim(); i++) {
        if (op1LeftCoords[i] >= op2LeftCoords[i] && op1LeftCoords[i] <= op2RightCoords[i]) {
            intersectionLeft->setCord(i, op1LeftCoords[i]);
            if (op2RightCoords[i] - op1LeftCoords[i] < tol) {
                intersectionRight->setCord(i, op1LeftCoords[i]);
                intersectionGrid->setAxisIndex(i, 0);
            }
            else {
                intersectionRight->setCord(i, op2RightCoords[i]);
                if (grid != nullptr) {
                    size_t buf;
                    grid->getAxisIndex(i, buf);
                    intersectionGrid->setAxisIndex(i, buf);
                }
                else
                    intersectionGrid->setAxisIndex(i, defaultGridSize_);
            }
        }
        else if (op2LeftCoords[i] >= op1LeftCoords[i] && op2LeftCoords[i] <= op1RightCoords[i]) {
            intersectionLeft->setCord(i, op2LeftCoords[i]);
            if (op1RightCoords[i] - op2LeftCoords[i] < tol) {
                intersectionRight->setCord(i, op2LeftCoords[i]);
                intersectionGrid->setAxisIndex(i, 0);
            }
            else {
                intersectionRight->setCord(i, op1RightCoords[i]);
                if (grid != nullptr) {
                    size_t buf;
                    grid->getAxisIndex(i, buf);
                    intersectionGrid->setAxisIndex(i, buf);
                }
                else
                    intersectionGrid->setAxisIndex(i, defaultGridSize_);
            }

        }
        else {
            isIntersect = false;
        }
    }

    delete op1Left;
    delete op2Left;
    delete op1Right;
    delete op2Right;

    ICompact* Intersection = nullptr;

    if (isIntersect)
        Intersection = ICompact::createCompact(intersectionLeft, intersectionRight, intersectionGrid);
    delete intersectionLeft;
    delete intersectionRight;
    delete intersectionGrid;

    return Intersection;
}

ICompact* ICompact::createCompactSpan(ICompact const* op1, ICompact const* op2, IMultiIndex const* const grid) {
    if (op1->getDim() != op2->getDim()) return nullptr;
    IVector* op1Left = nullptr;
    IVector* op1Right = nullptr;
    IVector* op2Left = nullptr;
    IVector* op2Right = nullptr;
    op1->getLeftBoundary(op1Left);
    op2->getLeftBoundary(op2Left);
    op1->getRightBoundary(op1Right);
    op2->getRightBoundary(op2Right);
    if (op1Left == nullptr || op1Right == nullptr || op2Left == nullptr || op2Right == nullptr) {
        delete op1Left;
        delete op2Left;
        delete op1Right;
        delete op2Right;
        return nullptr;
    }
    IVector* spanLeft = op1Left->clone();
    IVector* spanRight = op1Left->clone();
    IMultiIndex* spanGrid = op1->getGrid();
    if (spanLeft == nullptr || spanRight == nullptr || spanGrid == nullptr) {
        delete op1Left;
        delete op2Left;
        delete op1Right;
        delete op2Right;
        delete spanLeft;
        delete spanRight;
        delete spanGrid;
        return nullptr;
    }

    const double* op1LeftCoords = op1Left->getData();
    const double* op1RightCoords = op1Right->getData();
    const double* op2LeftCoords = op2Left->getData();
    const double* op2RightCoords = op2Right->getData();


    for (size_t i = 0; i < op1->getDim(); i++) {
        spanLeft->setCord(i, std::min(op1LeftCoords[i], op2LeftCoords[i]));
        spanRight->setCord(i, std::max(op1RightCoords[i], op2RightCoords[i]));
        spanGrid->setAxisIndex(i, defaultGridSize_);
    }
    delete op1Left;
    delete op2Left;
    delete op1Right;
    delete op2Right;


    ICompact* span = nullptr;
    if (grid == nullptr) {
        span = ICompact::createCompact(spanLeft, spanRight, spanGrid);
    }
    else {
        span = ICompact::createCompact(spanLeft, spanRight, grid);
    }
    delete spanLeft;
    delete spanRight;
    delete spanGrid;
    return span;
}
