#include "KoloskovSolver.h"
#include <iostream>

ILogger* Solver::logger_ = nullptr;

RC Solver::setLogger(ILogger* const logger) {
    logger_ = logger;
    return RC::SUCCESS;
}

RC ISolver::setLogger(ILogger* const logger) {
    Solver::setLogger(logger);
    return RC::SUCCESS;
}

Solver::Solver() {
    problem_ = nullptr;
    argsArea_ = nullptr;
    paramsArea_ = nullptr;
    solution_ = nullptr;
}

//Solver* Solver::createSolver() {

///IAA: What's for???
ISolver* ISolver::createSolver() {
    return (ISolver*)new(std::nothrow) Solver();
}

ISolver* Solver::clone() const
{
    ISolver* copy = ISolver::createSolver();
    if (copy == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }

    RC flag = RC::SUCCESS;
    if (problem_) {
        flag = copy->setProblem(problem_);
        if (flag != RC::SUCCESS) {
            delete copy;
            if (logger_) SendWarning(logger_, flag);
            return nullptr;
        }
    }
    if (argsArea_) {
        flag = copy->setArgsDomain(argsArea_, nullptr);
        if (flag != RC::SUCCESS) {
            delete copy;
            if (logger_) SendWarning(logger_, flag);
            return nullptr;
        }
    }
    if (paramsArea_) {
        flag = copy->setParamsDomain(paramsArea_);
        if (flag != RC::SUCCESS) {
            delete copy;
            if (logger_) SendWarning(logger_, flag);
            return nullptr;
        }
    }
    return copy;
}

RC Solver::setProblem(IDiffProblem const* const& pProblem)
{
    IDiffProblem* problemCopy = pProblem->clone();
    if (problemCopy == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    if (problem_)  delete problem_;
    problem_ = problemCopy;
    return RC::SUCCESS;
}

bool Solver::isValidArgsDomain(ICompact const* const& args) const
{
    if (problem_ == nullptr) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return false;
    }

    IVector* argsLeftBound = nullptr;
    RC flag = args->getLeftBoundary(argsLeftBound);
    if (flag != RC::SUCCESS) {
        if (logger_) SendWarning(logger_, flag);
        delete argsLeftBound;
        return false;
    }
    IVector* argsRightBound = nullptr;
    flag = args->getRightBoundary(argsRightBound);
    if (flag != RC::SUCCESS) {
        if (logger_) SendWarning(logger_, flag);
        delete argsRightBound;
        delete argsLeftBound;
        return false;
    }

    if (argsLeftBound->getDim() > maxProblemDim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        delete argsRightBound;
        delete argsLeftBound;
        return false;
    }

    bool isValid = problem_->isValidArgs(argsLeftBound) && problem_->isValidArgs(argsRightBound);
    delete argsRightBound;
    delete argsLeftBound;
    return isValid;
}

bool Solver::isValidParamsDomain(ICompact const* const& params) const
{
    IVector* paramsLeftBound = nullptr;
    RC flag = params->getLeftBoundary(paramsLeftBound);
    if (flag != RC::SUCCESS) {
        delete paramsLeftBound;
        return false;
    }
    IVector* paramsRightBound = nullptr;
    flag = params->getRightBoundary(paramsRightBound);
    if (flag != RC::SUCCESS) {
        delete paramsRightBound;
        delete paramsLeftBound;
        return false;
    }

    double leftCoord, rightCoord;
    for (size_t i = 0; i < paramsDim_; i++) {
        paramsLeftBound->getCord(i, leftCoord);
        paramsRightBound->getCord(i, rightCoord);
        if (leftCoord <= 0 || rightCoord >= 1) {
            delete paramsRightBound;
            delete paramsLeftBound;
            return false;
        }
    }

    delete paramsRightBound;
    delete paramsLeftBound;
    return true;
}

RC Solver::setArgsDomain(ICompact const* const& args, ILogger* logger)
{
    if (isValidArgsDomain(args) == false) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }
    ICompact* argsClone = args->clone();
    if (argsClone == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    if (argsArea_ != nullptr) {
        delete argsArea_;
    }
    argsArea_ = argsClone;
    problemDim_ = argsArea_->getDim();
    logger_ = logger;

    return RC::SUCCESS;
}

RC Solver::setParamsDomain(ICompact const* const& params)
{
    if (isValidParamsDomain(params) == false) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }
    ICompact* paramsClone = params->clone();
    if (paramsClone == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    if (paramsArea_ != nullptr) {
        delete paramsArea_;
    }
    paramsArea_ = paramsClone;
    return RC::SUCCESS;
}

RC Solver::solveByArgs(IVector const* const& initArg, IVector const* const& solverParams)
{
    if (problem_ == nullptr || argsArea_ == nullptr || paramsArea_ == nullptr) {
        if (logger_) SendWarning(logger_, RC::NULLPTR_ERROR);
        return RC::NULLPTR_ERROR;
    }
    if (argsArea_->isInside(initArg) == false || paramsArea_->isInside(solverParams) == false) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }
    return solve(initArg, solverParams);
}

RC Solver::solveByParams(IVector const* const& initParam, IVector const* const& solverParams)
{
    return solveByArgs(initParam, solverParams);
}

RC Solver::getSolution(IVector*& solution) const
{
    if(solution != nullptr){
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }
    if (solution_ == nullptr) {
        if (logger_) SendWarning(logger_, RC::NULLPTR_ERROR);
        return RC::NULLPTR_ERROR;
    }
    IVector* solutionCopy = solution_->clone();
    if (solutionCopy == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
///IAA: memory leaks!
    solution = solutionCopy;
    return RC::SUCCESS;
}

RC Solver::solve(IVector const* const& initArg, IVector const* const& solverParams) {
    double alphaInit, delta, lambda, eps;
    RC flag = solverParams->getCord(PRM::ALPH, alphaInit);
    if (flag != RC::SUCCESS) {
        if (logger_) SendWarning(logger_, flag);
        return flag;
    }
    flag = solverParams->getCord(PRM::DELT, delta);
    if (flag != RC::SUCCESS) {
        if (logger_) SendWarning(logger_, flag);
        return flag;
    }
    flag = solverParams->getCord(PRM::LAMB, lambda);
    if (flag != RC::SUCCESS) {
        if (logger_) SendWarning(logger_, flag);
        return flag;
    }
    flag = solverParams->getCord(PRM::EPS, eps);
    if (flag != RC::SUCCESS) {
        if (logger_) SendWarning(logger_, flag);
        return flag;
    }

    size_t* derData = (size_t*)calloc(problemDim_, sizeof(size_t));
    if (derData == NULL) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    IVector* grad = initArg->clone();
    IVector* x = initArg->clone();
    IVector* step = initArg->clone();
    IVector* xbuf = initArg->clone();
    IVector* stepbuf = initArg->clone();
    IMultiIndex* derOrder = IMultiIndex::createMultiIndex(problemDim_, derData);
    if (grad == nullptr || x == nullptr || derOrder == nullptr) {
        free(derData);
        delete(grad);
        delete(x);
        delete(derOrder);
        delete(step);
        delete(xbuf);
        delete(stepbuf);
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    double* xData = (double*)malloc(problemDim_ * sizeof(double));
    double* antigrad = (double*)malloc(problemDim_ * sizeof(double));
    double* stepData = (double*)malloc(problemDim_ * sizeof(double));
    double** hess = (double**)malloc(problemDim_ * sizeof(double*));
    if (xData == NULL || antigrad == NULL || stepData == NULL || hess == NULL) {
        free(xData);
        free(antigrad);
        free(stepData);
        free(hess);
        free(derData);
        delete(grad);
        delete(x);
        delete(derOrder);
        delete(step);
        delete(xbuf);
        delete(stepbuf);
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return RC::ALLOCATION_ERROR;
    }
    for (size_t i = 0; i < problemDim_; i++) {
        hess[i] = (double*)malloc(problemDim_ * sizeof(double));
        if (hess[i] == NULL) {
            for (size_t j = 0; j < i; j++)
                free(hess[j]);
            free(xData);
            free(antigrad);
            free(stepData);
            free(hess);
            free(derData);
            delete(grad);
            delete(x);
            delete(derOrder);
            delete(step);
            delete(xbuf);
            delete(stepbuf);
            if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
            return RC::ALLOCATION_ERROR;
        }
    }


    for (size_t i = 0; i < maxIterations_; i++) {
        problem_->evalGradientByArgs(x, grad);
        if (grad->norm(IVector::NORM::CHEBYSHEV) <= eps)
            break;

        double val;
        for (size_t j = 0; j < problemDim_; j++) {
            grad->getCord(j, val);
            antigrad[j] = -val;
        }
///IAA: What's the method supposed to be?
///KAO: Newton's method at each iteration uses the calculation of the target function's
///     Hessian in x_k point. This process is abstracted by this function.
        constructHessian(hess, x, derOrder);
        solveSystem(hess, antigrad, stepData);
        step->setData(problemDim_, stepData);

        double alpha = searchAlpha(alphaInit, lambda, delta * IVector::dot(grad, step), x, step, xbuf, stepbuf, hess, derOrder);
        step->scale(alpha);
        x->inc(step);
    }

    flag = RC::SUCCESS;
    if (grad->norm(IVector::NORM::CHEBYSHEV) > eps) {
        if (logger_) SendWarning(logger_, RC::VECTOR_NOT_FOUND);
        flag = RC::VECTOR_NOT_FOUND;
    }
    if (solution_ == nullptr) {
        solution_ = initArg->clone();
        if (solution_ == nullptr) {
            if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
            flag = RC::ALLOCATION_ERROR;
        }
    }
    solution_->setData(problemDim_, x->getData());

    free(xData);
    free(antigrad);
    free(stepData);
    free(derData);
    delete(grad);
    delete(x);
    delete(derOrder);
    delete(step);
    delete(xbuf);
    delete(stepbuf);
    for (size_t i = 0; i < problemDim_; i++)
        free(hess[i]);
    free(hess);

    return flag;
}

void Solver::constructHessian(double** hess, IVector const* const& x, IMultiIndex* derOrder) {
    for (size_t i = 0; i < problemDim_; i++)
        derOrder->setAxisIndex(i, 0);
    for (size_t i = 0; i < problemDim_; i++) {
        derOrder->incAxisIndex(i, 1);
        for (size_t j = 0; j < problemDim_; j++) {
            derOrder->incAxisIndex(j, 1);
            hess[i][j] = problem_->evalDerivativeByArgs(x, derOrder);
            if (i == j)
                derOrder->setAxisIndex(i, 1);
            else
                derOrder->setAxisIndex(j, 0);
        }
        derOrder->setAxisIndex(i, 0);
    }
    return;
}

double Solver::det(double** A) {
    if (problemDim_ == 1)
        return A[0][0];
    if (problemDim_ == 2)
        return A[0][0] * A[1][1] - A[0][1] * A[1][0];
    if (problemDim_ == 3)
        return A[0][0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) +
        A[1][0] * (A[0][2] * A[2][1] - A[0][1] * A[2][2]) +
        A[2][0] * (A[0][1] * A[1][2] - A[0][2] * A[1][1]);
    return 0.0;
}

bool Solver::isPositiveDef(double** matrix) {
    return det(matrix) > 0;
}

void Solver::solveSystem(double** A, double* b, double* xData) {
    if (problemDim_ == 1) {
        xData[0] = b[0] / A[0][0];
        return;
    }
    double determ = det(A);
    if (problemDim_ == 2) {
        xData[0] = b[0] * A[1][1] - b[1] * A[0][1];
        xData[1] = -b[0] * A[1][0] + b[1] * A[0][0];
    }
    else {
        xData[0] = b[0] * (A[1][1] * A[2][2] - A[1][2] * A[2][1]) +
            b[1] * (A[0][2] * A[2][1] - A[0][1] * A[2][2]) +
            b[2] * (A[0][1] * A[1][2] - A[0][2] * A[1][1]);
        xData[1] = b[0] * (A[1][2] * A[2][0] - A[1][0] * A[2][2]) +
            b[1] * (A[0][0] * A[2][2] - A[0][2] * A[2][0]) +
            b[2] * (A[0][2] * A[1][0] - A[0][0] * A[1][2]);
        xData[2] = b[0] * (A[1][0] * A[2][1] - A[1][1] * A[2][0]) +
            b[1] * (A[0][1] * A[2][0] - A[0][0] * A[2][1]) +
            b[2] * (A[0][0] * A[1][1] - A[0][1] * A[1][0]);
    }
    for (size_t i = 0; i < problemDim_; i++)
        xData[i] /= determ;
    return;
}

double Solver::searchAlpha(double alphaInit, double lambda, double eps, IVector* x, IVector* step, IVector* xbuf, IVector* stepbuf, double** hess, IMultiIndex* derOrder) {
    double alpha = alphaInit;
    double f_x = problem_->evalByArgs(x);
    stepbuf->setData(problemDim_, step->getData());
    stepbuf->scale(alpha);
    xbuf->setData(problemDim_, x->getData());
    xbuf->inc(stepbuf);
    constructHessian(hess, xbuf, derOrder);
    while (problem_->evalByArgs(xbuf) - f_x > alpha * eps && isPositiveDef(hess)) {
        alpha *= lambda;
        stepbuf->setData(problemDim_, step->getData());
        stepbuf->scale(alpha);
        xbuf->setData(problemDim_, x->getData());
        xbuf->inc(stepbuf);
        constructHessian(hess, xbuf, derOrder);
    }
    return alpha;
}


Solver::~Solver()
{
    delete problem_;
    delete paramsArea_;
    delete argsArea_;
    delete solution_;
}
ISolver::~ISolver(){}

