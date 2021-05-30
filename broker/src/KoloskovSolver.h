#pragma once
#include "..\include\ISolver.h"

#define SendWarning(Logger, Code) Logger->warning((Code), __FILE__, __func__, __LINE__)

class Solver : public ISolver {
private:
    IDiffProblem* problem_;
    ICompact* argsArea_;
    ICompact* paramsArea_;
    IVector* solution_;
    size_t argsDim_;
    static const size_t maxIterations_ = 1000;
    static const size_t paramsDim_ = 4;
    static const size_t maxProblemDim_ = 3;
    size_t problemDim_;

    RC solve(IVector const* const& initArg, IVector const* const& solverParams);
    void constructHessian(double** hess, IVector const* const& x, IMultiIndex* derOrder);
    double det(double** matrix);
    bool isPositiveDef(double** matrix);
    void solveSystem(double** matrix, double* b, double* xData);
    double searchAlpha(double alphaInit, double lambda, double eps, IVector* x, IVector* step, IVector* xbuf, IVector* stepbuf, double** hess, IMultiIndex* derOrder);
public:
    enum PRM {
        ALPH = 0,
        DELT = 1,
        LAMB = 2,
        EPS = 3,
    };
    Solver();
    static ILogger* logger_;
    //Newton method with step choice using splitting
    //With params: alpha_0 - Initial antigradient compression coefficient,
    //delta - Normative coefficient for the yield condition in splitting,
    //lambda - splitting coefficient
    //epsilon - yield condition constant
    static Solver* createSolver();
    static RC setLogger(ILogger* const pLogger);

    ISolver* clone() const override;

    RC setProblem(IDiffProblem const* const& pProblem) override;

    bool isValidArgsDomain(ICompact const* const& args) const override;
    bool isValidParamsDomain(ICompact const* const& params) const override;
    RC setArgsDomain(ICompact const* const& args, ILogger* logger) override;
    RC setParamsDomain(ICompact const* const& params) override;

    // initArg - starting point (x0)
    // solverParams - defined by developer of exact Solver, described in his header
    RC solveByArgs(IVector const* const& initArg, IVector const* const& solverParams) override;
    RC solveByParams(IVector const* const& initParam, IVector const* const& solverParams) override;
    RC getSolution(IVector*& solution) const override;

    ~Solver();
};
