#pragma once
#include "../include/IDiffProblem.h"

#define SendWarning(Logger, Code) Logger->warning((Code), __FILE__, __func__, __LINE__)

class Problem : public IDiffProblem {
private:
    static const size_t argsDim_ = 2;
    static const size_t paramsDim_ = 8;
    ICompact* paramsArea_;
    ICompact* argsArea_;
    IVector* tmpParams_;
    IVector* tmpArgs_;

    double eval(double const* params, double const* args) const;
public:
    enum ARG {
        X = 0,
        Y = 1,
    };
    enum PRM {
        A = 0,
        B = 1,
        C = 2,
        D = 3,
        E = 4,
        F = 5,
        G = 6,
        H = 7,
    };
    static ILogger* logger_;
    static RC setLogger(ILogger* const logger);

    //Ax^2 + By^2 + Cx + Dy + Esin(Fx + Gy) + H, A > 0, B > 0
    Problem();
    RC setArgsDomain(ICompact const* const& args, ILogger* logger) override;
    RC setParamsDomain(ICompact const* const& params) override;
    ~Problem();
    IDiffProblem* clone() const override;

    bool isValidParams(IVector const* const& params) const override;
    bool isValidArgs(IVector const* const& args) const override;

    RC setParams(IVector const* const& params) override;
    RC setArgs(IVector const* const& args) override;

    double evalByParams(IVector const* const& params) const override;
    double evalByArgs(IVector const* const& args) const override;

    double evalDerivativeByArgs(IVector const* const& args, IMultiIndex const* const& index) const override;
    double evalDerivativeByParams(IVector const* const& params, IMultiIndex const* const& index) const override;

    RC evalGradientByArgs(IVector const* const& args, IVector* const& val) const override;
    RC evalGradientByParams(IVector const* const& params, IVector* const& val) const override;

};
