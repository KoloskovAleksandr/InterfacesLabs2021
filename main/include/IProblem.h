#pragma once
#include "ICompact.h"
#include "ILogger.h"
#include "IVector.h"
#include "Interfacedllexport.h"


class IProblem {
public:
    static IProblem* createProblem();
    virtual RC setArgsDomain(ICompact const* const& args, ILogger* logger = nullptr) = 0;
    virtual RC setParamsDomain(ICompact const* const& params) = 0;
    virtual IProblem* clone() const = 0;

    static RC setLogger(ILogger* const logger);

    virtual bool isValidParams(IVector const* const& params) const = 0;
    virtual bool isValidArgs(IVector const* const& args) const = 0;

    virtual RC setParams(IVector const* const& params) = 0;
    virtual RC setArgs(IVector const* const& args) = 0;

    virtual double evalByArgs(IVector const* const& args) const = 0;
    virtual double evalByParams(IVector const* const& params) const = 0;

    virtual ~IProblem() = 0;

private:
    IProblem(const IProblem&) = delete;
    IProblem& operator=(const IProblem&) = delete;

protected:
    IProblem() = default;
};
