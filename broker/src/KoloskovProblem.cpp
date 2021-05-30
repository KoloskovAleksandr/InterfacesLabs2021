#include "KoloskovProblem.h"
#include <cmath>

ILogger* Problem::logger_ = nullptr;

RC Problem::setLogger(ILogger* const logger) {
    logger_ = logger;
    return RC::SUCCESS;
}

RC IProblem::setLogger(ILogger* const logger) {
    Problem::setLogger(logger);
    return RC::SUCCESS;
}

RC IDiffProblem::setLogger(ILogger* const logger) {
    Problem::setLogger(logger);
    return RC::SUCCESS;
}

Problem::Problem() {
    paramsArea_ = nullptr;
    argsArea_ = nullptr;
    tmpParams_ = nullptr;
    tmpArgs_ = nullptr;
}

RC Problem::setArgsDomain(ICompact const* const& args, ILogger* logger)
{
    if (args->getDim() != argsDim_) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }

    IVector* argsLeftBound = nullptr;
    RC flag = args->getLeftBoundary(argsLeftBound);
    ICompact* argsCopy = args->clone();
    if (flag != RC::SUCCESS) {
        delete argsLeftBound;
        delete argsCopy;
        if (logger_) SendWarning(logger_, flag);
        return flag;
    }
    if (argsArea_ != nullptr || tmpArgs_ != nullptr) {
        delete argsArea_;
        delete tmpArgs_;
    }
    argsArea_ = argsCopy;
    tmpArgs_ = argsLeftBound;
    logger_ = logger;
    return flag;
}

RC Problem::setParamsDomain(ICompact const* const& params)
{
    if (params->getDim() != paramsDim_) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }

    IVector* paramsLeftBound = nullptr;
    RC flag = params->getLeftBoundary(paramsLeftBound);
    if (flag != RC::SUCCESS) {
        delete paramsLeftBound;
        if (logger_) SendWarning(logger_, flag);
        return flag;
    }

    double a, b;
    paramsLeftBound->getCord(PRM::A, a);
    paramsLeftBound->getCord(PRM::B, b);
    if (a <= 0 || b <= 0) {
        delete paramsLeftBound;
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }

    ICompact* paramsCopy = params->clone();
    if (paramsArea_ != nullptr || tmpParams_ != nullptr) {
        delete paramsArea_;
        delete tmpParams_;
    }
    paramsArea_ = paramsCopy;
    tmpParams_ = paramsLeftBound;
    return RC::SUCCESS;
}

IDiffProblem* IDiffProblem::createDiffProblem() {
    return (IDiffProblem*) new(std::nothrow) Problem();
}

IProblem* IProblem::createProblem() {
    return (IProblem*) new(std::nothrow) Problem();
}

IDiffProblem* Problem::clone() const
{
    Problem* copy = new(std::nothrow) Problem();
    if (copy == nullptr) {
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    if (copy->setParamsDomain(paramsArea_) != RC::SUCCESS ||
        copy->setParams(tmpParams_) != RC::SUCCESS ||
        copy->setArgsDomain(argsArea_, nullptr) != RC::SUCCESS ||
        copy->setArgs(tmpArgs_) != RC::SUCCESS) {
        delete copy;
        if (logger_) SendWarning(logger_, RC::ALLOCATION_ERROR);
        return nullptr;
    }
    return (IDiffProblem*)copy;
}

bool Problem::isValidParams(IVector const* const& params) const
{
    return paramsArea_->isInside(params);
}

bool Problem::isValidArgs(IVector const* const& args) const
{
    return argsArea_->isInside(args);
}

RC Problem::setParams(IVector const* const& params)
{
    if (isValidParams(params) == false) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }
    return tmpParams_->setData(params->getDim(), params->getData());
}

RC Problem::setArgs(IVector const* const& args)
{
    if (isValidArgs(args) == false) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }
    return tmpArgs_->setData(args->getDim(), args->getData());
}

double Problem::eval(double const* params, double const* args) const
{
    return params[PRM::A] * args[ARG::X] * args[ARG::X] + params[PRM::B] * args[ARG::Y] * args[ARG::Y]
        + params[PRM::C] * args[ARG::X] + params[PRM::D] * args[ARG::Y]
        + params[PRM::E] * sin(params[PRM::F] * args[ARG::X] + params[PRM::G] * args[ARG::Y]) + params[PRM::H];
}

double Problem::evalByParams(IVector const* const& params) const
{
    if (params->getDim() != paramsDim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return NAN;
    }
    if (isValidParams(params) == false) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return NAN;
    }
    const double* paramsData = params->getData();
    const double* argsData = tmpArgs_->getData();
    return eval(paramsData, argsData);
}

double Problem::evalByArgs(IVector const* const& args) const
{
    if (args->getDim() != argsDim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return NAN;
    }
    if (isValidArgs(args) == false) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return NAN;
    }
    const double* paramsData = tmpParams_->getData();
    const double* argsData = args->getData();
    return eval(paramsData, argsData);
}

double Problem::evalDerivativeByArgs(IVector const* const& args, IMultiIndex const* const& index) const
{
    if (index->getDim() != argsDim_ || args->getDim() != argsDim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return NAN;
    }
    if (isValidArgs(args) == false) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return NAN;
    }
    const size_t* derOrder = index->getData();
    size_t dx = derOrder[ARG::X];
    size_t dy = derOrder[ARG::Y];
    const double* argsData = args->getData();
    const double* paramsData = tmpParams_->getData();

    if (dx == 0 && dy == 0) return evalByArgs(args);
    if (dx == 0) {
        if (dy == 1)
            return 2 * paramsData[PRM::B] * argsData[ARG::Y] + paramsData[PRM::D]
                + paramsData[PRM::E] * paramsData[PRM::G]
                    * std::cos(paramsData[PRM::F] * argsData[ARG::X] + paramsData[PRM::G] * argsData[ARG::Y]);
        if(dy == 2)
            return 2 * paramsData[PRM::B]
                - paramsData[PRM::E] * paramsData[PRM::G] * paramsData[PRM::G]
                    * sin(paramsData[PRM::F] * argsData[ARG::X] + paramsData[PRM::G] * argsData[ARG::Y]);
    }
    if (dy == 0) {
        if (dx == 1)
            return 2 * paramsData[PRM::A] * argsData[ARG::X] + paramsData[PRM::C]
                + paramsData[PRM::E] * paramsData[PRM::F]
                    * cos(paramsData[PRM::F] * argsData[ARG::X] + paramsData[PRM::G] * argsData[ARG::Y]);
        if (dx == 2)
            return 2 * paramsData[PRM::A]
                - paramsData[PRM::E] * paramsData[PRM::F] * paramsData[PRM::F]
                    * sin(paramsData[PRM::F] * argsData[ARG::X] + paramsData[PRM::G] * argsData[ARG::Y]);
    }

    double result = paramsData[PRM::E] * pow(paramsData[PRM::F], dx) * pow(paramsData[PRM::G], dy);
    if (((dx + dy) / 2) % 2 == 1)
        result *= -1.0f;
    if ((dx + dy) % 2 == 0)
        result *= sin(paramsData[PRM::F] * argsData[ARG::X] + paramsData[PRM::G] * argsData[ARG::Y]);
    else
        result *= cos(paramsData[PRM::F] * argsData[ARG::X] + paramsData[PRM::G] * argsData[ARG::Y]);

    return result;
}

double Problem::evalDerivativeByParams(IVector const* const& params, IMultiIndex const* const& index) const
{
    if (index->getDim() != paramsDim_ || params->getDim() != paramsDim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return NAN;
    }
    if (isValidParams(params) == false) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return NAN;
    }
    const size_t* derOrder = index->getData();
    const double* argsData = tmpArgs_->getData();
    const double* paramsData = params->getData();

    size_t sum = 0;
    bool isSolitaryInDer = false;
    size_t solitaryIndex = 0;
    for (size_t i = 0; i < paramsDim_; i++) {
        if (derOrder[i] > 0) {
            sum += derOrder[i];
            if (i < PRM::E || i == PRM::H) {
                isSolitaryInDer = true;
                solitaryIndex = i;
            }
        }
    }

    if (sum == 0) return eval(paramsData, argsData);
    if (isSolitaryInDer) {
        if (sum > 1)
            return 0.0f;
        else {
            switch (solitaryIndex) {
            case PRM::A:
                return argsData[ARG::X] * argsData[ARG::X];
            case PRM::B:
                return argsData[ARG::Y] * argsData[ARG::Y];
            case PRM::C:
                return argsData[ARG::X];
            case PRM::D:
                return argsData[ARG::Y];
            case PRM::H:
                return 1.0f;
            }
        }
    }
    if (derOrder[PRM::E] > 1) return 0.0f;

    double result = pow(argsData[ARG::X], derOrder[PRM::F]) * pow(argsData[ARG::Y], derOrder[PRM::G]);
    if (derOrder[PRM::E] == 0)
        result *= paramsData[PRM::E];
    if (((derOrder[PRM::F] + derOrder[PRM::G]) / 2) % 2 == 1)
        result *= -1.0f;
    if ((derOrder[PRM::F] + derOrder[PRM::G]) % 2 == 0)
        result *= sin(paramsData[PRM::F] * argsData[ARG::X] + paramsData[PRM::G] * argsData[ARG::Y]);
    else
        result *= cos(paramsData[PRM::F] * argsData[ARG::X] + paramsData[PRM::G] * argsData[ARG::Y]);

    return result;
}

RC Problem::evalGradientByArgs(IVector const* const& args, IVector* const& val) const
{
    if (val->getDim() != argsDim_ || args->getDim() != argsDim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return RC::MISMATCHING_DIMENSIONS;
    }
    if (isValidArgs(args) == false) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }
    const double* argsData = args->getData();
    const double* paramsData = tmpParams_->getData();
    val->setCord(ARG::X, 2 * paramsData[PRM::A] * argsData[ARG::X] + paramsData[PRM::C]
        + paramsData[PRM::E] * paramsData[PRM::F]
        * cos(paramsData[PRM::F] * argsData[ARG::X] + paramsData[PRM::G] * argsData[ARG::Y]));
    val->setCord(ARG::Y, 2 * paramsData[PRM::B] * argsData[ARG::Y] + paramsData[PRM::D]
        + paramsData[PRM::E] * paramsData[PRM::G]
        * cos(paramsData[PRM::F] * argsData[ARG::X] + paramsData[PRM::G] * argsData[ARG::Y]));
    return RC::SUCCESS;
}

RC Problem::evalGradientByParams(IVector const* const& params, IVector* const& val) const
{
    if (val->getDim() != paramsDim_ || params->getDim() != paramsDim_) {
        if (logger_) SendWarning(logger_, RC::MISMATCHING_DIMENSIONS);
        return RC::MISMATCHING_DIMENSIONS;
    }
    if (isValidParams(params) == false) {
        if (logger_) SendWarning(logger_, RC::INVALID_ARGUMENT);
        return RC::INVALID_ARGUMENT;
    }
    const double* argsData = tmpArgs_->getData();
    const double* paramsData = params->getData();

    val->setCord(PRM::A, argsData[ARG::X] * argsData[ARG::X]);
    val->setCord(PRM::B, argsData[ARG::Y] * argsData[ARG::Y]);
    val->setCord(PRM::C, argsData[ARG::X]);
    val->setCord(PRM::D, argsData[ARG::Y]);
    val->setCord(PRM::E, sin(paramsData[PRM::F] * argsData[ARG::X] + paramsData[PRM::G] * argsData[ARG::Y]));
    val->setCord(PRM::F, paramsData[PRM::E] * argsData[ARG::X]
                    * cos(paramsData[PRM::F] * argsData[ARG::X] + paramsData[PRM::G] * argsData[ARG::Y]));
    val->setCord(PRM::G, paramsData[PRM::E] * argsData[ARG::Y]
                    * cos(paramsData[PRM::F] * argsData[ARG::X] + paramsData[PRM::G] * argsData[ARG::Y]));
    val->setCord(PRM::H, 1.0f);
}

Problem::~Problem() {
    delete paramsArea_;
    delete argsArea_;
    delete tmpParams_;
    delete tmpArgs_;
}
IDiffProblem::~IDiffProblem() {}
IProblem::~IProblem(){}
