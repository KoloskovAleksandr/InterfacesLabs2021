#include "KoloskovBroker.h"
#include "KoloskovProblem.h"
#include "KoloskovSolver.h"

Broker* Broker::_instance = nullptr;

Broker* Broker::getInstance() {
    if (!_instance) {
        _instance = new (std::nothrow) Broker();
    }
    return _instance;
}

bool Broker::canCastTo(INTERFACE_IMPL impl) const {
    return impl == INTERFACE_IMPL::IPROBLEM || impl == INTERFACE_IMPL::IDIFFPROBLEM ||
           impl == INTERFACE_IMPL::ISOLVER;
}

void* Broker::getInterfaceImpl(INTERFACE_IMPL impl) const {
    switch (impl){
    case INTERFACE_IMPL::IPROBLEM:
    case INTERFACE_IMPL::IDIFFPROBLEM:
        return (void*) new (std::nothrow) Problem();
    case INTERFACE_IMPL::ISOLVER:
        return (void*) new (std::nothrow) Solver();
    default:
        return nullptr;
    }
}

void Broker::release() {
    delete _instance;
    _instance = nullptr;
}

Broker::~Broker() = default;

IBroker::~IBroker() = default;

void* getBroker() {
    return (void*) Broker::getInstance();
}
