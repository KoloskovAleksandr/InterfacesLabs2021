#pragma once
#include "Interfacedllexport.h"

class LIB_EXPORT IBroker {
public:
    enum class INTERFACE_IMPL {
        IPROBLEM,
        IDIFFPROBLEM,
        ISOLVER,
        TOTAL_INTERFACE_IMPL
    };

    IBroker() = default;

    virtual bool canCastTo(INTERFACE_IMPL impl) const = 0;
    virtual void* getInterfaceImpl(INTERFACE_IMPL impl) const = 0;
    virtual void release() = 0;

    virtual ~IBroker() = 0;
private:
    IBroker(const IBroker&) = delete;
    IBroker& operator=(const IBroker&) = delete;
};

extern "C" {
    LIB_EXPORT void* getBroker();
}
