#pragma once
#include "../include/IBroker.h"
#include "../include/IProblem.h"

namespace {
    class Broker : public IBroker {
    public:
        static Broker* _instance;
        static Broker* getInstance();

        Broker() = default;

        bool canCastTo(INTERFACE_IMPL impl) const override;
        void* getInterfaceImpl(INTERFACE_IMPL impl) const override;

        virtual void release();

        private:
        ~Broker() override;
    };
};
