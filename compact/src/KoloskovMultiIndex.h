#pragma once
#include "..\include\IMultiIndex.h"

#define SendWarning(Logger, Code) Logger->warning((Code), __FILE__, __func__, __LINE__)

class MultiIndex : public IMultiIndex {
private:
    size_t dim_;
    MultiIndex(size_t dim);
    size_t* getMutableData(size_t i = 0);
public:
    static ILogger* logger_;
    static RC setLogger(ILogger* const pLogger);

    static MultiIndex* MultiIndexFactory(size_t dim, const size_t* indices);

    IMultiIndex* clone() const override;
    size_t getDim() const override;
    const size_t* getData() const override;
    RC setData(size_t dim, size_t const* const& ptr_data) override;

    RC getAxisIndex(size_t axisIndex, size_t& val) const override;
    RC setAxisIndex(size_t axisIndex, size_t val) override;
    RC incAxisIndex(size_t axisIndex, size_t val) override;

    ~MultiIndex();   
};