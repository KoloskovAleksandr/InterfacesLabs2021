#pragma once
#include "..\include\IVector.h"
#include "..\include\ILogger.h"
#include "..\include\RC.h"

#define SendWarning(Logger, Code) Logger->warning((Code), __FILE__, __func__, __LINE__)

class Vector : public IVector {
private:
	size_t dim_;    
	Vector(size_t dim);
	double* getMutableData(size_t i = 0);
public:
    static ILogger* logger_;
	static Vector* VectorFactory(size_t dim, double const* const& ptr_data);
    static RC setLogger(ILogger* const logger);
    ~Vector() = default;

    double const* getData() const override;
    RC getCord(size_t index, double& val) const override;
    size_t getDim() const override;

    size_t sizeAllocated() const override;
    IVector* clone() const override;
    RC setData(size_t dim, double const* const& ptr_data) override;
    RC setCord(size_t index, double val) override;


    RC scale(double multiplier) override;
    RC inc(IVector const* const& op) override;
    RC dec(IVector const* const& op) override;
    double norm(NORM n) const override;
   
    RC applyFunction(const std::function<double(double)>& fun) override;
    RC foreach(const std::function<void(double)>& fun) const override;
};