#include <iostream>
#include <limits>
#include <windows.h>
#include "include/ILogger.h"
#include "include/IVector.h"
#include "include/ISet.h"
#include "include/ICompact.h"
#include "include/IDiffProblem.h"
#include "include/ISolver.h"
#include "include/IBroker.h"
using namespace std;


typedef void*(*BrokerProc)();
function<void(double)> print = [](double num) { std::cout << num << ' '; };

void testingVector(ILogger* logger){
    double a_data[] = {1, 2, 3, 4};
    double b_data[] = {-1, -2, -3, -4};
    double c_data[] = {0, 0, 0, 0};
    IVector* a = IVector::createVector(4, a_data);
    IVector* b = IVector::createVector(4, b_data);
    IVector* c = IVector::createVector(4, c_data);
    IVector::setLogger(logger);
    if(!a || !b || !c){
        cout << "Allocation error" << endl;
    }
    cout << "Vector a is:" << endl;
    a->foreach(print);
    cout << endl;
    cout << "Vector b is:" << endl;
    b->foreach(print);
    cout << endl;
    cout << "Vector c is:" << endl;
    c->foreach(print);
    cout << endl;

    IVector* sum = IVector::add(a, b);
    if(!sum) cout << "Allocation or sum error" << endl;
    cout << "Vector a + b is:" << endl;
    sum->foreach(print);
    cout << endl;
    cout << "Vectors adding using inc() function" << endl;

    cout << "Is a+b equal c by Chebyshev norm with tolerance = 0.0001:";
    cout << IVector::equals(sum, c, IVector::NORM::CHEBYSHEV, 0.0001) << endl;

    cout << "Is a+b equal a by Chebyshev norm with tolerance = 0.0001:";
    cout << IVector::equals(sum, a, IVector::NORM::CHEBYSHEV, 0.0001) << endl;

    cout << "Is a+b equal a by Chebyshev norm with tolerance = 1:";
    cout << IVector::equals(sum, a, IVector::NORM::CHEBYSHEV, 10) << endl;
    cout << "equal() using norm() function" << endl;

    IVector* sub = IVector::sub(sum, b);
    if(!sub) cout << "Allocation or sub error" << endl;
    cout << "Is (a+b)-b equal a by Chebyshev norm with tolerance = 0.0001:";
    cout << IVector::equals(sub, a, IVector::NORM::CHEBYSHEV, 0.0001) << endl;
    cout << "Vectors adding using dec() function" << endl;

    cout << "Multiplying b on -1" << endl;
    b->scale(-1);
    cout << "Is -b equal a by Chebyshev norm with tolerance = 0.0001:";
    cout << IVector::equals(b, a, IVector::NORM::CHEBYSHEV, 0.0001) << endl;


    IVector* buf = IVector::createVector(4, a_data);
    cout << "Moving a in buf" << endl;
    IVector::moveInstance(buf, a);
    if(a == nullptr)
        cout << "a is void" << endl;
    cout << "Vector buf is:" << endl;
    buf->foreach(print);
    cout << endl;
    cout << "Moving using copyInstance() function" << endl;
    cout << "Cloning buf in a" << endl;
    a = buf->clone();
    if(a == nullptr)
        cout << "a is void" << endl;
    cout << "Vector buf is:" << endl;
    buf->foreach(print);
    cout << endl;
    cout << "Vector a is:" << endl;
    a->foreach(print);
    cout << endl;
    cout << "Cloning using getDim() and getData functions" << endl;

    double cord;
    a->getCord(3, cord);
    cout << "4th coordinate in vector a is " << cord << endl;

    cout << "Setting Infinity in a 4th coordinate" << endl;
    RC isValid = a->setCord(3, numeric_limits<double>::infinity());

    cout << "Dot product of a and b: " << IVector::dot(a, b) << endl;
    cout << "Chebyshev norm of a: " << a->norm(IVector::NORM::CHEBYSHEV)<< endl;
    cout << "First norm of a: " << a->norm(IVector::NORM::FIRST)<< endl;
    cout << "Second norm of a: " << a->norm(IVector::NORM::SECOND)<< endl;

    delete a;
    delete b;
    delete c;
    delete sum;
    delete buf;
    delete sub;
}


/*** Stanislav Kirpichenko set testing section ***/

constexpr size_t vecNum = 10;
constexpr size_t dim = 3;

void freeVecs(IVector** arr) {
    for (size_t i = 0; i < vecNum; i++) {
        delete arr[i];
    }
}

bool printSetContent(ISet* set, function<void(double)>& printer) {
    if (set->getSize() == 0) {
        cout << "{ empty set }" << endl;
        return true;
    }
    auto iter = set->getBegin();
    if (!iter) {
        return false;
    }
    double coords[dim] = { 0 };
    IVector* vec = IVector::createVector(dim, coords);
    if (!vec) {
        delete iter;
        return false;
    }
    cout << '{' << endl;
    while (iter->isValid()) {
        iter->getVectorCoords(vec);
        cout << "{ ";
        vec->foreach(printer);
        cout << '}' << endl;
        iter->next();
    }
    cout << '}' << endl;
    delete vec;
    delete iter;
    return true;
}

void setTest(ILogger* logger) {
    IVector* vecArray[vecNum] = { nullptr };
    double coords[dim] = { 0 };
    function<void(double)> print = [](double num) { std::cout << num << ' '; };
    IVector::setLogger(logger);
    ISet::setLogger(logger);
    for (int i = 0; i < vecNum; i++) {
        for (int j = 0; j < dim; j++) {
            coords[j] = i - j;
        }
        vecArray[i] = IVector::createVector(dim, coords);
        if (!vecArray[i]) {
            freeVecs(vecArray);
            cout << "Test failed: not enough heap memory" << endl;
            return;
        }
        cout << "vec_" << i << " = { ";
        vecArray[i]->foreach(print);
        cout << "}" << endl;
    }
    cout << "Inserting vectors to set1" << endl;
    ISet* set1 = ISet::createSet();
    if (!set1) {
        cout << "Test failed: not enough heap memory" << endl;
        freeVecs(vecArray);
        return;
    }
    for (auto vec : vecArray) {
        RC code = set1->insert(vec, IVector::NORM::FIRST, 1.0);
        if (code != RC::SUCCESS) {
            cout << "Test failed with error" << endl;
            logger->severe(code);
            freeVecs(vecArray);
            delete set1;
            return;
        }
    }
    cout << "Set1 content: ";
    if (!printSetContent(set1, print)) {
        cout << "set printing error, test is failed" << endl;
        freeVecs(vecArray);
        delete set1;
        return;
    }
    cout << "This output was generated using set iterator" << endl;
    cout << "Cloning set1 to set2" << endl;
    ISet* set2 = set1->clone();
    if (!set2) {
        cout << "Test failed: not enough heap memory" << endl;
        freeVecs(vecArray);
        delete set1;
        return;
    }
    cout << "Result of the operation \"is set1 equals to set2\": " << ISet::equals(set1, set2, IVector::NORM::SECOND, 0.01) << endl;
    cout << "Deleting all vectors that are similar to {7, 6, 5} with infinite norm and tolerance 1" << endl;
    set1->remove(vecArray[vecNum - 3], IVector::NORM::CHEBYSHEV, 1.0);
    cout << "Set1 content: ";
    if (!printSetContent(set1, print)) {
        cout << "set printing error, test is failed" << endl;
        freeVecs(vecArray);
        delete set1;
        delete set2;
        return;
    }
    cout << "Set2 content: ";
    if (!printSetContent(set2, print)) {
        cout << "set printing error, test is failed" << endl;
        freeVecs(vecArray);
        delete set1;
        delete set2;
        return;
    }
    cout << "Union of set1 and set2: ";
    ISet* set3 = ISet::makeUnion(set1, set2, IVector::NORM::SECOND, 0.01);
    if (!set3) {
        cout << "Test failed: not enough heap memory" << endl;
        freeVecs(vecArray);
        delete set1;
        delete set2;
        return;
    }
    if (!printSetContent(set3, print)) {
        cout << "set printing error, test is failed" << endl;
        freeVecs(vecArray);
        delete set1;
        delete set2;
        return;
    }
    delete set3;
    cout << "Intersection of set1 and set2: ";
    set3 = ISet::makeIntersection(set1, set2, IVector::NORM::SECOND, 0.01);
    if (!set3) {
        cout << "Test failed: not enough heap memory" << endl;
        freeVecs(vecArray);
        delete set1;
        delete set2;
        return;
    }
    if (!printSetContent(set3, print)) {
        cout << "set printing error, test is failed" << endl;
        freeVecs(vecArray);
        delete set1;
        delete set2;
        return;
    }
    delete set3;
    cout << "Difference of set1 and set2: ";
    set3 = ISet::sub(set1, set2, IVector::NORM::SECOND, 0.01);
    if (!set3) {
        cout << "Test failed: not enough heap memory" << endl;
        freeVecs(vecArray);
        delete set1;
        delete set2;
        return;
    }
    if (!printSetContent(set3, print)) {
        cout << "set printing error, test is failed" << endl;
        freeVecs(vecArray);
        delete set1;
        delete set2;
        return;
    }
    delete set3;
    cout << "Symmetric difference of set1 and set2: ";
    set3 = ISet::symSub(set1, set2, IVector::NORM::SECOND, 0.01);
    if (!set3) {
        cout << "Test failed: not enough heap memory" << endl;
        freeVecs(vecArray);
        delete set1;
        delete set2;
        return;
    }
    if (!printSetContent(set3, print)) {
        cout << "set printing error, test is failed" << endl;
        freeVecs(vecArray);
        delete set1;
        delete set2;
        return;
    }
    delete set3;
    cout << "Result of the operation \"is set1 subset of set2\": " << ISet::subSet(set1, set2, IVector::NORM::SECOND, 0.01) << endl;
    cout << "Adding vector {0, 0, 0} to set1" << endl;
    coords[0] = coords[1] = coords[2] = 0;
    vecArray[0]->setData(dim, coords);
    if (set1->insert(vecArray[0], IVector::NORM::SECOND, 0.01) != RC::SUCCESS) {
        cout << "Test failed: troubles with inserting" << endl;
        freeVecs(vecArray);
        delete set1;
        delete set2;
        return;
    }
    cout << "Result of the operation \"is set1 subset of set2\": " << ISet::subSet(set1, set2, IVector::NORM::SECOND, 0.01) << endl;
   /* cout << "Inserting 10000 vectors to set2" << endl;
    for (int i = 0; i < 10000; i++) {
        coords[0] = (i * i - 7 * i + 3) % 1301;
        coords[1] = (-3 * i * i + 255 * i + 32) % 1553;
        coords[2] = (i * i - 7833 * i - 98) % 1217;
        vecArray[0]->setData(dim, coords);
        RC code = set2->insert(vecArray[0], IVector::NORM::SECOND, 0.01);
        if (code != RC::SUCCESS) {
            logger->severe(code);
            break;
        }
        if (i % 1000 == 0) {
            cout << '.';
            cout.flush();
        }
    }
    cout << endl;
    */
    cout << "Set2 size: " << set2->getSize() << endl;
    cout << "Deleting set2" << endl;
    delete set2;
    cout << "Set1 content: ";
    if (!printSetContent(set1, print)) {
        cout << "\nSet printing error, test is failed" << endl;
        freeVecs(vecArray);
        delete set1;
        return;
    }
    cout << "Getting iterator to the element with index 3" << endl;
    auto iterator = set1->getIterator(3);
    if (!iterator) {
        cout << "\nUnable to get iterator. Test error" << endl;
        delete set1;
        freeVecs(vecArray);
        return;
    }
    cout << "Vector inside iterator: { ";
    IVector* vec = nullptr;
    if (iterator->getVectorCopy(vec) != RC::SUCCESS) {
        cout << "\nUnable to get vector from the iterator. Test error" << endl;
        delete set1;
        freeVecs(vecArray);
        return;
    }
    vec->foreach(print);
    cout << "}\nDeleting vector with index 3 from the set" << endl;
    set1->remove(3);
    cout << "Vector inside iterator: { ";
    if (iterator->getVectorCoords(vec) != RC::SUCCESS) {
        cout << "\nUnable to get vector from the iterator. Test error" << endl;
        delete set1;
        freeVecs(vecArray);
        return;
    }
    vec->foreach(print);
    cout << "}\nMoving iterator forward" << endl;
    iterator->next();
    cout << "Vector inside iterator: { ";
    if (iterator->getVectorCoords(vec) != RC::SUCCESS) {
        cout << "Unable to get vector from the iterator. Test error" << endl;
        delete set1;
        freeVecs(vecArray);
        return;
    }
    vec->foreach(print);
    cout << "}\nMoving iterator back" << endl;
    iterator->previous();
    cout << "Vector inside iterator: { ";
    if (iterator->getVectorCoords(vec) != RC::SUCCESS) {
        cout << "Unable to get vector from the iterator. Test error" << endl;
        delete set1;
        freeVecs(vecArray);
        return;
    }
    vec->foreach(print);
    cout << '}' << endl;
    cout << "Set1 content: ";
    if (!printSetContent(set1, print)) {
        cout << "\nSet printing error, test is failed" << endl;
        freeVecs(vecArray);
        delete iterator;
        delete vec;
        delete set1;
        return;
    }
    cout << "Printing only vectors with even indexes in set1 (using iterators): {" << endl;
    iterator->makeBegin();
    while (iterator->isValid()) {
        iterator->getVectorCoords(vec);
        cout << "{ ";
        vec->foreach(print);
        cout << "}" << endl;
        iterator->next(2);
    }
    cout << "}\n\nTest is finished successfully" << endl;
    delete iterator;
    delete vec;
    delete set1;
    freeVecs(vecArray);
}

bool PrintCompact(ICompact* compact, IMultiIndex* order) {
    function<void(double)> print = [](double num) { std::cout << num << ' '; };
    if (compact->getDim() == 0) {
        cout << "{ empty set }" << endl;
        return true;
    }
    auto iter = compact->getBegin(order);
    if (!iter) return false;


    double coords[dim] = { 0 };
    IVector* vec = IVector::createVector(dim, coords);
    if (!vec) {
        delete iter;
        return false;
    }
    cout << '{' << endl;
    while (!iter->isAtEnd()) {
        iter->getVectorCoords(vec);
        cout << "{ ";
        vec->foreach(print);
        cout << '}' << endl;
        iter->moveForwardIterator();
    }
    cout << '}' << endl;
    delete vec;
    delete iter;
    return true;
}

void TestingCompact(ILogger* logger) {
    IVector::setLogger(logger);
    IMultiIndex::setLogger(logger);
    ICompact::setLogger(logger);
    ICompact::IIterator::setLogger(logger);

    double* leftData = (double*)malloc(dim * sizeof(double));
    double* rightData = (double*)malloc(dim * sizeof(double));
    size_t* gridData = (size_t*)malloc(dim * sizeof(size_t));
    size_t* orderData = (size_t*)malloc(dim * sizeof(size_t));
    if (leftData == NULL || rightData == NULL || gridData == NULL) {
        free(leftData);
        free(rightData);
        free(gridData);
        free(orderData);
        if (logger) logger->warning(RC::ALLOCATION_ERROR);
        return;
    }
    for (size_t i = 0; i < dim; i++) {
        leftData[i] = 0.0f;
        rightData[i] = (double)dim;
        gridData[i] = (size_t)dim + 1;
        orderData[i] = (size_t)i;
    }
    IVector* left = IVector::createVector(dim, leftData);
    IVector* right = IVector::createVector(dim, rightData);
    IMultiIndex* grid = IMultiIndex::createMultiIndex(dim, gridData);
    IMultiIndex* order = IMultiIndex::createMultiIndex(dim, orderData);
    ICompact* compact1 = ICompact::createCompact(left, right, grid);

    cout << "Print Compact1 with leftBound = {0,0,0}, rightBound = {3,3,3}, grid = {4, 4, 4} and bypass order = {x->y->z}" << endl;
    PrintCompact(compact1, order);

    for (size_t i = 0; i < dim; i++)
        orderData[i] = (size_t)dim - i - 1;
    order->setData(dim, orderData);
    cout << "Print Compact1 with leftBound = {0,0,0}, rightBound = {3,3,3}, grid = {4, 4, 4} and bypass order = {z->y->x}" << endl;
    PrintCompact(compact1, order);

    double* isInData = (double*)malloc(dim * sizeof(double));
    for (size_t i = 0; i < dim; i++)
        isInData[i] = 2.0f;
    IVector* isIn = IVector::createVector(dim, isInData);
    cout << "is vector = {2,2,2} in compact: " << compact1->isInside(isIn) << endl;

    for (size_t i = 0; i < dim; i++)
        isInData[i] = 5.0f;
    isIn->setData(dim, isInData);
    cout << "is vector = {5,5,5} in compact: " << compact1->isInside(isIn) << endl;

    for (size_t i = 0; i < dim; i++)
        isInData[i] = 3.0f;
    isIn->setData(dim, isInData);
    cout << "is vector = {3,3,3} in compact: " << compact1->isInside(isIn) << endl;
    free(isInData);
    delete isIn;

    for (size_t i = 0; i < dim; i++) {
        leftData[i] = 2.0f;
        rightData[i] = (double)dim + 2.0f;
        orderData[i] = (size_t)i;
    }
    left->setData(dim, leftData);
    right->setData(dim, rightData);
    order->setData(dim, orderData);
    ICompact* compact2 = ICompact::createCompact(left, right, grid);
    cout << "Print Compact2 with leftBound = {2,2,2}, rightBound = {5,5,5}, grid = {4, 4, 4} and bypass order = {x->y->z}" << endl;
    PrintCompact(compact2, order);

    ICompact* intersection = ICompact::createIntersection(compact1, compact2, grid, 0.0001);
    cout << "Print intersection between Compact1 and Compact2 with grid = {4, 4, 4} and tolerance = 0.0001, and bypass order = {x->y->z}" << endl;
    PrintCompact(intersection, order);
    delete intersection;
    delete compact2;

    leftData[0] = 3.0f;
    left->setData(dim, leftData);
    compact2 = ICompact::createCompact(left, right, grid);
    cout << "Create Compact2 with leftBound = {3,2,2}, rightBound = {5,5,5}, grid = {4, 4, 4}..." << endl;
    intersection = ICompact::createIntersection(compact1, compact2, grid, 0.0001);
    cout << "Print intersection between Compact1 and Compact2 with grid = {4, 4, 4} and tolerance = 0.0001, and bypass order = {x->y->z}" << endl;
    PrintCompact(intersection, order);
    delete intersection;
    delete compact2;

    for (size_t i = 0; i < dim; i++)
        leftData[i] = 3.0f;
    left->setData(dim, leftData);
    compact2 = ICompact::createCompact(left, right, grid);
    cout << "Create Compact2 with leftBound = {4,4,4}, rightBound = {5,5,5}, grid = {4, 4, 4}..." << endl;
    intersection = ICompact::createIntersection(compact1, compact2, grid, 0.0001);
    cout << "Print intersection between Compact1 and Compact2 with grid = {4, 4, 4} and tolerance = 0.0001, and bypass order = {x->y->z}" << endl;
    PrintCompact(intersection, order);
    delete intersection;
    delete compact2;

    for (size_t i = 0; i < dim; i++)
        leftData[i] = 4.0f;
    left->setData(dim, leftData);
    compact2 = ICompact::createCompact(left, right, grid);
    cout << "Create Compact2 with leftBound = {4,4,4}, rightBound = {5,5,5}, grid = {4, 4, 4}..." << endl;
    intersection = ICompact::createIntersection(compact1, compact2, grid, 0.0001);
    cout << "Print intersection between Compact1 and Compact2 with grid = {4, 4, 4} and tolerance = 0.0001, and bypass order = {x->y->z}" << endl;
    if (intersection == nullptr)
        cout << "Intersection in nullptr" << endl;

    ICompact* span = ICompact::createCompactSpan(compact1, compact2, grid);
    cout << "Print intersection between Compact1 and Compact2 with grid = {4, 4, 4} and bypass order = {x->y->z}" << endl;
    PrintCompact(span, order);

    cout << "Printing is done using iterators" << endl;
    cout << "Interator methods uses getLeftBoundary(), getRightBoundary(), getVectorCopy(), getVectorCoords() methods" << endl;

    free(leftData);
    free(rightData);
    free(gridData);
    free(orderData);
    delete left;
    delete right;
    delete grid;
    delete order;
    delete compact1;
    delete compact2;
    delete intersection;
    delete span;
}


void TestingProblem(ILogger* logger) {
    const size_t paramsDim_ = 8;
    const size_t argsDim_ = 2;

    IVector::setLogger(logger);
    IMultiIndex::setLogger(logger);
    ICompact::setLogger(logger);
    ICompact::IIterator::setLogger(logger);

    HMODULE dll = LoadLibrary("../dlls/broker.dll");
    if (!dll) {
        cout << "Unable to load library" << endl;
        return;
    }
    BrokerProc bp = (BrokerProc)GetProcAddress(dll, "getBroker");
    if (!bp) {
        cout << "Unable to load procedure" << endl;
        FreeLibrary(dll);
        return;
    }
    IBroker* broker = (IBroker*)bp();



    double leftParamsData[paramsDim_] = { 1, 1, 0, 0, 0, 0, 0, 0 };
    double rightParamsData[paramsDim_] = { 5, 5, 5, 5, 5, 5, 5, 5 };
    size_t derOrderParamsData[paramsDim_] = { 0, 0, 0, 0, 0, 0, 0, 0 };
    double leftArgsData[argsDim_] = { -5, -5 };
    double rightArgsData[argsDim_] = { 5, 5 };
    size_t derOrderArgsData[argsDim_] = { 0, 0 };
    size_t paramsGridData[paramsDim_] = { 10, 10, 10, 10, 10, 10, 10, 10 };
    size_t argsGridData[paramsDim_] = { 10, 10 };

    IVector* paramsLeft = IVector::createVector(paramsDim_, leftParamsData);
    IVector* paramsRight = IVector::createVector(paramsDim_, rightParamsData);
    IMultiIndex* paramsGrid = IMultiIndex::createMultiIndex(paramsDim_, paramsGridData);
    IMultiIndex* derOrderParams = IMultiIndex::createMultiIndex(paramsDim_, derOrderParamsData);
    IVector* argsLeft = IVector::createVector(argsDim_, leftArgsData);
    IVector* argsRight = IVector::createVector(argsDim_, rightArgsData);
    IMultiIndex* argsGrid = IMultiIndex::createMultiIndex(argsDim_, argsGridData);
    IMultiIndex* derOrderArgs = IMultiIndex::createMultiIndex(argsDim_, derOrderArgsData);

    cout << "Creating problem Ax^2 + By^2 + Cx + Dy + Esin(Fx + Gy) + H, A > 0, B > 0" << endl;
    cout << "On compact [-5, 5] x [-5, 5]" << endl;
    ICompact* compactParams = ICompact::createCompact(paramsLeft, paramsRight, paramsGrid);
    ICompact* compactArgs = ICompact::createCompact(argsLeft, argsRight, argsGrid);
    IDiffProblem* problem = (IDiffProblem*)broker->getInterfaceImpl(IBroker::INTERFACE_IMPL::IDIFFPROBLEM);
    problem->setParamsDomain(compactParams);
    problem->setArgsDomain(compactArgs, logger);

    double paramsForSolveData[paramsDim_] = {2, 3, 0, 0, 0, 0, 0, 0};
    double argsForSolveData[argsDim_] = { 1, 1 };

    IVector* paramsForSolve = IVector::createVector(paramsDim_, paramsForSolveData);
    IVector* argsForSolve = IVector::createVector(argsDim_, argsForSolveData);
    cout << "Stating f(x, y) = 2x^2 + 3y^2" << endl;
    problem->setParams(paramsForSolve);
    cout << "f(1,1) = " << problem->evalByArgs(argsForSolve) << endl;

    argsForSolve->setCord(0, 5.0f);
    argsForSolve->setCord(1, 5.0f);
    cout << "f(5,5) = " << problem->evalByArgs(argsForSolve) << endl;

    argsForSolve->setCord(0, 6.0f);
    argsForSolve->setCord(1, 6.0f);
    cout << "f(6,6) = ";
    problem->evalByArgs(argsForSolve);
    cout << endl;

    argsForSolve->setCord(0, 5.0f);
    argsForSolve->setCord(1, 5.0f);
    derOrderArgs->setAxisIndex(0, 1);
    cout << "df/dx(5,5) = " << problem->evalDerivativeByArgs(argsForSolve, derOrderArgs) << endl;

    derOrderArgs->setAxisIndex(0, 2);
    cout << "d^2f/dx^2(5,5) = " << problem->evalDerivativeByArgs(argsForSolve, derOrderArgs) << endl;

    derOrderArgs->setAxisIndex(0, 1);
    derOrderArgs->setAxisIndex(1, 1);
    cout << "d^2f/dxdy(5,5) = " << problem->evalDerivativeByArgs(argsForSolve, derOrderArgs) << endl;

    derOrderArgs->setAxisIndex(0, 0);
    derOrderArgs->setAxisIndex(1, 2);
    cout << "d^2f/dy^2(5,5) = " << problem->evalDerivativeByArgs(argsForSolve, derOrderArgs) << endl;

    argsForSolve->setCord(0, 6.0f);
    argsForSolve->setCord(1, 6.0f);
    cout << "d^2f/dy^2(6,6) = ";
    problem->evalDerivativeByArgs(argsForSolve, derOrderArgs);
    cout << endl;

    argsForSolve->setCord(0, 5.0f);
    argsForSolve->setCord(1, 5.0f);
    IVector* gradientArgs = argsForSolve->clone();
    function<void(double)> print = [](double num) { std::cout << num << ' '; };
    problem->evalGradientByArgs(argsForSolve, gradientArgs);
    cout << "GradF(5,5) = {";
    gradientArgs->foreach(print);
    cout << "}" << endl;

    paramsForSolve->setData(paramsDim_, rightParamsData);
    IVector* gradientParams = paramsForSolve->clone();
    problem->evalGradientByParams(paramsForSolve, gradientParams);
    cout << "GradbyParamsF(5,5,5,5,5,5,5,5) = {";
    gradientParams->foreach(print);
    cout << "}" << endl;

    delete paramsLeft;
    delete paramsRight;
    delete paramsGrid;
    delete derOrderParams;
    delete compactParams;
    delete argsLeft;
    delete argsRight;
    delete argsGrid;
    delete derOrderArgs;
    delete compactArgs;
    delete paramsForSolve;
    delete argsForSolve;
    delete problem;
    delete gradientArgs;
    delete gradientParams;
    broker->release();
    FreeLibrary(dll);
}

void TestingSolver(ILogger* logger) {
    const size_t solverParamsDim_ = 4;
    const size_t paramsDim_ = 8;
    const size_t argsDim_ = 2;
    function<void(double)> print = [](double num) { std::cout << num << ' '; };

    IVector::setLogger(logger);
    IMultiIndex::setLogger(logger);
    ICompact::setLogger(logger);
    ICompact::IIterator::setLogger(logger);


    HMODULE dll = LoadLibrary("../dlls/broker.dll");
    if (!dll) {
        cout << "Unable to load library" << endl;
        return;
    }
    BrokerProc bp = (BrokerProc)GetProcAddress(dll, "getBroker");
    if (!bp) {
        cout << "Unable to load procedure" << endl;
        FreeLibrary(dll);
        return;
    }
    IBroker* broker = (IBroker*)bp();


    double leftParamsData[paramsDim_] = { 1, 1, 0, 0, 0, 0, 0, 0 };
    double rightParamsData[paramsDim_] = { 5, 5, 5, 5, 5, 5, 5, 5 };
    double leftArgsData[argsDim_] = { -5, -5 };
    double rightArgsData[argsDim_] = { 5, 5 };
    size_t paramsGridData[paramsDim_] = { 10, 10, 10, 10, 10, 10, 10, 10 };
    size_t argsGridData[paramsDim_] = { 10, 10 };

    double leftSolverParamsData[solverParamsDim_] = { 0.001, 0.001, 0.001, 0.000001 };
    double rightSolverParamsData[solverParamsDim_] = { 0.9, 0.9, 0.9, 0.9 };
    size_t solverParamsGridData[solverParamsDim_] = { 0, 0, 0, 0 };

    IVector* paramsLeft = IVector::createVector(paramsDim_, leftParamsData);
    IVector* paramsRight = IVector::createVector(paramsDim_, rightParamsData);
    IMultiIndex* paramsGrid = IMultiIndex::createMultiIndex(paramsDim_, paramsGridData);
    IVector* argsLeft = IVector::createVector(argsDim_, leftArgsData);
    IVector* argsRight = IVector::createVector(argsDim_, rightArgsData);
    IMultiIndex* argsGrid = IMultiIndex::createMultiIndex(argsDim_, argsGridData);
    IVector* solverParamsLeft = IVector::createVector(solverParamsDim_, leftSolverParamsData);
    IVector* solverParamsRight = IVector::createVector(solverParamsDim_, rightSolverParamsData);
    IMultiIndex* solverParamsGrid = IMultiIndex::createMultiIndex(solverParamsDim_, solverParamsGridData);

    cout << "Creating problem Ax^2 + By^2 + Cx + Dy + Esin(Fx + Gy) + H, A > 0, B > 0" << endl;
    cout << "On compact [-5, 5] x [-5, 5]" << endl;
    ICompact* compactParams = ICompact::createCompact(paramsLeft, paramsRight, paramsGrid);
    ICompact* compactArgs = ICompact::createCompact(argsLeft, argsRight, argsGrid);
    ICompact* compactSolverParams = ICompact::createCompact(solverParamsLeft, solverParamsRight, solverParamsGrid);
    IDiffProblem* problem = (IDiffProblem*)broker->getInterfaceImpl(IBroker::INTERFACE_IMPL::IDIFFPROBLEM);
    problem->setParamsDomain(compactParams);
    problem->setArgsDomain(compactArgs, logger);

    double paramsForProblemData[paramsDim_] = { 1, 3, 3, 2, 1, 3, 1, 0 };
    IVector* paramsForProblem = IVector::createVector(paramsDim_, paramsForProblemData);
    cout << "Stating f(x, y) = x^2 + 3y^2 + 3x + 2y + sin(3x + y)" << endl;
    problem->setParams(paramsForProblem);

    ISolver* solver = (ISolver*)broker->getInterfaceImpl(IBroker::INTERFACE_IMPL::ISOLVER);
    solver->setProblem(problem);
    solver->setArgsDomain(compactArgs, logger);
    solver->setParamsDomain(compactSolverParams);


    double x0_Data[argsDim_] = { -3.0f, 2.63 };
    double InitParamsData[solverParamsDim_] = { 0.5, 0.5, 0.5, 0.0001 };
    IVector* x0 = IVector::createVector(argsDim_, x0_Data);
    IVector* InitParams = IVector::createVector(solverParamsDim_, InitParamsData);
    solver->solveByArgs(x0, InitParams);
    IVector* solution = nullptr;
    solver->getSolution(solution);
    cout << "Solution = {";
    solution->foreach(print);
    cout << "}" << endl;

    delete paramsLeft;
    delete paramsRight;
    delete paramsGrid;
    delete argsLeft;
    delete argsRight;
    delete argsGrid;
    delete solverParamsLeft;
    delete solverParamsRight;
    delete solverParamsGrid;
    delete compactParams;
    delete compactArgs;
    delete compactSolverParams;
    delete problem;
    delete paramsForProblem;
    delete solver;
    delete x0;
    delete InitParams;
    delete solution;
    broker->release();
    FreeLibrary(dll);
}

int main() {
    ILogger* logger = ILogger::createLogger();
    testingVector(logger);
    setTest(logger);
    TestingCompact(logger);
    TestingProblem(logger);
    TestingSolver(logger);
    delete logger;
    return 0;
}
