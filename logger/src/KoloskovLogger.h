#pragma once
#include <cstdio>
#include "../include/ILogger.h"

class Logger : public ILogger {
private:
    FILE* stream;
    const char* rcMessage(RC code) {
        switch (code) {
        case RC::UNKNOWN:
            return "Unknown error";
        case RC::SUCCESS:
            return "Operation successfull";
        case RC::INVALID_ARGUMENT:
            return "Invalid argument";
        case RC::MISMATCHING_DIMENSIONS:
            return "Mismatching dimensions";
        case RC::INDEX_OUT_OF_BOUND:
            return "Index is out of bounds";
        case RC::INFINITY_OVERFLOW:
            return "Result is greater than infinity";
        case RC::NOT_NUMBER:
            return "Result is undefined";
        case RC::ALLOCATION_ERROR:
            return "Memory allocation error";
        case RC::NULLPTR_ERROR:
            return "Nullptr as argument";
        case RC::FILE_NOT_FOUND:
            return "File not found";
        case RC::VECTOR_NOT_FOUND:
            return "Vector not found";
        case RC::IO_ERROR:
            return "Reading or writing error";
        case RC::MEMORY_INTERSECTION:
            return "Found intersecting memory while copying instance";
        default:
            return "";
        }
    }
    const char* rcLevel(Level level) {
        switch (level) {
        case Level::SEVERE:
            return "SEVERE";
        case Level::WARNING:
            return "WARNING";
        case Level::INFO:
            return "INFO";
        default:
            return "";
        }
    }
    
public:
    Logger(const char* const& filename, bool overwrite = true);
    RC log(RC code, Level level, const char* const& srcfile, const char* const& function, int line) override;
    RC log(RC code, Level level) override;
    ~Logger();
};