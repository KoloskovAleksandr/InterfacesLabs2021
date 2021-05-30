#include "KoloskovLogger.h"
#include <new>

Logger::Logger(const char* const& filename, bool overwrite) {
	if (!filename)	stream = stdout;
	else {
		stream = fopen(filename, overwrite ? "w" : "a");
		if (!stream)	stream = stdout;
	}
}

ILogger* ILogger::createLogger() {
	return (ILogger*)new(std::nothrow) Logger(nullptr);
}

ILogger* ILogger::createLogger(const char* const& filename, bool overwrite) {
	return (ILogger*)new(std::nothrow) Logger(filename, overwrite);
}

RC Logger::log(RC code, Level level) {
	fprintf(stream, "%s: %s\n", rcLevel(level), rcMessage(code));
	return RC::SUCCESS;
}

RC Logger::log(RC code, Level level, const char* const& srcfile, const char* const& function, int line){
	log(code, level);
	fprintf(stream, "File: %s, Function: %s, Line: %i\n", srcfile, function, line);
	return RC::SUCCESS;
}

Logger::~Logger() {
	fclose(stream);
}

ILogger::~ILogger() {};
