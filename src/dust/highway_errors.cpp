/// Highway errors as C++ exceptions (instead of a message on stderr and abort)

#include <cstdarg>
#include <cstdio>
#include <stdexcept>
#include <string>

#include <hwy/base.h>

namespace hwy {

HWY_DLLEXPORT void Warn(const char*, int, const char*, ...) {}

HWY_DLLEXPORT HWY_NORETURN void Abort(const char*, int, const char* format, ...)
{
  char message[800];
  va_list args;
  va_start(args, format);
  vsnprintf(message, sizeof(message), format, args);
  va_end(args);
  throw std::runtime_error(std::string("Highway: ") + message);
}

} // namespace hwy
