#include "mvardie.h"

#include <stdlib.h>
#include <stdarg.h>
#include <stdio.h>

void mvar_die(const char *format, ...)
{
    va_list ap;

    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);

    fputc('\n', stderr);

    abort();
}
