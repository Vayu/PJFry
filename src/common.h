/*
 * common.h - common includes for PJFry library
 *
 * Valery Yundin <valery.yundin@desy.de>
 */

#ifndef QUL_COMMON_H
#define QUL_COMMON_H

#ifdef HAVE_CONFIG_H
#   include <config.h>
#endif /* HAVE_CONFIG_H */

#include <complex>
#include <limits>
#include <inttypes.h>
// #include <stdint.h>

#include <cassert>
#ifndef NDEBUG
#include <cstdio> // DEBUG
#endif

typedef std::complex<double> ncomplex;
// typedef long long int int64_t;

// Forward declarations
class ICache;

class MCache;
class MCache5;

class Minor5;

class Kinem5;
class Kinem4;
class Kinem3;
class Kinem2;

#define CONST __attribute__ ((const))
#define PURE __attribute__ ((pure))

#endif /* QUL_COMMON_H */

