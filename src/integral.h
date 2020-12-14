/*
 * integral.h - scalar integrals
 *
 * this file is part of PJFry library
 * Copyright 2011 Valery Yundin
 */

#ifndef QUL_INTEGRAL_H
#define QUL_INTEGRAL_H

#include "common.h"
#include "kinem.h"
#include "cache.h"

class Initialize {
public:
  Initialize();
  ~Initialize();
};

ICache::Ival qlI1(const Kinem1& k);
ICache::Ival qlI2(const Kinem2& k);
ICache::Ival qlI3(const Kinem3& k);
ICache::Ival qlI4(const Kinem4& k);

#if !defined(DPKIND)
// Fortran DOUBLE PRECISION / REAL*8 / KIND=8 / kind(1d0) / selected_real_kind(15)
#define DPKIND 8
#endif /* !defined(DPKIND) */

#if DPKIND == 4
// Fortran REAL / REAL*4 / KIND=4 / kind(1.0) / selected_real_kind(6) / -freal-8-real-4
#define DPKIND_FLOAT float
#define DPKIND_COMPLEX std::complex<float>
#elif DPKIND == 8
// Fortran DOUBLE PRECISION / REAL*8 / KIND=8 / kind(1d0) / selected_real_kind(15)
#define DPKIND_FLOAT double
#if defined(USE_QCDLOOP2)
#define DPKIND_COMPLEX ql::complex
#else /* defined(USE_QCDLOOP2) */
#define DPKIND_COMPLEX std::complex<double>
#endif /* defined(USE_QCDLOOP2) */
#elif DPKIND == 10
// Fortran REAL*10 / KIND=10 / selected_real_kind(18) / -freal-8-real-10
#define DPKIND_FLOAT long double
#if defined(USE_QCDLOOP) || defined(USE_ONELOOP)
#include <complex.h> // ISO C99 _Complex
#define DPKIND_COMPLEX long double _Complex
#else /* defined(USE_QCDLOOP) || defined(USE_ONELOOP) */
#define DPKIND_COMPLEX std::complex<long double>
#endif /* defined(USE_QCDLOOP) || defined(USE_ONELOOP) */
#elif DPKIND == 16
// Fortran REAL*16 / KIND=16 / kind(1q0) / selected_real_kind(33) / -freal-8-real-16
#if defined(USE_QCDLOOP2)
#define DPKIND_FLOAT ql::qdouble
#define DPKIND_COMPLEX ql::qcomplex
#else /* defined(USE_QCDLOOP2) */
extern "C" { /* for gcc4.7 compatibility */
#include <quadmath.h>
} /* extern "C" */
#define DPKIND_FLOAT __float128
#define DPKIND_COMPLEX __complex128
// namespace std { inline __float128 sqrt(__float128 x) { return sqrtq(x); } }
// namespace std { inline ostream& operator<<(ostream& o, __float128 x) { /* a brutal hack */ o << ((long double)x); return o; } }
// namespace std { inline ostream& operator<<(ostream& o, __complex128 z) { o << "(" << crealq(z) << "," << cimagq(z) << ")"; return o; } }
#endif /* defined(USE_QCDLOOP2) */
#endif /* DPKIND == 16 */

//
// QCDLoop 1.x / FF or QCDLoop1 / FF
//
#ifdef USE_QCDLOOP
#define qlprec F77_FUNC(ffprec,FFPREC)
#define qlflag F77_FUNC(ffflag,FFFLAG)
extern "C" {
  void F77_FUNC(qlinit,QLINIT)();
  void F77_FUNC(ffexi,FFEXI)();
  
#ifdef USE_F2C
  void F77_FUNC(qli1,QLI1)(DPKIND_COMPLEX *rslt, DPKIND_FLOAT *m1, DPKIND_FLOAT *mu2, int *ep);
  void F77_FUNC(qli2,QLI2)(DPKIND_COMPLEX *rslt, DPKIND_FLOAT *p1, DPKIND_FLOAT *m1, DPKIND_FLOAT *m2, DPKIND_FLOAT *mu2, int *ep);
  void F77_FUNC(qli3,QLI3)(DPKIND_COMPLEX *rslt, DPKIND_FLOAT *p1, DPKIND_FLOAT *p2, DPKIND_FLOAT *p3, DPKIND_FLOAT *m1, DPKIND_FLOAT *m2, DPKIND_FLOAT *m3, DPKIND_FLOAT *mu2, int *ep);
  void F77_FUNC(qli4,QLI4)(DPKIND_COMPLEX *rslt, DPKIND_FLOAT *p1, DPKIND_FLOAT *p2, DPKIND_FLOAT *p3, DPKIND_FLOAT *p4, DPKIND_FLOAT *s12, DPKIND_FLOAT *s23, DPKIND_FLOAT *m1, DPKIND_FLOAT *m2, DPKIND_FLOAT *m3, DPKIND_FLOAT *m4, DPKIND_FLOAT *mu2, int *ep);
#else
  DPKIND_COMPLEX F77_FUNC(qli1,QLI1)(DPKIND_FLOAT *m1, DPKIND_FLOAT *mu2, int *ep);
  DPKIND_COMPLEX F77_FUNC(qli2,QLI2)(DPKIND_FLOAT *p1, DPKIND_FLOAT *m1, DPKIND_FLOAT *m2, DPKIND_FLOAT *mu2, int *ep);
  DPKIND_COMPLEX F77_FUNC(qli3,QLI3)(DPKIND_FLOAT *p1, DPKIND_FLOAT *p2, DPKIND_FLOAT *p3, DPKIND_FLOAT *m1, DPKIND_FLOAT *m2, DPKIND_FLOAT *m3, DPKIND_FLOAT *mu2, int *ep);
  DPKIND_COMPLEX F77_FUNC(qli4,QLI4)(DPKIND_FLOAT *p1, DPKIND_FLOAT *p2, DPKIND_FLOAT *p3, DPKIND_FLOAT *p4, DPKIND_FLOAT *s12, DPKIND_FLOAT *s23, DPKIND_FLOAT *m1, DPKIND_FLOAT *m2, DPKIND_FLOAT *m3, DPKIND_FLOAT *m4, DPKIND_FLOAT *mu2, int *ep);
#endif
  
  extern struct {
    DPKIND_FLOAT xloss, precx, precc, xalogm, xclogm, xalog2, xclog2, reqprc;
  } qlprec;
  
  extern struct {
    int lwrite, ltest, l4also, ldc3c4, lmem, lwarn, ldot,
      nevent, ner, id, idsub, nwidth, nschem, onshel, idot;
  } qlflag;
} /* extern "C" */
#endif /* USE_QCDLOOP */

//
// OneLOop
//
#ifdef USE_ONELOOP
extern "C" {
  void F77_FUNC_(avh_olo_mu_set,AVH_OLO_MU_SET)(DPKIND_FLOAT *mu);
  void F77_FUNC_(avh_olo_onshell,AVH_OLO_ONSHELL)(DPKIND_FLOAT *thrs);
  /* all below functions return three values in "DPKIND_COMPLEX rslt[3]" */
  void F77_FUNC_(avh_olo_a0m,AVH_OLO_A0M)(DPKIND_COMPLEX *rslt, DPKIND_FLOAT *m1);
  void F77_FUNC_(avh_olo_b0m,AVH_OLO_B0M)(DPKIND_COMPLEX *rslt, DPKIND_FLOAT *p1, DPKIND_FLOAT *m1, DPKIND_FLOAT *m2);
  void F77_FUNC_(avh_olo_c0m,AVH_OLO_C0M)(DPKIND_COMPLEX *rslt, DPKIND_FLOAT *p1, DPKIND_FLOAT *p2, DPKIND_FLOAT *p3, DPKIND_FLOAT *m1, DPKIND_FLOAT *m2, DPKIND_FLOAT *m3);
  void F77_FUNC_(avh_olo_d0m,AVH_OLO_D0M)(DPKIND_COMPLEX *rslt, DPKIND_FLOAT *p1, DPKIND_FLOAT *p2, DPKIND_FLOAT *p3, DPKIND_FLOAT *p4, DPKIND_FLOAT *s12, DPKIND_FLOAT *s23, DPKIND_FLOAT *m1, DPKIND_FLOAT *m2, DPKIND_FLOAT *m3, DPKIND_FLOAT *m4);
} /* extern "C" */
#endif /* USE_ONELOOP */

#endif /* QUL_INTEGRAL_H */
