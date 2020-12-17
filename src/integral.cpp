/*
 * integral.cpp - scalar integrals wrappers
 *
 * this file is part of PJFry library
 * Copyright 2011 Valery Yundin
 *
 * QCDLoop 2.x interface by Jacek M. Holeczek 2020.01.31
 */

#include "integral.h"
#include "minor.h"
#include "cache.h"

#include <cmath> // std::sqrt
#include <iostream>

static Initialize pjinit;

//
// QCDLoop 2.x
//
#ifdef USE_QCDLOOP2
// the same identifier is defined in "PJFry/config.h" and "qcdloop/config.h"
// #if defined(VERSION)
// #undef VERSION
// #endif /* defined(VERSION) */

#include "qcdloop/qcdloop.h" // needs "-std=c++11" or newer
#include <stdexcept>

Initialize::Initialize() {
#ifndef NDEBUG
  std::cout << std::endl << "PJFRY init" << std::endl;
#endif
}

Initialize::~Initialize() { /* NOOP */ }

ICache::Ival qlI1(const Kinem1 &k) {
  ICache::Ival ivalue;
  
  static ql::TadPole<DPKIND_COMPLEX, DPKIND_FLOAT, DPKIND_FLOAT> td;
  static std::vector<DPKIND_FLOAT> mI1(1);
  static std::vector<DPKIND_COMPLEX> r(3);
  
  try {
    mI1[0] = k.m1();
    td.integral(r, ICache::getMu2(), mI1);
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
  
  ivalue.val[0] = r[0];
  ivalue.val[1] = r[1];
  ivalue.val[2] = r[2];
  return ivalue;
}

ICache::Ival qlI2(const Kinem2 &k) {
  ICache::Ival ivalue;
  
  static ql::Bubble<DPKIND_COMPLEX, DPKIND_FLOAT, DPKIND_FLOAT> bb;
  static std::vector<DPKIND_FLOAT> mI2(2);
  static std::vector<DPKIND_FLOAT> pI2(1);
  static std::vector<DPKIND_COMPLEX> r(3);
  
  try {
    mI2[0] = k.m1();
    mI2[1] = k.m2();
    pI2[0] = k.p1();
    bb.integral(r, ICache::getMu2(), mI2, pI2);
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
  
  ivalue.val[0] = r[0];
  ivalue.val[1] = r[1];
  ivalue.val[2] = r[2];
  return ivalue;
}

ICache::Ival qlI3(const Kinem3 &k) {
  ICache::Ival ivalue;
  
  static ql::Triangle<DPKIND_COMPLEX, DPKIND_FLOAT, DPKIND_FLOAT> tr;
  static std::vector<DPKIND_FLOAT> mI3(3);
  static std::vector<DPKIND_FLOAT> pI3(3);
  static std::vector<DPKIND_COMPLEX> r(3);
  
  try {
    mI3[0] = k.m1();
    mI3[1] = k.m2();
    mI3[2] = k.m3();
    pI3[0] = k.p1();
    pI3[1] = k.p2();
    pI3[2] = k.p3();
    tr.integral(r, ICache::getMu2(), mI3, pI3);
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
  
  ivalue.val[0] = r[0];
  ivalue.val[1] = r[1];
  ivalue.val[2] = r[2];
  return ivalue;
}

ICache::Ival qlI4(const Kinem4 &k) {
  ICache::Ival ivalue;
  
  static ql::Box<DPKIND_COMPLEX, DPKIND_FLOAT, DPKIND_FLOAT> bo;
  static std::vector<DPKIND_FLOAT> mI4(4);
  static std::vector<DPKIND_FLOAT> pI4(6);
  static std::vector<DPKIND_COMPLEX> r(3);
  
  try {
    if (k.s12() == 0. || k.s23() == 0.) {
      if (k.p1() != 0. && k.p3() != 0.) {
        mI4[0] = k.m1();
        mI4[1] = k.m3();
        mI4[2] = k.m2();
        mI4[3] = k.m4();
        pI4[0] = k.s12();
        pI4[1] = k.p2();
        pI4[2] = k.s23();
        pI4[3] = k.p4();
        pI4[4] = k.p1();
        pI4[5] = k.p3();
      } else if (k.p2() != 0. && k.p4() != 0.) {
        mI4[0] = k.m2();
        mI4[1] = k.m4();
        mI4[2] = k.m3();
        mI4[3] = k.m1();
        pI4[0] = k.s23();
        pI4[1] = k.p3();
        pI4[2] = k.s12();
        pI4[3] = k.p1();
        pI4[4] = k.p2();
        pI4[5] = k.p4();
      } else { assert(0); }
    } else {
      mI4[0] = k.m1();
      mI4[1] = k.m2();
      mI4[2] = k.m3();
      mI4[3] = k.m4();
      pI4[0] = k.p1();
      pI4[1] = k.p2();
      pI4[2] = k.p3();
      pI4[3] = k.p4();
      pI4[4] = k.s12();
      pI4[5] = k.s23();
    }
    bo.integral(r, ICache::getMu2(), mI4, pI4);
  } catch (std::exception &e) {
    std::cout << e.what() << std::endl;
    exit(-1);
  }
  
  ivalue.val[0] = r[0];
  ivalue.val[1] = r[1];
  ivalue.val[2] = r[2];
  return ivalue;
}
#endif /* USE_QCDLOOP2 */

//
// QCDLoop 1.x / FF or QCDLoop1 / FF
//
#ifdef USE_QCDLOOP
Initialize::Initialize() {
#ifndef NDEBUG
  std::cout << "PJFRY init" << std::endl;
#endif
  
  F77_FUNC(qlinit,QLINIT)();
  
  // QCDLoop 1.x / FF (and thus also QCDLoop1) include some
  // CERNLIB routines so, no need to go beyond "long double"
#if DPKIND == 10 || DPKIND == 16
  // const long double dbl_min = std::numeric_limits<double>::min();
  const long double dbl_min = std::numeric_limits<long double>::min();
#else /* DPKIND == 10 || DPKIND == 16 */
  const DPKIND_FLOAT dbl_min = std::numeric_limits<DPKIND_FLOAT>::min();
#endif /* DPKIND == 10 || DPKIND == 16 */
  
  if (qlprec.xalogm < dbl_min) {
    qlprec.xalogm = dbl_min;
    qlprec.xalog2 = std::sqrt(dbl_min);
#ifndef NDEBUG
    std::cout << "Set xalogm to normalized value " << dbl_min << std::endl;
    std::cout << "Set xalog2 to normalized value " << std::sqrt(dbl_min) << std::endl;
#endif
  }
  if (qlprec.xclogm < dbl_min) {
    qlprec.xclogm = dbl_min;
    qlprec.xclog2 = std::sqrt(dbl_min);
#ifndef NDEBUG
    std::cout << "Set xclogm to normalized value " << dbl_min << std::endl;
    std::cout << "Set xclog2 to normalized value " << std::sqrt(dbl_min) << std::endl;
#endif
  }
  if (qlflag.lwarn) {
    qlflag.lwarn = 0;
#ifndef NDEBUG
    std::cout << "Disable FF warnings " << qlflag.lwarn << std::endl;
#endif
  }
}

Initialize::~Initialize() {
  F77_FUNC(ffexi,FFEXI)();
}

ICache::Ival qlI1(const Kinem1& k) {
  ICache::Ival ivalue;
  
  int ep;
  DPKIND_FLOAT m1  = k.m1();
  DPKIND_FLOAT mu2 = ICache::getMu2();
  
#ifdef USE_F2C
  DPKIND_COMPLEX rslt;
  ep = 0;
  F77_FUNC(qli1,QLI1)(&rslt, &m1, &mu2, &ep);
  ivalue.val[0] = rslt;
  ep = -1;
  F77_FUNC(qli1,QLI1)(&rslt, &m1, &mu2, &ep);
  ivalue.val[1] = rslt;
  ep = -2;
  F77_FUNC(qli1,QLI1)(&rslt, &m1, &mu2, &ep);
  ivalue.val[2] = rslt;
#else
  ep = 0;
  ivalue.val[0] = F77_FUNC(qli1,QLI1)(&m1, &mu2, &ep);
  ep = -1;
  ivalue.val[1] = F77_FUNC(qli1,QLI1)(&m1, &mu2, &ep);
  ep = -2;
  ivalue.val[2] = F77_FUNC(qli1,QLI1)(&m1, &mu2, &ep);
#endif
  
  return ivalue;
}

ICache::Ival qlI2(const Kinem2& k) {
  ICache::Ival ivalue;
  
  int ep;
  DPKIND_FLOAT p1  = k.p1();
  DPKIND_FLOAT m1  = k.m1();
  DPKIND_FLOAT m2  = k.m2();
  DPKIND_FLOAT mu2 = ICache::getMu2();
  
#ifdef USE_F2C
  DPKIND_COMPLEX rslt;
  ep = 0;
  F77_FUNC(qli2,QLI2)(&rslt, &p1, &m1, &m2, &mu2, &ep);
  ivalue.val[0] = rslt;
  ep = -1;
  F77_FUNC(qli2,QLI2)(&rslt, &p1, &m1, &m2, &mu2, &ep);
  ivalue.val[1] = rslt;
  ep = -2;
  F77_FUNC(qli2,QLI2)(&rslt, &p1, &m1, &m2, &mu2, &ep);
  ivalue.val[2] = rslt;
#else
  ep = 0;
  ivalue.val[0] = F77_FUNC(qli2,QLI2)(&p1, &m1, &m2, &mu2, &ep);
  ep = -1;
  ivalue.val[1] = F77_FUNC(qli2,QLI2)(&p1, &m1, &m2, &mu2, &ep);
  ep = -2;
  ivalue.val[2] = F77_FUNC(qli2,QLI2)(&p1, &m1, &m2, &mu2, &ep);
#endif
  
  return ivalue;
}

ICache::Ival qlI3(const Kinem3& k) {
  ICache::Ival ivalue;
  
  int ep;
  DPKIND_FLOAT p1  = k.p1();
  DPKIND_FLOAT p2  = k.p2();
  DPKIND_FLOAT p3  = k.p3();
  DPKIND_FLOAT m1  = k.m1();
  DPKIND_FLOAT m2  = k.m2();
  DPKIND_FLOAT m3  = k.m3();
  DPKIND_FLOAT mu2 = ICache::getMu2();
  
#ifdef USE_F2C
  DPKIND_COMPLEX rslt;
  ep = 0;
  F77_FUNC(qli3,QLI3)(&rslt, &p1, &p2, &p3, &m1, &m2, &m3, &mu2, &ep);
  ivalue.val[0] = rslt;
  ep = -1;
  F77_FUNC(qli3,QLI3)(&rslt, &p1, &p2, &p3, &m1, &m2, &m3, &mu2, &ep);
  ivalue.val[1] = rslt;
  ep = -2;
  F77_FUNC(qli3,QLI3)(&rslt, &p1, &p2, &p3, &m1, &m2, &m3, &mu2, &ep);
  ivalue.val[2] = rslt;
#else
  ep = 0;
  ivalue.val[0] = F77_FUNC(qli3,QLI3)(&p1, &p2, &p3, &m1, &m2, &m3, &mu2, &ep);
  ep = -1;
  ivalue.val[1] = F77_FUNC(qli3,QLI3)(&p1, &p2, &p3, &m1, &m2, &m3, &mu2, &ep);
  ep = -2;
  ivalue.val[2] = F77_FUNC(qli3,QLI3)(&p1, &p2, &p3, &m1, &m2, &m3, &mu2, &ep);
#endif
  
  return ivalue;
}

ICache::Ival qlI4(const Kinem4& k) {
  ICache::Ival ivalue;
  
  int ep;
  DPKIND_FLOAT p1  = k.p1();
  DPKIND_FLOAT p2  = k.p2();
  DPKIND_FLOAT p3  = k.p3();
  DPKIND_FLOAT p4  = k.p4();
  DPKIND_FLOAT s12 = k.s12();
  DPKIND_FLOAT s23 = k.s23();
  DPKIND_FLOAT m1  = k.m1();
  DPKIND_FLOAT m2  = k.m2();
  DPKIND_FLOAT m3  = k.m3();
  DPKIND_FLOAT m4  = k.m4();
  DPKIND_FLOAT mu2 = ICache::getMu2();
  
  if (s12 == 0. || s23 == 0.) {
    if (p1 != 0. && p3 != 0.) {
#ifdef USE_F2C
      DPKIND_COMPLEX rslt;
      ep = 0;
      F77_FUNC(qli4,QLI4)(&rslt, &s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4, &mu2, &ep);
      ivalue.val[0] = rslt;
      ep = -1;
      F77_FUNC(qli4,QLI4)(&rslt, &s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4, &mu2, &ep);
      ivalue.val[1] = rslt;
      ep = -2;
      F77_FUNC(qli4,QLI4)(&rslt, &s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4, &mu2, &ep);
      ivalue.val[2] = rslt;
#else
      ep = 0;
      ivalue.val[0] = F77_FUNC(qli4,QLI4)(&s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4, &mu2, &ep);
      ep = -1;
      ivalue.val[1] = F77_FUNC(qli4,QLI4)(&s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4, &mu2, &ep);
      ep = -2;
      ivalue.val[2] = F77_FUNC(qli4,QLI4)(&s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4, &mu2, &ep);
#endif
    } else if (p2 != 0. && p4 != 0.) {
#ifdef USE_F2C
      DPKIND_COMPLEX rslt;
      ep = 0;
      F77_FUNC(qli4,QLI4)(&rslt, &s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1, &mu2, &ep);
      ivalue.val[0] = rslt;
      ep = -1;
      F77_FUNC(qli4,QLI4)(&rslt, &s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1, &mu2, &ep);
      ivalue.val[1] = rslt;
      ep = -2;
      F77_FUNC(qli4,QLI4)(&rslt, &s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1, &mu2, &ep);
      ivalue.val[2] = rslt;
#else
      ep = 0;
      ivalue.val[0] = F77_FUNC(qli4,QLI4)(&s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1, &mu2, &ep);
      ep = -1;
      ivalue.val[1] = F77_FUNC(qli4,QLI4)(&s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1, &mu2, &ep);
      ep = -2;
      ivalue.val[2] = F77_FUNC(qli4,QLI4)(&s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1, &mu2, &ep);
#endif
    } else { assert(0); }
  } else {
#ifdef USE_F2C
    DPKIND_COMPLEX rslt;
    ep = 0;
    F77_FUNC(qli4,QLI4)(&rslt, &p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4, &mu2, &ep);
    ivalue.val[0] = rslt;
    ep = -1;
    F77_FUNC(qli4,QLI4)(&rslt, &p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4, &mu2, &ep);
    ivalue.val[1] = rslt;
    ep = -2;
    F77_FUNC(qli4,QLI4)(&rslt, &p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4, &mu2, &ep);
    ivalue.val[2] = rslt;
#else
    ep = 0;
    ivalue.val[0] = F77_FUNC(qli4,QLI4)(&p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4, &mu2, &ep);
    ep = -1;
    ivalue.val[1] = F77_FUNC(qli4,QLI4)(&p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4, &mu2, &ep);
    ep = -2;
    ivalue.val[2] = F77_FUNC(qli4,QLI4)(&p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4, &mu2, &ep);
#endif
  }
  
  return ivalue;
}
#endif /* USE_QCDLOOP */

//
// OneLOop
//
#ifdef USE_ONELOOP
Initialize::Initialize() {
#ifndef NDEBUG
  std::cout << "PJFRY init" << std::endl;
#endif
  
  DPKIND_FLOAT mu = std::sqrt(ICache::getMu2());
  F77_FUNC_(avh_olo_mu_set,AVH_OLO_MU_SET)(&mu);
  
  DPKIND_FLOAT thrs = Minor5::getmeps();
  F77_FUNC_(avh_olo_onshell,AVH_OLO_ONSHELL)(&thrs);
#ifndef NDEBUG
  std::cout << "Set avh_olo_onshell threshold value to " << Minor5::getmeps() << std::endl;
#endif
}

Initialize::~Initialize() { /* NOOP */ }

ICache::Ival qlI1(const Kinem1& k) {
  ICache::Ival ivalue;
  
  DPKIND_FLOAT m1 = k.m1();
  DPKIND_COMPLEX rslt[3];
  F77_FUNC_(avh_olo_a0m,AVH_OLO_A0M)(rslt, &m1);
  
  ivalue.val[0] = rslt[0];
  ivalue.val[1] = rslt[1];
  ivalue.val[2] = rslt[2];
  return ivalue;
}

ICache::Ival qlI2(const Kinem2& k) {
  ICache::Ival ivalue;
  
  DPKIND_FLOAT p1 = k.p1();
  DPKIND_FLOAT m1 = k.m1();
  DPKIND_FLOAT m2 = k.m2();
  DPKIND_COMPLEX rslt[3];
  F77_FUNC_(avh_olo_b0m,AVH_OLO_B0M)(rslt, &p1, &m1, &m2);
  
  ivalue.val[0] = rslt[0];
  ivalue.val[1] = rslt[1];
  ivalue.val[2] = rslt[2];
  return ivalue;
}

ICache::Ival qlI3(const Kinem3& k) {
  ICache::Ival ivalue;
  
  DPKIND_FLOAT p1 = k.p1();
  DPKIND_FLOAT p2 = k.p2();
  DPKIND_FLOAT p3 = k.p3();
  DPKIND_FLOAT m1 = k.m1();
  DPKIND_FLOAT m2 = k.m2();
  DPKIND_FLOAT m3 = k.m3();
  DPKIND_COMPLEX rslt[3];
  F77_FUNC_(avh_olo_c0m,AVH_OLO_C0M)(rslt, &p1, &p2, &p3, &m1, &m2, &m3);
  
  ivalue.val[0] = rslt[0];
  ivalue.val[1] = rslt[1];
  ivalue.val[2] = rslt[2];
  return ivalue;
}

ICache::Ival qlI4(const Kinem4& k) {
  ICache::Ival ivalue;
  
  DPKIND_FLOAT p1  = k.p1();
  DPKIND_FLOAT p2  = k.p2();
  DPKIND_FLOAT p3  = k.p3();
  DPKIND_FLOAT p4  = k.p4();
  DPKIND_FLOAT s12 = k.s12();
  DPKIND_FLOAT s23 = k.s23();
  DPKIND_FLOAT m1  = k.m1();
  DPKIND_FLOAT m2  = k.m2();
  DPKIND_FLOAT m3  = k.m3();
  DPKIND_FLOAT m4  = k.m4();
  DPKIND_COMPLEX rslt[3];
  if (s12 == 0. || s23 == 0.) {
    if (p1 != 0. && p3 != 0.) {
      F77_FUNC_(avh_olo_d0m,AVH_OLO_D0M)(rslt, &s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4);
    } else if (p2 != 0. && p4 != 0.) {
      F77_FUNC_(avh_olo_d0m,AVH_OLO_D0M)(rslt, &s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1);
    } else { assert(0); }
  } else {
    F77_FUNC_(avh_olo_d0m,AVH_OLO_D0M)(rslt, &p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4);
  }
  
  ivalue.val[0] = rslt[0];
  ivalue.val[1] = rslt[1];
  ivalue.val[2] = rslt[2];
  return ivalue;
}
#endif /* USE_ONELOOP */
