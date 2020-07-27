/*
 * integral.h - scalar integrals wrappers
 *
 * this file is part of PJFry library
 * Copyright 2011 Valery Yundin
 *
 * QCDLoop 2.x interface by Jacek M. Holeczek 2020.01.31
 */

#include "integral.h"
#include "minor.h"
#include "cache.h"

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
#include <iostream>

Initialize::Initialize() {
  std::cout << std::endl << "PJFRY init" << std::endl;
}

Initialize::~Initialize() {}

ICache::Ival qlI1(const Kinem1 &k) {
  static ql::TadPole<ql::complex, double, double> td;
  static std::vector<double> mI1(1);
  static std::vector<ql::complex> r(3);
  ICache::Ival ivalue;
  
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
  static ql::Bubble<ql::complex, double, double> bb;
  static std::vector<double> mI2(2);
  static std::vector<double> pI2(1);
  static std::vector<ql::complex> r(3);
  ICache::Ival ivalue;
  
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
  static ql::Triangle<ql::complex, double, double> tr;
  static std::vector<double> mI3(3);
  static std::vector<double> pI3(3);
  static std::vector<ql::complex> r(3);
  ICache::Ival ivalue;
  
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
  static ql::Box<ql::complex, double, double> bo;
  static std::vector<double> mI4(4);
  static std::vector<double> pI4(6);
  static std::vector<ql::complex> r(3);
  ICache::Ival ivalue;
  
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
// QCDLoop1 and QCDLoop 1.x / FF
//
#ifdef USE_QCDLOOP
Initialize::Initialize()
{
  printf("PJFRY init\n");
  const double dbl_min=std::numeric_limits<double>::min();

  F77_FUNC(qlinit,QLINIT)();

  if (qlprec.xalogm < dbl_min) {
    qlprec.xalogm=dbl_min;
    qlprec.xalog2=sqrt(dbl_min);
    #ifndef NDEBUG
    printf("Set xalogm to normalized value %e\n", qlprec.xalogm);
    printf("Set xalog2 to normalized value %e\n", qlprec.xalog2);
    #endif
  }
  if (qlprec.xclogm < dbl_min) {
    qlprec.xclogm=dbl_min;
    qlprec.xclog2=sqrt(dbl_min);
    #ifndef NDEBUG
    printf("Set xclogm to normalized value %e\n", qlprec.xclogm);
    printf("Set xclog2 to normalized value %e\n", qlprec.xclog2);
    #endif
  }
  if (qlflag.lwarn) {
    qlflag.lwarn=0;
    #ifndef NDEBUG
    printf("Disable FF warnings %d\n",qlflag.lwarn);
    #endif
  }
}

Initialize::~Initialize()
{
  F77_FUNC(ffexi,FFEXI)();
}

ICache::Ival qlI1(const Kinem1& k)
{
  int ep;
  ICache::Ival ivalue;

  double m1=k.m1();
  double mu2=ICache::getMu2();

#ifdef USE_F2C
  std::complex<double> rslt;
  ep=0;
  F77_FUNC(qli1,QLI1)(&rslt, &m1, &mu2, &ep);
  ivalue.val[0]=rslt;
  ep=-1;
  F77_FUNC(qli1,QLI1)(&rslt, &m1, &mu2, &ep);
  ivalue.val[1]=rslt;
  ep=-2;
  F77_FUNC(qli1,QLI1)(&rslt, &m1, &mu2, &ep);
  ivalue.val[2]=rslt;
#else
  ep=0;
  ivalue.val[0]=F77_FUNC(qli1,QLI1)(&m1, &mu2, &ep);
  ep=-1;
  ivalue.val[1]=F77_FUNC(qli1,QLI1)(&m1, &mu2, &ep);
  ep=-2;
  ivalue.val[2]=F77_FUNC(qli1,QLI1)(&m1, &mu2, &ep);
#endif

  return ivalue;
}

ICache::Ival qlI2(const Kinem2& k)
{
  int ep;
  ICache::Ival ivalue;

  double p1=k.p1();
  double m1=k.m1();
  double m2=k.m2();
  double mu2=ICache::getMu2();

#ifdef USE_F2C
  std::complex<double> rslt;
  ep=0;
  F77_FUNC(qli2,QLI2)(&rslt, &p1, &m1, &m2, &mu2, &ep);
  ivalue.val[0]=rslt;
  ep=-1;
  F77_FUNC(qli2,QLI2)(&rslt, &p1, &m1, &m2, &mu2, &ep);
  ivalue.val[1]=rslt;
  ep=-2;
  F77_FUNC(qli2,QLI2)(&rslt, &p1, &m1, &m2, &mu2, &ep);
  ivalue.val[2]=rslt;
#else
  ep=0;
  ivalue.val[0]=F77_FUNC(qli2,QLI2)(&p1, &m1, &m2, &mu2, &ep);
  ep=-1;
  ivalue.val[1]=F77_FUNC(qli2,QLI2)(&p1, &m1, &m2, &mu2, &ep);
  ep=-2;
  ivalue.val[2]=F77_FUNC(qli2,QLI2)(&p1, &m1, &m2, &mu2, &ep);
#endif

  return ivalue;
}

ICache::Ival qlI3(const Kinem3& k)
{
  int ep;
  ICache::Ival ivalue;

  double p1=k.p1();
  double p2=k.p2();
  double p3=k.p3();
  double m1=k.m1();
  double m2=k.m2();
  double m3=k.m3();
  double mu2=ICache::getMu2();

#ifdef USE_F2C
  std::complex<double> rslt;
  ep=0;
  F77_FUNC(qli3,QLI3)(&rslt, &p1, &p2, &p3, &m1, &m2, &m3, &mu2, &ep);
  ivalue.val[0]=rslt;
  ep=-1;
  F77_FUNC(qli3,QLI3)(&rslt, &p1, &p2, &p3, &m1, &m2, &m3, &mu2, &ep);
  ivalue.val[1]=rslt;
  ep=-2;
  F77_FUNC(qli3,QLI3)(&rslt, &p1, &p2, &p3, &m1, &m2, &m3, &mu2, &ep);
  ivalue.val[2]=rslt;
#else
  ep=0;
  ivalue.val[0]=F77_FUNC(qli3,QLI3)(&p1, &p2, &p3, &m1, &m2, &m3, &mu2, &ep);
  ep=-1;
  ivalue.val[1]=F77_FUNC(qli3,QLI3)(&p1, &p2, &p3, &m1, &m2, &m3, &mu2, &ep);
  ep=-2;
  ivalue.val[2]=F77_FUNC(qli3,QLI3)(&p1, &p2, &p3, &m1, &m2, &m3, &mu2, &ep);
#endif

  return ivalue;
}

ICache::Ival qlI4(const Kinem4& k)
{
  int ep;
  ICache::Ival ivalue;

  double p1=k.p1();
  double p2=k.p2();
  double p3=k.p3();
  double p4=k.p4();
  double s12=k.s12();
  double s23=k.s23();
  double m1=k.m1();
  double m2=k.m2();
  double m3=k.m3();
  double m4=k.m4();
  double mu2=ICache::getMu2();

  if (s12==0. || s23==0.) {
    if (p1!=0. && p3!=0.) {
#   ifdef USE_F2C
      std::complex<double> rslt;
      ep=0;
      F77_FUNC(qli4,QLI4)(&rslt, &s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4, &mu2, &ep);
      ivalue.val[0]=rslt;
      ep=-1;
      F77_FUNC(qli4,QLI4)(&rslt, &s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4, &mu2, &ep);
      ivalue.val[1]=rslt;
      ep=-2;
      F77_FUNC(qli4,QLI4)(&rslt, &s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4, &mu2, &ep);
      ivalue.val[2]=rslt;
#   else
      ep=0;
      ivalue.val[0]=F77_FUNC(qli4,QLI4)(&s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4, &mu2, &ep);
      ep=-1;
      ivalue.val[1]=F77_FUNC(qli4,QLI4)(&s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4, &mu2, &ep);
      ep=-2;
      ivalue.val[2]=F77_FUNC(qli4,QLI4)(&s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4, &mu2, &ep);
#   endif
    } else if (p2!=0. && p4!=0.) {
#   ifdef USE_F2C
      std::complex<double> rslt;
      ep=0;
      F77_FUNC(qli4,QLI4)(&rslt, &s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1, &mu2, &ep);
      ivalue.val[0]=rslt;
      ep=-1;
      F77_FUNC(qli4,QLI4)(&rslt, &s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1, &mu2, &ep);
      ivalue.val[1]=rslt;
      ep=-2;
      F77_FUNC(qli4,QLI4)(&rslt, &s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1, &mu2, &ep);
      ivalue.val[2]=rslt;
#   else
      ep=0;
      ivalue.val[0]=F77_FUNC(qli4,QLI4)(&s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1, &mu2, &ep);
      ep=-1;
      ivalue.val[1]=F77_FUNC(qli4,QLI4)(&s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1, &mu2, &ep);
      ep=-2;
      ivalue.val[2]=F77_FUNC(qli4,QLI4)(&s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1, &mu2, &ep);
#   endif
    } else { assert(0); }
  } else {
# ifdef USE_F2C
    std::complex<double> rslt;
    ep=0;
    F77_FUNC(qli4,QLI4)(&rslt, &p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4, &mu2, &ep);
    ivalue.val[0]=rslt;
    ep=-1;
    F77_FUNC(qli4,QLI4)(&rslt, &p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4, &mu2, &ep);
    ivalue.val[1]=rslt;
    ep=-2;
    F77_FUNC(qli4,QLI4)(&rslt, &p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4, &mu2, &ep);
    ivalue.val[2]=rslt;
# else
    ep=0;
    ivalue.val[0]=F77_FUNC(qli4,QLI4)(&p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4, &mu2, &ep);
    ep=-1;
    ivalue.val[1]=F77_FUNC(qli4,QLI4)(&p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4, &mu2, &ep);
    ep=-2;
    ivalue.val[2]=F77_FUNC(qli4,QLI4)(&p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4, &mu2, &ep);
# endif
  }

  return ivalue;
}
#endif /* USE_QCDLOOP */

//
// OneLOop
//
#ifdef USE_ONELOOP
Initialize::Initialize()
{
  printf("PJFRY init\n");

  double mu=sqrt(ICache::getMu2());
  F77_FUNC_(avh_olo_mu_set,AVH_OLO_MU_SET)(&mu);

  double thrs=Minor5::getmeps();
  F77_FUNC_(avh_olo_onshell,AVH_OLO_ONSHELL)(&thrs);
  #ifndef NDEBUG
  printf("Set avh_olo_onshell threshold value to %e\n", thrs);
  #endif
}

Initialize::~Initialize()
{
  /* NOOP */
}

ICache::Ival qlI1(const Kinem1& k)
{
  ICache::Ival ivalue;

  double m1=k.m1();
  std::complex<double> rslt[3];
  F77_FUNC_(avh_olo_a0m,AVH_OLO_A0M)(rslt, &m1);

  ivalue.val[0]=rslt[0];
  ivalue.val[1]=rslt[1];
  ivalue.val[2]=rslt[2];
  return ivalue;
}

ICache::Ival qlI2(const Kinem2& k)
{
  ICache::Ival ivalue;

  double p1=k.p1();
  double m1=k.m1();
  double m2=k.m2();
  std::complex<double> rslt[3];
  F77_FUNC_(avh_olo_b0m,AVH_OLO_B0M)(rslt, &p1, &m1, &m2);

  ivalue.val[0]=rslt[0];
  ivalue.val[1]=rslt[1];
  ivalue.val[2]=rslt[2];
  return ivalue;
}

ICache::Ival qlI3(const Kinem3& k)
{
  ICache::Ival ivalue;

  double p1=k.p1();
  double p2=k.p2();
  double p3=k.p3();
  double m1=k.m1();
  double m2=k.m2();
  double m3=k.m3();
  std::complex<double> rslt[3];
  F77_FUNC_(avh_olo_c0m,AVH_OLO_C0M)(rslt, &p1, &p2, &p3, &m1, &m2, &m3);

  ivalue.val[0]=rslt[0];
  ivalue.val[1]=rslt[1];
  ivalue.val[2]=rslt[2];
  return ivalue;
}

ICache::Ival qlI4(const Kinem4& k)
{
  ICache::Ival ivalue;

  double p1=k.p1();
  double p2=k.p2();
  double p3=k.p3();
  double p4=k.p4();
  double s12=k.s12();
  double s23=k.s23();
  double m1=k.m1();
  double m2=k.m2();
  double m3=k.m3();
  double m4=k.m4();
  std::complex<double> rslt[3];
  if (s12==0. || s23==0.) {
    if (p1!=0. && p3!=0.) {
      F77_FUNC_(avh_olo_d0m,AVH_OLO_D0M)(rslt, &s12, &p2, &s23, &p4, &p1, &p3, &m1, &m3, &m2, &m4);
    } else if (p2!=0. && p4!=0.) {
      F77_FUNC_(avh_olo_d0m,AVH_OLO_D0M)(rslt, &s23, &p3, &s12, &p1, &p2, &p4, &m2, &m4, &m3, &m1);
    } else { assert(0); }
  } else {
    F77_FUNC_(avh_olo_d0m,AVH_OLO_D0M)(rslt, &p1, &p2, &p3, &p4, &s12, &s23, &m1, &m2, &m3, &m4);
  }

  ivalue.val[0]=rslt[0];
  ivalue.val[1]=rslt[1];
  ivalue.val[2]=rslt[2];
  return ivalue;
}
#endif /* USE_ONELOOP */
