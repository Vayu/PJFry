//
// This C++ source code reproduces some results from:
//
// "A complete reduction of one-loop tensor 5- and 6-point integrals",
// Th. Diakonidis, J. Fleischer, J. Gluza, K. Kajda, T. Riemann, J.B. Tausk
// arXiv:0812.2134 [hep-ph] https://arxiv.org/abs/0812.2134
// https://journals.aps.org/prd/abstract/10.1103/PhysRevD.80.036003
//
// It can be used as a "Standalone Application" or as a "ROOT Script".
// Note: ROOT libraries are required in both cases.
// For ROOT, see: https://root.cern.ch/
// "ROOT Script" usage example: root -b -l -q arXiv_0812_2134.cxx
//
// Last update: 2020.07.15
//
// Changes:
// 2020.07.15 - (initial version)
//

#include <complex> // std::complex
#include <cmath> // std::abs
#include <iostream> // std::cout, std::defaultfloat, std::scientific
#include <iomanip> // std::setprecision

#include "Math/Vector4D.h" // ROOT::Math::PxPyPzEVector

#if defined(__CLING__)
// disable some specific clang diagnostics in ROOT
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
#endif /* defined(__CLING__) */

#include "pjfry.h"

// Kahan–Babushka-Neumaier summation algorithm
class KBNSum {
public:
  KBNSum() { s = c = std::complex<double>(0, 0); }
  void Add(const std::complex<double>& v) {
    std::complex<double> t = s + v;
    std::complex<double> cs = (s - t) + v; // if sum is bigger, low-order digits of v are lost
    std::complex<double> cv = (v - t) + s; // else low-order digits of sum are lost
    c += std::complex<double>( ((std::abs(s.real()) >= std::abs(v.real())) ? cs.real() : cv.real()),
                               ((std::abs(s.imag()) >= std::abs(v.imag())) ? cs.imag() : cv.imag()) );
    s = t;
  }
  KBNSum& operator+=(const std::complex<double>& v) { Add(v); return *this; }
  std::complex<double> Sum() const { return s + c; } // correction only applied once in the very end
  operator std::complex<double>() const { return Sum(); }
  std::complex<double> s; // a running sum
  std::complex<double> c; // a running compensation for lost low-order bits
};

// "ROOT Script" entry point
void arXiv_0812_2134(int on_shell = 1, int verbose = 0, int parameters = 0) {
  // set the format used for floating-point values on output operations
  // std::cout << std::defaultfloat; // defaultfloat notation
  // std::cout << std::setprecision(17); // std::defaultfloat 17 decimal digits
  std::cout << std::scientific; // scientific notation
  std::cout << std::setprecision(16); // std::scientific 17 decimal digits
  
  {
    // "Table 5.1" a.k.a. "TABLE I"
    // notes:
    //   1. it describes a hexagon diagram for g g -> t t_bar q q_bar, where
    //      the light quarks q, q_bar are considered to be exactly massless,
    //   2. the internal lines with masses m1, ..., m6 are set to some quite
    //      arbitrary values, for which there is no physical motivation at all
    //      (a test case with as many different masses as possible)
    // Bas Tausk, private communication (note: 2 * 0.217745543033608E+03 = 0.435491086067216E+03)
    ROOT::Math::PxPyPzEVector p1( 0.0,                    0.0,                    0.217745543033608E+03,  0.217745543033608E+03); // g
    ROOT::Math::PxPyPzEVector p2( 0.0,                    0.0,                   -0.217745543033608E+03,  0.217745543033608E+03); // g
    ROOT::Math::PxPyPzEVector p3(-0.475795122438478E+02,  0.421268225217969E+02,  0.840971806502753E+02, -0.203694145606766E+03); // t
    ROOT::Math::PxPyPzEVector p4( 0.552159614718998E+02, -0.466920343021516E+02, -0.900100867242416E+02, -0.209072365894312E+03); // t_bar
    ROOT::Math::PxPyPzEVector p5( 0.530631952055453E+01,  0.296982674372730E+01, -0.314568707927384E+01, -0.684633076491974E+01); // q
    ROOT::Math::PxPyPzEVector p6(-0.129427687486066E+02,  0.159538503662741E+01,  0.905859315324015E+01, -0.158782438012184E+02); // q_bar
    if (parameters == 1) p6 = -(p1 + p2 + p3 + p4 + p5);
    if ((parameters == 2) || (parameters == 3)) { // low precision values "published" in papers
      p1 = ROOT::Math::PxPyPzEVector( 0.0,             0.0,             0.21774554E+03,  0.21774554E+03); // g
      p2 = ROOT::Math::PxPyPzEVector( 0.0,             0.0,            -0.21774554E+03,  0.21774554E+03); // g
      p3 = ROOT::Math::PxPyPzEVector(-0.47579512E+02,  0.42126823E+02,  0.84097181E+02, -0.20369415E+03); // t
      p4 = ROOT::Math::PxPyPzEVector( 0.55215961E+02, -0.46692034E+02, -0.90010087E+02, -0.20907237E+03); // t_bar
      p5 = ROOT::Math::PxPyPzEVector( 0.53063195E+01,  0.29698267E+01, -0.31456871E+01, -0.68463308E+01); // q
      p6 = ROOT::Math::PxPyPzEVector(-0.12942769E+02,  0.15953850E+01,  0.90585932E+01, -0.15878244E+02); // q_bar
      if (parameters == 3) p6 = -(p1 + p2 + p3 + p4 + p5);
    }
    // ROOT::Math::PxPyPzEVector p0 = p1 + p2 + p3 + p4 + p5 + p6;
    double m1(110.0), m2(120.0), m3(130.0), m4(140.0), m5(150.0), m6(160.0);
    
    // note: for five-point functions, we shrink line 2, drop m2 and, in
    //       order to retain momentum conservation, fix p1 + p2 -> p1
    double ps1 = (p1 + p2).M2();
    double ps2 = p3.M2();
    double ps3 = p4.M2();
    double ps4 = p5.M2();
    double ps5 = p6.M2();
    double s12 = ((p1 + p2) + p3).M2();
    double s23 = (p3 + p4).M2();
    double s34 = (p4 + p5).M2();
    double s45 = (p5 + p6).M2();
    double s15 = ((p1 + p2) + p6).M2();
    double ms1 = m1 * m1; // note: m2 is not used at all
    double ms2 = m3 * m3;
    double ms3 = m4 * m4;
    double ms4 = m5 * m5;
    double ms5 = m6 * m6;
    
    if (on_shell > 0) {
      // put the momenta exactly on shell
      ps2 = ps3 = 174.30 * 174.30; // t, t_bar
      ps4 = ps5 = 0; // q, q_bar (massless)
    }
    
    // set the regularization scale squared
    double mu2 = 1.0;
    PJFry::SetMu2(mu2); // note: no massless internal particles so it does NOT matter here
    // std::cout << std::endl << "regularization scale squared = " << PJFry::GetMu2() << std::endl;
    
    if (verbose > 0) {
      std::cout << std::endl;
      std::cout << "... Table 5.1 (a.k.a. TABLE I) ..." << std::endl;
      std::cout << std::endl;
      std::cout << "mu2 = " << mu2 << std::endl;
      std::cout << "ps1 = " << ps1 << std::endl;
      std::cout << "ps2 = " << ps2 << std::endl;
      std::cout << "ps3 = " << ps3 << std::endl;
      std::cout << "ps4 = " << ps4 << std::endl;
      std::cout << "ps5 = " << ps5 << std::endl;
      std::cout << "s12 = " << s12 << std::endl;
      std::cout << "s23 = " << s23 << std::endl;
      std::cout << "s34 = " << s34 << std::endl;
      std::cout << "s45 = " << s45 << std::endl;
      std::cout << "s15 = " << s15 << std::endl;
      std::cout << "ms1 = " << ms1 << std::endl;
      std::cout << "ms2 = " << ms2 << std::endl;
      std::cout << "ms3 = " << ms3 << std::endl;
      std::cout << "ms4 = " << ms4 << std::endl;
      std::cout << "ms5 = " << ms5 << std::endl;
    }
    
    // "Table 5.7" a.k.a. "TABLE VII"
    // note: no massless internal particles so all "E*(-1)" and "E*(-2)" coefficients should be 0
    std::cout << std::endl;
    std::cout << "... Table 5.7 (a.k.a. TABLE VII) ..." << std::endl;
    std::cout << std::endl;
#if 1 /* 0 or 1 */
    std::cout << "E_0 ( 0)    = " << PJFry::E0v0(ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_0 (-1)    = " << PJFry::E0v0(ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_0 (-2)    = " << PJFry::E0v0(ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
#else /* 0 or 1 */
    std::cout << "E_0 ( 0)    = " << PJFry::E0v1(0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_0 (-1)    = " << PJFry::E0v1(0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_0 (-2)    = " << PJFry::E0v1(0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
#endif /* 0 or 1 */
    std::cout << "E_1 ( 0)    = " << PJFry::E0v1(1, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_1 (-1)    = " << PJFry::E0v1(1, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_1 (-2)    = " << PJFry::E0v1(1, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
    std::cout << "E_12 ( 0)   = " << PJFry::E0v2(1, 2, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_12 (-1)   = " << PJFry::E0v2(1, 2, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_12 (-2)   = " << PJFry::E0v2(1, 2, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
    std::cout << "E_00 ( 0)   = " << PJFry::E0v2(0, 0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_00 (-1)   = " << PJFry::E0v2(0, 0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_00 (-2)   = " << PJFry::E0v2(0, 0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
    std::cout << "E_123 ( 0)  = " << PJFry::E0v3(1, 2, 3, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_123 (-1)  = " << PJFry::E0v3(1, 2, 3, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_123 (-2)  = " << PJFry::E0v3(1, 2, 3, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
    std::cout << "E_001 ( 0)  = " << PJFry::E0v3(0, 0, 1, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_001 (-1)  = " << PJFry::E0v3(0, 0, 1, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_001 (-2)  = " << PJFry::E0v3(0, 0, 1, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
  }
  
  {
    // "Table 5.2" a.k.a. "TABLE II"
    // note: all internal particles are massless
    // "golem95/demos/demo_6point.f90" and "golem95/demos/demo_cmplx_masses.f90"
    ROOT::Math::PxPyPzEVector p1( 0.0,                  0.0,                  0.5,                  0.5);
    ROOT::Math::PxPyPzEVector p2( 0.0,                  0.0,                 -0.5,                  0.5);
    ROOT::Math::PxPyPzEVector p3(-0.12741179719516801, -0.08262476614744381, -0.11713105190921771, -0.19178191094778038);
    ROOT::Math::PxPyPzEVector p4( 0.06648281097623857,  0.3189378514746887,   0.08471424069583446, -0.33662712284553753);
    ROOT::Math::PxPyPzEVector p5( 0.20363139428835617, -0.044157623555325,   -0.0571065672034082,  -0.21604814388379073);
    ROOT::Math::PxPyPzEVector p6(-0.1427024080694266,  -0.19215546177191994,  0.08952337841679145, -0.2555428223228916);
    if (parameters == 1) p6 = -(p1 + p2 + p3 + p4 + p5);
    if ((parameters == 2) || (parameters == 3)) { // low precision values "published" in papers
      p1 = ROOT::Math::PxPyPzEVector( 0.0,         0.0,         0.5,         0.5);
      p2 = ROOT::Math::PxPyPzEVector( 0.0,         0.0,        -0.5,         0.5);
      p3 = ROOT::Math::PxPyPzEVector(-0.12741180, -0.08262477, -0.11713105, -0.19178191);
      p4 = ROOT::Math::PxPyPzEVector( 0.06648281,  0.31893785,  0.08471424, -0.33662712);
      p5 = ROOT::Math::PxPyPzEVector( 0.20363139, -0.04415762, -0.05710657, -0.21604814);
      p6 = -(p1 + p2 + p3 + p4 + p5);
    }
    // ROOT::Math::PxPyPzEVector p0 = p1 + p2 + p3 + p4 + p5 + p6;
    double m1(0.0), m2(0.0), m3(0.0), m4(0.0), m5(0.0), m6(0.0);
    
    // note: for five-point functions, we shrink line 2, drop m2 and, in
    //       order to retain momentum conservation, fix p1 + p2 -> p1
    double ps1 = (p1 + p2).M2();
    double ps2 = p3.M2();
    double ps3 = p4.M2();
    double ps4 = p5.M2();
    double ps5 = p6.M2();
    double s12 = ((p1 + p2) + p3).M2();
    double s23 = (p3 + p4).M2();
    double s34 = (p4 + p5).M2();
    double s45 = (p5 + p6).M2();
    double s15 = ((p1 + p2) + p6).M2();
    double ms1 = m1 * m1; // note: m2 is not used at all
    double ms2 = m3 * m3;
    double ms3 = m4 * m4;
    double ms4 = m5 * m5;
    double ms5 = m6 * m6;
    
    if (on_shell > 0) {
      // put the momenta exactly on shell
      ps2 = ps3 = ps4 = ps5 = 0; // all are massless
    }
    
    // set the regularization scale squared
    double mu2 = 1.0;
    PJFry::SetMu2(mu2); // note: all internal particles are massless so it DOES matter here
    // std::cout << std::endl << "regularization scale squared = " << PJFry::GetMu2() << std::endl;
    
    if (verbose > 0) {
      std::cout << std::endl;
      std::cout << "... Table 5.2 (a.k.a. TABLE II) ..." << std::endl;
      std::cout << std::endl;
      std::cout << "mu2 = " << mu2 << std::endl;
      std::cout << "ps1 = " << ps1 << std::endl;
      std::cout << "ps2 = " << ps2 << std::endl;
      std::cout << "ps3 = " << ps3 << std::endl;
      std::cout << "ps4 = " << ps4 << std::endl;
      std::cout << "ps5 = " << ps5 << std::endl;
      std::cout << "s12 = " << s12 << std::endl;
      std::cout << "s23 = " << s23 << std::endl;
      std::cout << "s34 = " << s34 << std::endl;
      std::cout << "s45 = " << s45 << std::endl;
      std::cout << "s15 = " << s15 << std::endl;
      std::cout << "ms1 = " << ms1 << std::endl;
      std::cout << "ms2 = " << ms2 << std::endl;
      std::cout << "ms3 = " << ms3 << std::endl;
      std::cout << "ms4 = " << ms4 << std::endl;
      std::cout << "ms5 = " << ms5 << std::endl;
    }
    
    // "Table 5.9" a.k.a. "TABLE IX"
    std::cout << std::endl;
    std::cout << "... Table 5.9 (a.k.a. TABLE IX) ..." << std::endl;
    std::cout << std::endl;
#if 1 /* 0 or 1 */
    std::cout << "E_0 ( 0)    = " << PJFry::E0v0(ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_0 (-1)    = " << PJFry::E0v0(ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_0 (-2)    = " << PJFry::E0v0(ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
#else /* 0 or 1 */
    std::cout << "E_0 ( 0)    = " << PJFry::E0v1(0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_0 (-1)    = " << PJFry::E0v1(0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_0 (-2)    = " << PJFry::E0v1(0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
#endif /* 0 or 1 */
    std::cout << "E_2 ( 0)    = " << PJFry::E0v1(2, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_2 (-1)    = " << PJFry::E0v1(2, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_2 (-2)    = " << PJFry::E0v1(2, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
    std::cout << "E_33 ( 0)   = " << PJFry::E0v2(3, 3, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_33 (-1)   = " << PJFry::E0v2(3, 3, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_33 (-2)   = " << PJFry::E0v2(3, 3, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
    std::cout << "E_00 ( 0)   = " << PJFry::E0v2(0, 0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_00 (-1)   = " << PJFry::E0v2(0, 0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_00 (-2)   = " << PJFry::E0v2(0, 0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
    std::cout << "E_444 ( 0)  = " << PJFry::E0v3(4, 4, 4, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_444 (-1)  = " << PJFry::E0v3(4, 4, 4, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_444 (-2)  = " << PJFry::E0v3(4, 4, 4, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
    std::cout << "E_001 ( 0)  = " << PJFry::E0v3(0, 0, 1, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
    std::cout << "E_001 (-1)  = " << PJFry::E0v3(0, 0, 1, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
    std::cout << "E_001 (-2)  = " << PJFry::E0v3(0, 0, 1, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
    std::cout << std::endl;
  }
  
  {
    // "Table 5.3" a.k.a. "TABLE III"
    // "Table 5.1" a.k.a. "TABLE I" scaled down by a factor 100
    // notes:
    //   1. it describes a hexagon diagram for g g -> t t_bar q q_bar, where
    //      the light quarks q, q_bar are considered to be exactly massless,
    //   2. the internal lines with masses m1, ..., m6 are g, g, g, t, g, q
    // Bas Tausk, private communication (note: 2 * 0.217745543033608E+03 = 0.435491086067216E+03)
    ROOT::Math::PxPyPzEVector p1( 0.0,                    0.0,                    0.217745543033608E+03,  0.217745543033608E+03); // g
    ROOT::Math::PxPyPzEVector p2( 0.0,                    0.0,                   -0.217745543033608E+03,  0.217745543033608E+03); // g
    ROOT::Math::PxPyPzEVector p3(-0.475795122438478E+02,  0.421268225217969E+02,  0.840971806502753E+02, -0.203694145606766E+03); // t
    ROOT::Math::PxPyPzEVector p4( 0.552159614718998E+02, -0.466920343021516E+02, -0.900100867242416E+02, -0.209072365894312E+03); // t_bar
    ROOT::Math::PxPyPzEVector p5( 0.530631952055453E+01,  0.296982674372730E+01, -0.314568707927384E+01, -0.684633076491974E+01); // q
    ROOT::Math::PxPyPzEVector p6(-0.129427687486066E+02,  0.159538503662741E+01,  0.905859315324015E+01, -0.158782438012184E+02); // q_bar
    if (parameters == 1) p6 = -(p1 + p2 + p3 + p4 + p5);
    if ((parameters == 2) || (parameters == 3)) { // low precision values "published" in papers
      p1 = ROOT::Math::PxPyPzEVector( 0.0,             0.0,             0.21774554E+03,  0.21774554E+03); // g
      p2 = ROOT::Math::PxPyPzEVector( 0.0,             0.0,            -0.21774554E+03,  0.21774554E+03); // g
      p3 = ROOT::Math::PxPyPzEVector(-0.47579512E+02,  0.42126823E+02,  0.84097181E+02, -0.20369415E+03); // t
      p4 = ROOT::Math::PxPyPzEVector( 0.55215961E+02, -0.46692034E+02, -0.90010087E+02, -0.20907237E+03); // t_bar
      p5 = ROOT::Math::PxPyPzEVector( 0.53063195E+01,  0.29698267E+01, -0.31456871E+01, -0.68463308E+01); // q
      p6 = ROOT::Math::PxPyPzEVector(-0.12942769E+02,  0.15953850E+01,  0.90585932E+01, -0.15878244E+02); // q_bar
      if (parameters == 3) p6 = -(p1 + p2 + p3 + p4 + p5);
    }
    p1 /= 100; p2 /= 100; p3 /= 100; p4 /= 100; p5 /= 100; p6 /= 100;
    // ROOT::Math::PxPyPzEVector p0 = p1 + p2 + p3 + p4 + p5 + p6;
    double m1(0.0), m2(0.0), m3(0.0), m4(1.7430), m5(0.0), m6(0.0);
    
    // note: for five-point functions, we shrink line 2, drop m2 and, in
    //       order to retain momentum conservation, fix p1 + p2 -> p1
    double ps1 = (p1 + p2).M2();
    double ps2 = p3.M2();
    double ps3 = p4.M2();
    double ps4 = p5.M2();
    double ps5 = p6.M2();
    double s12 = ((p1 + p2) + p3).M2();
    double s23 = (p3 + p4).M2();
    double s34 = (p4 + p5).M2();
    double s45 = (p5 + p6).M2();
    double s15 = ((p1 + p2) + p6).M2();
    double ms1 = m1 * m1; // note: m2 is not used at all
    double ms2 = m3 * m3;
    double ms3 = m4 * m4;
    double ms4 = m5 * m5;
    double ms5 = m6 * m6;
    
    if (on_shell > 0) {
      // put the momenta exactly on shell
      ps2 = ps3 = 1.7430 * 1.7430; // t, t_bar (= ms3)
      ps4 = ps5 = 0; // q, q_bar (massless)
    }
    
#if 1 /* 0 or 1 */
    ROOT::Math::PxPyPzEVector q[5] = {
      ROOT::Math::PxPyPzEVector(0, 0, 0, 0), // q0 = q4 + p6 = 0
      p1 + p2, // q1 = q0 + p1 + p2
      (p1 + p2) + p3, // q2 = q1 + p3
      -(p5 + p6), // q3 = q2 + p4 = -(p5 + p6)
      -p6 // q4 = q3 + p5 = -p6
    };
#else /* 0 or 1 */
    ROOT::Math::PxPyPzEVector q[5] = {
      ROOT::Math::PxPyPzEVector(0, 0, 0, 0), // -q0 = -q4 + p6 = 0
      -(p1 + p2), // -q1 = -q0 + p1 + p2
      -((p1 + p2) + p3), // -q2 = -q1 + p3
      p5 + p6, // -q3 = -q2 + p4 = -(p5 + p6)
      p6 // -q4 = -q3 + p5 = -p6
    };
#endif /* 0 or 1 */
    
    // set the regularization scale squared
    double mu2 = 1.0;
    PJFry::SetMu2(mu2); // note: some internal particles are massless so it DOES matter here
    // std::cout << std::endl << "regularization scale squared = " << PJFry::GetMu2() << std::endl;
    
    if (verbose > 0) {
      std::cout << std::endl;
      std::cout << "... Table 5.3 (a.k.a. TABLE III) ..." << std::endl;
      std::cout << std::endl;
      std::cout << "mu2 = " << mu2 << std::endl;
      std::cout << "ps1 = " << ps1 << std::endl;
      std::cout << "ps2 = " << ps2 << std::endl;
      std::cout << "ps3 = " << ps3 << std::endl;
      std::cout << "ps4 = " << ps4 << std::endl;
      std::cout << "ps5 = " << ps5 << std::endl;
      std::cout << "s12 = " << s12 << std::endl;
      std::cout << "s23 = " << s23 << std::endl;
      std::cout << "s34 = " << s34 << std::endl;
      std::cout << "s45 = " << s45 << std::endl;
      std::cout << "s15 = " << s15 << std::endl;
      std::cout << "ms1 = " << ms1 << std::endl;
      std::cout << "ms2 = " << ms2 << std::endl;
      std::cout << "ms3 = " << ms3 << std::endl;
      std::cout << "ms4 = " << ms4 << std::endl;
      std::cout << "ms5 = " << ms5 << std::endl;
    }
    
    if (verbose > 1) {
      for (int i = 0; i < 5; i++) std::cout << "q" << i << "  = " << q[i] << std::endl;
    }
    
    // "Table 5.10" a.k.a. "TABLE X"
    std::cout << std::endl;
    std::cout << "... Table 5.10 (a.k.a. TABLE X) ..." << std::endl;
    std::cout << std::endl;
    
    std::complex<double> E0[3]; // note: pj_complex is just a std::complex<double>
#if 1 /* 0 or 1 */
    for (int n = 0; n < 3; n++) {
      E0[n] = PJFry::E0v0(ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, n);
    }
#else /* 0 or 1 */
    for (int n = 0; n < 3; n++) {
      E0[n] = PJFry::E0v1(0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, n);
    }
#endif /* 0 or 1 */
    std::cout << "E_0 ( 0)    = " << E0[0] << std::endl;
    std::cout << "E_0 (-1)    = " << E0[1] << std::endl;
    std::cout << "E_0 (-2)    = " << E0[2] << std::endl;
    std::cout << std::endl;
    
    std::complex<double> E3[3]; // note: pj_complex is just a std::complex<double>
    for (int n = 0; n < 3; n++) {
      KBNSum sum;
      for (int i = 1; i <= 4; i++) {
        sum +=
          (q[i]).Z() * // 3
          PJFry::E0v1(i, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, n);
        if (verbose > 2) std::cout << sum.s << "   \t" << sum.c << std::endl;
      }
      E3[n] = sum;
      if (verbose > 2) std::cout << E3[n] << std::endl;
    }
    std::cout << "E^3 ( 0)    = " << E3[0] << std::endl;
    std::cout << "E^3 (-1)    = " << E3[1] << std::endl;
    std::cout << "E^3 (-2)    = " << E3[2] << std::endl;
    std::cout << std::endl;
    
    std::complex<double> E23[3]; // note: pj_complex is just a std::complex<double>
    for (int n = 0; n < 3; n++) {
      KBNSum sum;
      for (int i = 1; i <= 4; i++) {
        for (int j = 1; j <= 4; j++) {
          sum +=
            (q[i]).Y() * // 2.
            (q[j]).Z() * // .3
            PJFry::E0v2(i, j, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, n);
          if (verbose > 2) std::cout << sum.s << "   \t" << sum.c << std::endl;
        }
      }
      E23[n] = sum;
      if (verbose > 2) std::cout << E23[n] << std::endl;
    }
    std::cout << "E^23 ( 0)   = " << E23[0] << std::endl;
    std::cout << "E^23 (-1)   = " << E23[1] << std::endl;
    std::cout << "E^23 (-2)   = " << E23[2] << std::endl;
    std::cout << std::endl;
    
    std::complex<double> E012[3]; // note: pj_complex is just a std::complex<double>
    for (int n = 0; n < 3; n++) {
      KBNSum sum;
      for (int i = 1; i <= 4; i++) {
        for (int j = 1; j <= 4; j++) {
          for (int k = 1; k <= 4; k++) {
            sum +=
              (q[i]).E() * // 0..
              (q[j]).X() * // .1.
              (q[k]).Y() * // ..2
              PJFry::E0v3(i, j, k, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, n);
            if (verbose > 2) std::cout << sum.s << "   \t" << sum.c << std::endl;
          }
        }
      }
      E012[n] = sum;
      if (verbose > 2) std::cout << E012[n] << std::endl;
    }
    std::cout << "E^012 ( 0)  = " << E012[0] << std::endl;
    std::cout << "E^012 (-1)  = " << E012[1] << std::endl;
    std::cout << "E^012 (-2)  = " << E012[2] << std::endl;
    std::cout << std::endl;
    
    std::complex<double> E2130[3]; // note: pj_complex is just a std::complex<double>
    for (int n = 0; n < 3; n++) {
      KBNSum sum;
      for (int i = 1; i <= 4; i++) {
        for (int j = 1; j <= 4; j++) {
          for (int k = 1; k <= 4; k++) {
            for (int l = 1; l <= 4; l++) {
              sum +=
                (q[i]).Y() * // 2...
                (q[j]).X() * // .1..
                (q[k]).Z() * // ..3.
                (q[l]).E() * // ...0
                PJFry::E0v4(i, j, k, l, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, n);
              if (verbose > 2) std::cout << sum.s << "   \t" << sum.c << std::endl;
            }
          }
        }
      }
      E2130[n] = sum;
      if (verbose > 2) std::cout << E2130[n] << std::endl;
    }
    std::cout << "E^2130 ( 0) = " << E2130[0] << std::endl;
    std::cout << "E^2130 (-1) = " << E2130[1] << std::endl;
    std::cout << "E^2130 (-2) = " << E2130[2] << std::endl;
    std::cout << std::endl;
  }
}

#if !defined(__CINT__) && !defined(__CLING__) && !defined(__ACLIC__)
// "Standalone Application" entry point
int main(int /*argc*/, char** /*argv*/) {
  arXiv_0812_2134();
  return 0;
}
#endif /* !defined(__CINT__) && !defined(__CLING__) && !defined(__ACLIC__) */

// end of file arXiv_0812_2134.cxx by Jacek M. Holeczek

