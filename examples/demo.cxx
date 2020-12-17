//
// This C++ source code has been converted from the original
// Fortran "PJFry/examples/demo.f" source code.
//
// It can be used as a "Standalone Application" or as a "ROOT Script".
// Note: ROOT libraries are not needed for the "Standalone Application".
// For ROOT, see: https://root.cern.ch/
// "ROOT Script" usage example: root -b -l -q demo.cxx
//
// Last update: Sat Jan 18 20:55:00 UTC 2020
//
// Changes:
// 2020.01.18 - (initial version)
//

#include <iostream> // std::cout, std::scientific
#include <iomanip> // std::setprecision

#if defined(__CLING__)
// disable some specific clang diagnostics in ROOT
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
#endif /* defined(__CLING__) */

#include "pjfry.h"

// "ROOT Script" entry point
void demo() {
  double ps1, ps2, ps3, ps4, ps5;
  double ms1, ms2, ms3, ms4, ms5;
  double s12, s23, s34, s45, s15;
  
  // set the format used for floating-point values on output operations
  std::cout << std::scientific; // scientific notation
  std::cout << std::setprecision(17); // 18 decimal digits
  
  // set renormalization scale mu^2
  PJFry::SetMu2(36.);
  std::cout << "mu^2 = " << PJFry::GetMu2() << std::endl;
  
  std::cout << "-----------------------------------------------" << std::endl;
  
  std::cout << "2-mass point (m1 m1 -> m2 m2 0)" << std::endl;
  ps1 = 4;
  ps2 = 4;
  ps3 = 25;
  ps4 = 25;
  ps5 = 0;
  s12 =  0.10000000000000000e+05;
  s23 = -0.20042636597603382e+04;
  s34 =  0.10413130839415544e+04;
  s45 =  0.61341040415443358e+04;
  s15 = -0.33860497881996525e+04;
  ms1 = 0;
  ms2 = 4;
  ms3 = 0;
  ms4 = 25;
  ms5 = 0;
  
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
  
  std::cout << "E0(0) " << PJFry::E0v0(ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  std::cout << "E0(-1) " << PJFry::E0v0(ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 1) << std::endl;
  std::cout << "E0(-2) " << PJFry::E0v0(ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 2) << std::endl;
  std::cout << "E1(0) " << PJFry::E0v1(1, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  std::cout << "E2(0) " << PJFry::E0v1(2, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  
  std::cout << "E00(0) " << PJFry::E0v2(0, 0, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  std::cout << "E11(0) " << PJFry::E0v2(1, 1, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  std::cout << "E12(0) " << PJFry::E0v2(1, 2, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  
  std::cout << "E001(0) " << PJFry::E0v3(0, 0, 1, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  std::cout << "E112(0) " << PJFry::E0v3(1, 1, 2, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  std::cout << "E12(0)3 " << PJFry::E0v3(1, 2, 3, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  
  std::cout << "D111(0) " << PJFry::D0v3(1, 1, 1, s12, ps3, ps4, ps5, s45, s34, ms1, ms3, ms4, ms5, 0) << std::endl;
  std::cout << "D111(-1) " << PJFry::D0v3(1, 1, 1, s12, ps3, ps4, ps5, s45, s34, ms1, ms3, ms4, ms5, 1) << std::endl;
  std::cout << "D111(-2) " << PJFry::D0v3(1, 1, 1, s12, ps3, ps4, ps5, s45, s34, ms1, ms3, ms4, ms5, 2) << std::endl;
  std::cout << "D1233(0) " << PJFry::D0v4(1, 2, 3, 3, s12, ps3, ps4, ps5, s45, s34, ms1, ms3, ms4, ms5, 0) << std::endl;
  
  std::cout << "C122(0) " << PJFry::C0v3(1, 2, 2, s45, ps4, ps5, ms1, ms4, ms5, 0) << std::endl;
  
  std::cout << "-----------------------------------------------" << std::endl;
  
  std::cout << "0 mass, small ()_4 point" << std::endl;
  ps1 = 0;
  ps2 = 0;
  ps3 = 0;
  ps4 = 0;
  ps5 = 0;
  s12 = 25;
  s23 = -24.73938641761069;
  s34 =  24.52647702751698;
  s45 =  0.2605610444840831;
  s15 = -0.2145177937292013;
  ms1 = 0;
  ms2 = 0;
  ms3 = 0;
  ms4 = 0;
  ms5 = 0;
  
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
  
  std::cout << "E " << PJFry::E0v0(ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  std::cout << "E1 " << PJFry::E0v1(1, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  std::cout << "E12 " << PJFry::E0v2(1, 2, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  std::cout << "E123 " << PJFry::E0v3(1, 2, 3, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  std::cout << "E1234 " << PJFry::E0v4(1, 2, 3, 4, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
  std::cout << "E00234 " << PJFry::E0v5(0, 0, 2, 3, 4, ps1, ps2, ps3, ps4, ps5, s12, s23, s34, s45, s15, ms1, ms2, ms3, ms4, ms5, 0) << std::endl;
}

#if !defined(__CINT__) && !defined(__CLING__) && !defined(__ACLIC__)
// "Standalone Application" entry point
int main(int /*argc*/, char** /*argv*/) {
  demo();
  return 0;
}
#endif /* !defined(__CINT__) && !defined(__CLING__) && !defined(__ACLIC__) */

// end of file demo.cxx by Jacek M. Holeczek

