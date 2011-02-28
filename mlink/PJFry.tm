/* To launch this program from within Mathematica use:
*   In[1]:= link = Install["lhapdf"]
*
* Or, launch this program from a shell and establish a
* peer-to-peer connection.  When given the prompt Create Link:
* type a port name. ( On Unix platforms, a port name is a
* number less than 65536.  On Mac or Windows platforms,
* it's an arbitrary word.)
* Then, from within Mathematica use:
*   In[1]:= link = Install["portname", LinkMode->Connect]
*/

#include <iostream>
#include <mathlink.h>
#include "../src/pjfry.h"

// ====================================================================

:Evaluate:  Print["PJFry MathLink"];
:Evaluate:  Print["Type Names[\"PJFry`*\"] to show exported names"];

:Evaluate:  BeginPackage["PJFry`"]

:Evaluate:  GetMu2::usage = "GetMu2[] returns the renormalization scale squared"
:Evaluate:  SetMu2::usage = "SetMu2[mu^2] sets the renormalization scale squared"

:Evaluate:  A0v0::usage  = "A0v0[m1,ep=0]"

:Evaluate:  B0v0::usage  = "B0v0[p1,m1,m2,ep=0]"
:Evaluate:  B0v1::usage  = "B0v1[i,p1,m1,m2,ep=0]"
:Evaluate:  B0v2::usage  = "B0v2[i,j,p1,m1,m2,ep=0]"

:Evaluate:  C0v0::usage  = "C0v0[p1,p2,p3,m1,m2,m3,ep=0]"
:Evaluate:  C0v1::usage  = "C0v1[i,p1,p2,p3,m1,m2,m3,ep=0]"
:Evaluate:  C0v2::usage  = "C0v2[i,j,p1,p2,p3,m1,m2,m3,ep=0]"
:Evaluate:  C0v3::usage  = "C0v3[i,j,k,p1,p2,p3,m1,m2,m3,ep=0]"

:Evaluate:  D0v0::usage  = "D0v0[p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,ep=0]"
:Evaluate:  D0v1::usage  = "D0v1[i,p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,ep=0]"
:Evaluate:  D0v2::usage  = "D0v2[i,j,p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,ep=0]"
:Evaluate:  D0v3::usage  = "D0v3[i,j,k,p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,ep=0]"
:Evaluate:  D0v4::usage  = "D0v4[i,j,k,l,p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,ep=0]"

:Evaluate:  E0v0::usage  = "E0v0[p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,m1,m2,m3,m4,m5,ep=0]"
:Evaluate:  E0v1::usage  = "E0v1[i,p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,m1,m2,m3,m4,m5,ep=0]"
:Evaluate:  E0v2::usage  = "E0v2[i,j,p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,m1,m2,m3,m4,m5,ep=0]"
:Evaluate:  E0v3::usage  = "E0v3[i,j,k,p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,m1,m2,m3,m4,m5,ep=0]"
:Evaluate:  E0v4::usage  = "E0v4[i,j,k,l,p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,m1,m2,m3,m4,m5,ep=0]"
:Evaluate:  E0v5::usage  = "E0v5[i,j,k,l,m,p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,m1,m2,m3,m4,m5,ep=0]"

:Evaluate:  Begin["`Private`"]

:Evaluate:  SetAttributes[GetMu2, NumericFunction]
:Evaluate:  SetAttributes[SetMu2, NumericFunction]

:Evaluate:  SetAttributes[A0v0, NumericFunction]

:Evaluate:  SetAttributes[B0v0, NumericFunction]
:Evaluate:  SetAttributes[B0v1, NumericFunction]
:Evaluate:  SetAttributes[B0v2, NumericFunction]

:Evaluate:  SetAttributes[C0v0, NumericFunction]
:Evaluate:  SetAttributes[C0v1, NumericFunction]
:Evaluate:  SetAttributes[C0v2, NumericFunction]
:Evaluate:  SetAttributes[C0v3, NumericFunction]

:Evaluate:  SetAttributes[D0v0, NumericFunction]
:Evaluate:  SetAttributes[D0v1, NumericFunction]
:Evaluate:  SetAttributes[D0v2, NumericFunction]
:Evaluate:  SetAttributes[D0v3, NumericFunction]
:Evaluate:  SetAttributes[D0v4, NumericFunction]

:Evaluate:  SetAttributes[E0v0, NumericFunction]
:Evaluate:  SetAttributes[E0v1, NumericFunction]
:Evaluate:  SetAttributes[E0v2, NumericFunction]
:Evaluate:  SetAttributes[E0v3, NumericFunction]
:Evaluate:  SetAttributes[E0v4, NumericFunction]
:Evaluate:  SetAttributes[E0v5, NumericFunction]

// ==================== My wrappers =========================

:Begin:
:Function: fGetMu2
:Pattern: GetMu2[]
:Arguments: {}
:ArgumentTypes: {}
:ReturnType: Real64
:End:

:Begin:
:Function: fSetMu2
:Pattern: SetMu2[mu2_?NumericQ]
:Arguments: {N[mu2]}
:ArgumentTypes: {Real64}
:ReturnType: Real64
:End:

:Begin:
:Function: fA0v0
:Pattern: A0v0[m1_?NumericQ,ep_Integer:0]
:Arguments: {N[m1],ep}
:ArgumentTypes: {Real64,Integer}
:ReturnType: Manual
:End:

// ==================== BUBBLES =========================

:Begin:
:Function: fB0v0
:Pattern: B0v0[p1_?NumericQ,m1_?NumericQ,m2_?NumericQ,ep_Integer:0]
:Arguments: {N[p1],N[m1],N[m2],ep}
:ArgumentTypes: {Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fB0v1
:Pattern: B0v1[i_Integer,p1_?NumericQ,m1_?NumericQ,m2_?NumericQ,ep_Integer:0]
:Arguments: {i,N[p1],N[m1],N[m2],ep}
:ArgumentTypes: {Integer,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fB0v2
:Pattern: B0v2[i_Integer,j_Integer,p1_?NumericQ,m1_?NumericQ,m2_?NumericQ,ep_Integer:0]
:Arguments: {i,j,N[p1],N[m1],N[m2],ep}
:ArgumentTypes: {Integer,Integer,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

// ==================== TRIANGLES =========================

:Begin:
:Function: fC0v0
:Pattern: C0v0[p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,ep_Integer:0]
:Arguments: {N[p1],N[p2],N[p3],N[m1],N[m2],N[m3],ep}
:ArgumentTypes: {Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fC0v1
:Pattern: C0v1[i_Integer,p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,ep_Integer:0]
:Arguments: {i,N[p1],N[p2],N[p3],N[m1],N[m2],N[m3],ep}
:ArgumentTypes: {Integer,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fC0v2
:Pattern: C0v2[i_Integer,j_Integer,p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,ep_Integer:0]
:Arguments: {i,j,N[p1],N[p2],N[p3],N[m1],N[m2],N[m3],ep}
:ArgumentTypes: {Integer,Integer,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fC0v3
:Pattern: C0v3[i_Integer,j_Integer,k_Integer,p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,ep_Integer:0]
:Arguments: {i,j,k,N[p1],N[p2],N[p3],N[m1],N[m2],N[m3],ep}
:ArgumentTypes: {Integer,Integer,Integer,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

// ==================== BOXES =========================

:Begin:
:Function: fD0v0
:Pattern: D0v0[p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,p4_?NumericQ,s12_?NumericQ,s23_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,m4_?NumericQ,ep_Integer:0]
:Arguments: {N[p1],N[p2],N[p3],N[p4],N[s12],N[s23],N[m1],N[m2],N[m3],N[m4],ep}
:ArgumentTypes: {Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fD0v1
:Pattern: D0v1[i_Integer,p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,p4_?NumericQ,s12_?NumericQ,s23_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,m4_?NumericQ,ep_Integer:0]
:Arguments: {i,N[p1],N[p2],N[p3],N[p4],N[s12],N[s23],N[m1],N[m2],N[m3],N[m4],ep}
:ArgumentTypes: {Integer,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fD0v2
:Pattern: D0v2[i_Integer,j_Integer,p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,p4_?NumericQ,s12_?NumericQ,s23_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,m4_?NumericQ,ep_Integer:0]
:Arguments: {i,j,N[p1],N[p2],N[p3],N[p4],N[s12],N[s23],N[m1],N[m2],N[m3],N[m4],ep}
:ArgumentTypes: {Integer,Integer,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fD0v3
:Pattern: D0v3[i_Integer,j_Integer,k_Integer,p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,p4_?NumericQ,s12_?NumericQ,s23_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,m4_?NumericQ,ep_Integer:0]
:Arguments: {i,j,k,N[p1],N[p2],N[p3],N[p4],N[s12],N[s23],N[m1],N[m2],N[m3],N[m4],ep}
:ArgumentTypes: {Integer,Integer,Integer,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fD0v4
:Pattern: D0v4[i_Integer,j_Integer,k_Integer,l_Integer,p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,p4_?NumericQ,s12_?NumericQ,s23_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,m4_?NumericQ,ep_Integer:0]
:Arguments: {i,j,k,l,N[p1],N[p2],N[p3],N[p4],N[s12],N[s23],N[m1],N[m2],N[m3],N[m4],ep}
:ArgumentTypes: {Integer,Integer,Integer,Integer,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

// ==================== PENTAGONS =========================

:Begin:
:Function: fE0v0
:Pattern: E0v0[p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,p4_?NumericQ,p5_?NumericQ,s12_?NumericQ,s23_?NumericQ,s34_?NumericQ,s45_?NumericQ,s15_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,m4_?NumericQ,m5_?NumericQ,ep_Integer:0]
:Arguments: {N[p1],N[p2],N[p3],N[p4],N[p5],N[s12],N[s23],N[s34],N[s45],N[s15],N[m1],N[m2],N[m3],N[m4],N[m5],ep}
:ArgumentTypes: {Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fE0v1
:Pattern: E0v1[i_Integer,p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,p4_?NumericQ,p5_?NumericQ,s12_?NumericQ,s23_?NumericQ,s34_?NumericQ,s45_?NumericQ,s15_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,m4_?NumericQ,m5_?NumericQ,ep_Integer:0]
:Arguments: {i,N[p1],N[p2],N[p3],N[p4],N[p5],N[s12],N[s23],N[s34],N[s45],N[s15],N[m1],N[m2],N[m3],N[m4],N[m5],ep}
:ArgumentTypes: {Integer,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fE0v2
:Pattern: E0v2[i_Integer,j_Integer,p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,p4_?NumericQ,p5_?NumericQ,s12_?NumericQ,s23_?NumericQ,s34_?NumericQ,s45_?NumericQ,s15_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,m4_?NumericQ,m5_?NumericQ,ep_Integer:0]
:Arguments: {i,j,N[p1],N[p2],N[p3],N[p4],N[p5],N[s12],N[s23],N[s34],N[s45],N[s15],N[m1],N[m2],N[m3],N[m4],N[m5],ep}
:ArgumentTypes: {Integer,Integer,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fE0v3
:Pattern: E0v3[i_Integer,j_Integer,k_Integer,p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,p4_?NumericQ,p5_?NumericQ,s12_?NumericQ,s23_?NumericQ,s34_?NumericQ,s45_?NumericQ,s15_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,m4_?NumericQ,m5_?NumericQ,ep_Integer:0]
:Arguments: {i,j,k,N[p1],N[p2],N[p3],N[p4],N[p5],N[s12],N[s23],N[s34],N[s45],N[s15],N[m1],N[m2],N[m3],N[m4],N[m5],ep}
:ArgumentTypes: {Integer,Integer,Integer,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fE0v4
:Pattern: E0v4[i_Integer,j_Integer,k_Integer,l_Integer,p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,p4_?NumericQ,p5_?NumericQ,s12_?NumericQ,s23_?NumericQ,s34_?NumericQ,s45_?NumericQ,s15_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,m4_?NumericQ,m5_?NumericQ,ep_Integer:0]
:Arguments: {i,j,k,l,N[p1],N[p2],N[p3],N[p4],N[p5],N[s12],N[s23],N[s34],N[s45],N[s15],N[m1],N[m2],N[m3],N[m4],N[m5],ep}
:ArgumentTypes: {Integer,Integer,Integer,Integer,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:

:Begin:
:Function: fE0v5
:Pattern: E0v5[i_Integer,j_Integer,k_Integer,l_Integer,m_Integer,p1_?NumericQ,p2_?NumericQ,p3_?NumericQ,p4_?NumericQ,p5_?NumericQ,s12_?NumericQ,s23_?NumericQ,s34_?NumericQ,s45_?NumericQ,s15_?NumericQ,m1_?NumericQ,m2_?NumericQ,m3_?NumericQ,m4_?NumericQ,m5_?NumericQ,ep_Integer:0]
:Arguments: {i,j,k,l,m,N[p1],N[p2],N[p3],N[p4],N[p5],N[s12],N[s23],N[s34],N[s45],N[s15],N[m1],N[m2],N[m3],N[m4],N[m5],ep}
:ArgumentTypes: {Integer,Integer,Integer,Integer,Integer,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Real64,Integer}
:ReturnType: Manual
:End:


:Evaluate: End[]
:Evaluate: EndPackage[]

// ====================================================================

using namespace std;

double fGetMu2() {
    return PJFry::GetMu2();
}

double fSetMu2(double newmu2) {
    return PJFry::SetMu2(newmu2);
}

void fA0v0(double m1, int ep) {
    complex<double> result = PJFry::A0v0(m1,ep);
    MLPutFunction(stdlink, "Complex", 2);
    MLPutReal64(stdlink, result.real());
    MLPutReal64(stdlink, result.imag());
}

void fB0v0(double p1,
           double m1,  double m2, int ep) {
    complex<double> result = PJFry::B0v0(p1,m1,m2,ep);
    MLPutFunction(stdlink, "Complex", 2);
    MLPutReal64(stdlink, result.real());
    MLPutReal64(stdlink, result.imag());
}

void fB0v1(int i,
           double p1,
           double m1,  double m2, int ep) {
    complex<double> result = PJFry::B0v1(i,p1,m1,m2,ep);
    MLPutFunction(stdlink, "Complex", 2);
    MLPutReal64(stdlink, result.real());
    MLPutReal64(stdlink, result.imag());
}

void fB0v2(int i, int j,
           double p1,
           double m1,  double m2, int ep) {
    complex<double> result = PJFry::B0v2(i,j,p1,m1,m2,ep);
    MLPutFunction(stdlink, "Complex", 2);
    MLPutReal64(stdlink, result.real());
    MLPutReal64(stdlink, result.imag());
}

void fC0v0(double p1,  double p2,  double p3,
           double m1,  double m2,  double m3, int ep) {
    complex<double> result = PJFry::C0v0(p1,p2,p3,m1,m2,m3,ep);
    MLPutFunction(stdlink, "Complex", 2);
    MLPutReal64(stdlink, result.real());
    MLPutReal64(stdlink, result.imag());
}

void fC0v1(int i,
           double p1,  double p2,  double p3,
           double m1,  double m2,  double m3, int ep) {
    complex<double> result = PJFry::C0v1(i,p1,p2,p3,m1,m2,m3,ep);
    MLPutFunction(stdlink, "Complex", 2);
    MLPutReal64(stdlink, result.real());
    MLPutReal64(stdlink, result.imag());
}

void fC0v2(int i, int j,
           double p1,  double p2,  double p3,
           double m1,  double m2,  double m3, int ep) {
    complex<double> result = PJFry::C0v2(i,j,p1,p2,p3,m1,m2,m3,ep);
    MLPutFunction(stdlink, "Complex", 2);
    MLPutReal64(stdlink, result.real());
    MLPutReal64(stdlink, result.imag());
}

void fC0v3(int i, int j, int k,
           double p1,  double p2,  double p3,
           double m1,  double m2,  double m3, int ep) {
    complex<double> result = PJFry::C0v3(i,j,k,p1,p2,p3,m1,m2,m3,ep);
    MLPutFunction(stdlink, "Complex", 2);
    MLPutReal64(stdlink, result.real());
    MLPutReal64(stdlink, result.imag());
}

void fD0v0(double p1,  double p2,  double p3,  double p4,
           double s12, double s23,
           double m1,  double m2,  double m3,  double m4,  int ep) {
  complex<double> result = PJFry::D0v0(p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,ep);
  MLPutFunction(stdlink, "Complex", 2);
  MLPutReal64(stdlink, result.real());
  MLPutReal64(stdlink, result.imag());
}

void fD0v1(int i,
           double p1,  double p2,  double p3,  double p4,
           double s12, double s23,
           double m1,  double m2,  double m3,  double m4,  int ep) {
  complex<double> result = PJFry::D0v1(i,p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,ep);
  MLPutFunction(stdlink, "Complex", 2);
  MLPutReal64(stdlink, result.real());
  MLPutReal64(stdlink, result.imag());
}

void fD0v2(int i, int j,
           double p1,  double p2,  double p3,  double p4,
           double s12, double s23,
           double m1,  double m2,  double m3,  double m4,  int ep) {
  complex<double> result = PJFry::D0v2(i,j,p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,ep);
  MLPutFunction(stdlink, "Complex", 2);
  MLPutReal64(stdlink, result.real());
  MLPutReal64(stdlink, result.imag());
}

void fD0v3(int i, int j, int k,
           double p1,  double p2,  double p3,  double p4,
           double s12, double s23,
           double m1,  double m2,  double m3,  double m4,  int ep) {
  complex<double> result = PJFry::D0v3(i,j,k,p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,ep);
  MLPutFunction(stdlink, "Complex", 2);
  MLPutReal64(stdlink, result.real());
  MLPutReal64(stdlink, result.imag());
}

void fD0v4(int i, int j, int k, int l,
           double p1,  double p2,  double p3,  double p4,
           double s12, double s23,
           double m1,  double m2,  double m3,  double m4,  int ep) {
  complex<double> result = PJFry::D0v4(i,j,k,l,p1,p2,p3,p4,s12,s23,m1,m2,m3,m4,ep);
  MLPutFunction(stdlink, "Complex", 2);
  MLPutReal64(stdlink, result.real());
  MLPutReal64(stdlink, result.imag());
}

void fE0v0(double p1,  double p2,  double p3,  double p4,  double p5,
           double s12, double s23, double s34, double s45, double s15,
           double m1,  double m2,  double m3,  double m4,  double m5, int ep) {
  complex<double> result = PJFry::E0v0(p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,m1,m2,m3,m4,m5,ep);
  MLPutFunction(stdlink, "Complex", 2);
  MLPutReal64(stdlink, result.real());
  MLPutReal64(stdlink, result.imag());
}

void fE0v1(int i,
           double p1,  double p2,  double p3,  double p4,  double p5,
           double s12, double s23, double s34, double s45, double s15,
           double m1,  double m2,  double m3,  double m4,  double m5, int ep) {
  complex<double> result = PJFry::E0v1(i,p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,m1,m2,m3,m4,m5,ep);
  MLPutFunction(stdlink, "Complex", 2);
  MLPutReal64(stdlink, result.real());
  MLPutReal64(stdlink, result.imag());
}

void fE0v2(int i, int j,
           double p1,  double p2,  double p3,  double p4,  double p5,
           double s12, double s23, double s34, double s45, double s15,
           double m1,  double m2,  double m3,  double m4,  double m5, int ep) {
  complex<double> result = PJFry::E0v2(i,j,p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,m1,m2,m3,m4,m5,ep);
  MLPutFunction(stdlink, "Complex", 2);
  MLPutReal64(stdlink, result.real());
  MLPutReal64(stdlink, result.imag());
}

void fE0v3(int i, int j, int k,
           double p1,  double p2,  double p3,  double p4,  double p5,
           double s12, double s23, double s34, double s45, double s15,
           double m1,  double m2,  double m3,  double m4,  double m5, int ep) {
  complex<double> result = PJFry::E0v3(i,j,k,p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,m1,m2,m3,m4,m5,ep);
  MLPutFunction(stdlink, "Complex", 2);
  MLPutReal64(stdlink, result.real());
  MLPutReal64(stdlink, result.imag());
}

void fE0v4(int i, int j, int k, int l,
           double p1,  double p2,  double p3,  double p4,  double p5,
           double s12, double s23, double s34, double s45, double s15,
           double m1,  double m2,  double m3,  double m4,  double m5, int ep) {
  complex<double> result = PJFry::E0v4(i,j,k,l,p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,m1,m2,m3,m4,m5,ep);
  MLPutFunction(stdlink, "Complex", 2);
  MLPutReal64(stdlink, result.real());
  MLPutReal64(stdlink, result.imag());
}

void fE0v5(int i, int j, int k, int l, int m,
           double p1,  double p2,  double p3,  double p4,  double p5,
           double s12, double s23, double s34, double s45, double s15,
           double m1,  double m2,  double m3,  double m4,  double m5, int ep) {
  complex<double> result = PJFry::E0v5(i,j,k,l,m,p1,p2,p3,p4,p5,s12,s23,s34,s45,s15,m1,m2,m3,m4,m5,ep);
  MLPutFunction(stdlink, "Complex", 2);
  MLPutReal64(stdlink, result.real());
  MLPutReal64(stdlink, result.imag());
}


// ====================================================================

// extern "C" void rrinit_();
int main(int argc, char* argv[]) {

  return MLMain(argc, argv);
}
