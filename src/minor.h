/*
 * minor.h - signed minors classes for PJFry library
 *
 * Valery Yundin <valery.yundin@desy.de>
 */

#ifndef QUL_MINOR_H
#define QUL_MINOR_H

#include "common.h"
#include "kinem.h"
#include "pointer.h"
#include <bitset>

#define tswap(x,y,t) \
    if (x > y) { \
    t=y; \
    y=x; \
    x=t; \
    }

class MinorBase : public SRefCnt
{
  public:
    MinorBase() {};

    // Symmetric index - start from 1
    inline static int ns(int i, int j) CONST
    {
      return ( i<=j ? (i-1)+((j-1)*j)/2 : (j-1)+((i-1)*i)/2 );
    };

    inline static int nss(int i, int j) CONST // ordered
    {
      return (i-1)+((j-1)*j)/2;
    };

    // Symmetric index - generic
    inline static int is(int i, int j) CONST
    {
      return ( i<=j ? i+j*(j+1)/2 : j+i*(i+1)/2 );
    };

    inline static int is(int i, int j, int k) CONST
    {
      if (i <= j) {
        return (j <= k ? i+ti2[j]+ti3[k] : is(i,k)+ti3[j]);
      }
      else {
        return (i >  k ? is(j,k)+ti3[i] : j+ti2[i]+ti3[k]);
      }
    };

    inline static int is(int i, int j, int k, int l) CONST
    {
      if (i <= j) {
        if (j <= k) {
          return (k <= l ? i+ti2[j]+ti3[k]+ti4[l]
                         : is(i,j,l)+ti4[k]  );
        }
        else {
          return (j >  l ? is(i,k,l)+ti4[j]
                         : is(i,k)+ti3[j]+ti4[l]  );
        }
      }
      else {
        if (i > k) {
          return (i >  l ? is(j,k,l)+ti4[i]
                         : is(j,k)+ti3[i]+ti4[l]  );
        }
        else {
          return (k <= l ? j+ti2[i]+ti3[k]+ti4[l]
                         : is(i,j,l)+ti4[k]  );
        }
      }
    };

    inline static int iss(int i, int j) CONST // ordered
    {
      assert(i<=j);
      return i+j*(j+1)/2;
    }

    inline static int iss(int i, int j, int k) CONST // ordered
    {
      assert(i <= j && j <= k);
      return i+ti2[j]+ti3[k];
    }

    inline static int iss(int i, int j, int k, int l) CONST // ordered
    {
      assert(i <= j && j <= k && k <= l);
      return i+ti2[j]+ti3[k]+ti4[l];
    }

    inline static int iss(int i, int j, int k, int l, int m) CONST // ordered
    {
      assert(i <= j && j <= k && k <= l && l <= m);
      return i+ti2[j]+ti3[k]+ti4[l]+ti5[m];
    }

    inline static double getmeps()
    {
      return meps;
    }

    // Utility functions
    static int im3(int i, int j, int k) CONST; // Antisymmetric index for "3 out of 6" minors
    static int im2(int i, int j) CONST;        // Antisymmetric index for "2 out of 6" minors
    static int signM3ud(int i, int j, int k, int l, int m, int n) CONST; // Signature[{i,j,k}]*Signature[{l,m,n}]
    static int signM2ud(int i, int j, int l, int m) CONST;               // Signature[{i,j}]*Signature[{l,m}]

    // fill 'free' array with indices which are not occupied by 'set' array
    static void freeidxM3(int set[], int free[]);
  private:
    static const unsigned char ti2[8];
    static const unsigned char ti3[8];
    static const unsigned char ti4[8];
    static const unsigned char ti5[8];

  protected:
    static const unsigned char idxtbl[64];

    static const double teps; // expansion target accuracy

    static const double ceps;

    static const double deps1;
    static const double deps2;
    static const double deps3;

    static const double deps;
    static const double meps; // onshell cutoff
};

template <int N>
class Minor : public MinorBase
{
  public:
    Minor() {};

    inline double Kay(int i, int j) PURE
    {
      if (i==0) {
        return j==0 ? 0 : 1;
      } else {
        return j==0 ? 1 : Cay[ns(i,j)];
      }
    };

  protected:
    // Cayley matrix (symmetric)
    static const int DCay=N+1;
    double Cay[(DCay-1)*(DCay)/2];
};

class Minor5 : public Minor<5>
{
  public:
    friend class SPtr<Minor5>;
    typedef SPtr<Minor5> Ptr;
    static Ptr create(const Kinem5 &k) { return Ptr(new Minor5(k)); };
    static Ptr create(const Kinem4 &k) { return Ptr(new Minor5(k)); };

    ncomplex evalE(int ep);
    ncomplex evalE(int ep, int i);
    ncomplex evalE(int ep, int i, int j);
    ncomplex evalE(int ep, int i, int j, int k);
    ncomplex evalE(int ep, int i, int j, int k, int l);
    ncomplex evalE(int ep, int i, int j, int k, int l, int m);

    ncomplex I4s(int ep, int s);

    ncomplex I3st(int ep, int s, int t);
    ncomplex I4Ds(int ep, int s);
    ncomplex I4Dsi(int ep, int s, int i);

    ncomplex I2stu(int ep, int s, int t, int u);
    ncomplex I3Dst(int ep, int s, int t);
    ncomplex I4D2s(int ep, int s);

    ncomplex I4D2si(int ep, int s, int i);
    ncomplex I3Dsti(int ep, int s, int t, int i);
    ncomplex I3D2st(int ep, int s, int t);
    ncomplex I4D3s(int ep, int s);
    ncomplex I4D2sij(int ep, int s, int i, int j);

    ncomplex I2Dstu(int ep, int s, int t, int u);
    ncomplex I3D2sti(int ep, int s, int t, int i);
    ncomplex I4D3si(int ep, int s, int i);
    ncomplex I4D3sij(int ep, int s, int i, int j);

    ncomplex I2Dstui(int ep, int s, int t, int u, int i);
    ncomplex I3D2stij(int ep, int s, int t, int i, int j);
    ncomplex I4D3sijk(int ep, int s, int i, int j, int k);

    ncomplex I4D4s(int ep, int s);
    ncomplex I4D4si(int ep, int s, int i);
    ncomplex I3D3sti(int ep, int s, int t, int i);
    ncomplex I4D4sij(int ep, int s, int i, int j);
    ncomplex I2D2stui(int ep, int s, int t, int u, int i);
    ncomplex I3D3stij(int ep, int s, int t, int i, int j);
    ncomplex I4D4sijk(int ep, int s, int i, int j, int k);
    ncomplex I2D2stuij(int ep, int s, int t, int u, int i, int j);
    ncomplex I3D3stijk(int ep, int s, int t, int i, int j, int k);
    ncomplex I4D4sijkl(int ep, int s, int i, int j, int k, int l);

    // Gram4
    ncomplex I2D2stu(int ep, int s, int t, int u);
    ncomplex I3D3st(int ep, int s, int t);
    ncomplex I2D3stu(int ep, int s, int t, int u);
    ncomplex I3D4st(int ep, int s, int t);
    ncomplex I2D4stu(int ep, int s, int t, int u);
    ncomplex I3D5st(int ep, int s, int t);
    ncomplex I2D5stu(int ep, int s, int t, int u);
    ncomplex I3D6st(int ep, int s, int t);
    ncomplex I2D6stu(int ep, int s, int t, int u);
    ncomplex I3D7st(int ep, int s, int t);

    //Aux

    double M1(int i, int l) PURE;
    double M2(int i, int j, int l, int m) PURE;
    double M3(int i, int j, int k, int l, int m, int n) PURE;

    double maxS4(int s) PURE;
    double maxS3(int s, int t) PURE;

    double gram3(double p1, double p2, double p3) PURE;

  private:
    Minor5(const Kinem5 &k);                        // prevent direct creation
    Minor5(const Kinem4 &k);                        // prevent direct creation
    Minor5(const Minor5 &m) { assert(0); };             // prevent copy-constructing
    Minor5& operator= (const Minor5& m) { assert(0); }; // prevent reassignment

    Kinem5 kinem;
    int smax;

    // Maximal elements of sub-Cayley's
    double pmaxS4[DCay-1];
    double pmaxS3[10];     // symm 5x5 w/o diag els

    // flags marking evuation steps
    enum EvalM {E_None=0,
                E_M1, E_M2, E_M3,
                                      E_I4D4sijkl,
                           E_I4D4sijk=E_I4D4sijkl+3,
                E_I3D3stij=E_I4D4sijk+3,
                                                              E_I4D4sij=E_I3D3stij+3,
                                                    E_I3D3sti=E_I4D4sij+3,
                                           E_I4D4si=E_I3D3sti+3,
                                   E_I4D4s=E_I4D4si+3,
                         E_I2D6stu=E_I4D4s+3,
                E_I3D7st=E_I2D6stu+3,
                                                               E_I2D5stu=E_I3D7st+3,
                                                      E_I3D6st=E_I2D5stu+3,
                                            E_I2D4stu=E_I3D6st+3,
                                   E_I3D5st=E_I2D4stu+3,
                         E_I2D3stu=E_I3D5st+3,
                E_I3D4st=E_I2D3stu+3,
                                                          E_I3D3stijk=E_I3D4st+3,
                                               E_I2D2stui=E_I3D3stijk+3,
                                   E_I2D2stuij=E_I2D2stui+3,
                         E_I2D2stu=E_I2D2stuij+3,
                E_I3D3st=E_I2D2stu+3,
                                                      E_I4D3sijk=E_I3D3st+3,
                                           E_I3D2stij=E_I4D3sijk+3,
                                 E_I2Dstui=E_I3D2stij+3,
                        E_I3D2st=E_I2Dstui+2,
                E_I4D3s=E_I3D2st+3,
                                                               E_I4D3sij=E_I4D3s+3,
                                                      E_I4D3si=E_I4D3sij+3,
                                            E_I3D2sti=E_I4D3si+3,
                                   E_I2Dstu=E_I3D2sti+3,
                          E_I3Dsti=E_I2Dstu+3,
                E_I4D2sij=E_I3Dsti+3,
                                                        E_I2stu=E_I4D2sij+3,
                                                E_I4D2s=E_I2stu+3,
                                       E_I4D2si=E_I4D2s+3,
                               E_I3Dst=E_I4D2si+3,
                       E_I4Dsi=E_I3Dst+3,
                E_I4Ds=E_I4Dsi+3,
                             E_I4s=E_I4Ds+3,
                      E_I3st=E_I4s+3,
                E_DUM=E_I3st+3,
                E_LEN};
    std::bitset<E_LEN> fEval;

    // number of unique ways you can scratch out 3 rows in 6x6 matrix
    static const int DM3=20; // Binomial[6,3] = 6*5*4/3! = 20
    double pM3[DM3*(DM3+1)/2];
    static const int DM2=15; // Binomial[6,2] = 15
    double pM2[DM2*(DM2+1)/2];
    static const int DM1=6; // Binomial[6,1] = 6
    double pM1[DM1*(DM1+1)/2];

    // Integral values
    ncomplex pI4s[3][DCay-1];           // s{1..5}

    ncomplex pI3st[3][10];              // symm 5x5 w/o diag els
    ncomplex pI4Ds[1][DCay-1];          // s{1..5}                        // finite
    ncomplex pI4Dsi[3][DCay-2][DCay-1]; // i{1..4}, s{1..5}

    ncomplex pI2stu[2][10];             // symm 5x5x5 w/o diag els        // (0,1) parts only
    ncomplex pI3Dst[1][10];             // symm 5x5 w/o diag els 5*6-5=10 // finite part only
    ncomplex pI4D2s[1][DCay-1];         // s{1..5}                        // finite part only

    ncomplex pI4D2si[1][DCay-2][DCay-1]; // i{1..4} s{1..5}               // finite
    ncomplex pI3Dsti[3][DCay-2][10];     // i{1..4} + symm 5x5 w/o diag els
    ncomplex pI4D2sij[3][(DCay-2)*(DCay-1)/2][DCay-1];  // symm 4x4, s{1..5}

    ncomplex pI2Dstu[2][10];              // symm 5x5x5 w/o diag els      // (0,1) parts only
    ncomplex pI3D2st[2][10];              // symm 5x5 w/o diag els        // (0,1) parts only
    ncomplex pI4D3s[2][DCay-1];           // s{1..5}                      // (0,1) parts only
    ncomplex pI3D2sti[2][DCay-2][10];     // i{1..4} + symm 5x5 w/o diag els // (0,1) parts only
    ncomplex pI4D3si[2][DCay-2][DCay-1];  // i{1..4} s{1..5}              // (0,1) parts only
    ncomplex pI4D3sij[1][(DCay-2)*(DCay-1)/2][DCay-1];  // symm 4x4, s{1..5} // finite

    ncomplex pI2Dstui[2][DCay-2][DCay-1];     // ~~ 16 elements           // finite part only
    ncomplex pI3D2stij[3][(DCay-2)*(DCay-1)/2][10];  // 'symm 4x4' x 'symm 5x5 w/o diag els'
    ncomplex pI4D3sijk[3][20][DCay-1];     // 4x4x4 symm(20 els), s{1..5}

    ncomplex pI2D2stu[2][10];             // symm 5x5x5 w/o diag els      // (0,1) parts only
    ncomplex pI3D3st[2][10];              // symm 5x5 w/o diag els        // (0,1) parts only
    ncomplex pI4D4s[2][DCay-1];           // s{1..5}                      // (0,1) parts only
    ncomplex pI4D4si[2][DCay-2][DCay-1];  // i{1..4} s{1..5}              // (0,1) parts only
    ncomplex pI3D3sti[2][DCay-2][10];     // i{1..4} + symm 5x5 w/o diag els // (0,1) parts only
    ncomplex pI4D4sij[2][(DCay-2)*(DCay-1)/2][DCay-1];  // symm 4x4, s{1..5} // (0,1) parts only
    ncomplex pI2D2stui[2][DCay-2][DCay-1];     // ~~ 16 elements           // (0,1) parts only
    ncomplex pI3D3stij[2][(DCay-2)*(DCay-1)/2][10];  // 'symm 4x4' x 'symm 5x5 w/o diag els' // (0,1) parts only
    ncomplex pI4D4sijk[2][20][DCay-1];     // 4x4x4 symm(20 els), s{1..5} // finite
    ncomplex pI2D2stuij[2][DCay-2][DCay-1][2]; // as stui but with 0:i==j and 1:i!=j
    ncomplex pI3D3stijk[3][20][10];            // 4x4x4 symm(20 els), x 'symm 5x5 w/o diag els'
    ncomplex pI4D4sijkl[3][35][DCay-1];    // 4x4x4x4 symm(35 els), s{1..5} // ????

    // Gram4
    // TODO cleanup these funcs
    ncomplex pI2D3stu[2][10];             // symm 5x5x5 w/o diag els      // (0,1) parts only
    ncomplex pI3D4st[2][10];              // symm 5x5 w/o diag els        // (0,1) parts only
    ncomplex pI2D4stu[2][10];             // symm 5x5x5 w/o diag els      // (0,1) parts only
    ncomplex pI3D5st[2][10];              // symm 5x5 w/o diag els        // (0,1) parts only
    ncomplex pI2D5stu[2][10];             // symm 5x5x5 w/o diag els      // (0,1) parts only
    ncomplex pI3D6st[2][10];              // symm 5x5 w/o diag els        // (0,1) parts only
    ncomplex pI2D6stu[2][10];             // symm 5x5x5 w/o diag els      // (0,1) parts only
    ncomplex pI3D7st[2][10];              // symm 5x5 w/o diag els        // (0,1) parts only

    // Aux

    // evaluate and save for future use M1,M2,M3 minors
    void evalM1();
    void evalM2();
    void evalM3();

    // find and save maximal elements of subcayley matrices
    void maxCay();

    // evaluate and save for the future use scratched integrals
    void I4sEval(int ep);

    void I3stEval(int ep);
    void I4DsEval(int ep);
    void I4DsiEval(int ep);

    void I2stuEval(int ep);
    void I3DstEval(int ep);
    void I4D2sEval(int ep);

    void I4D2siEval(int ep);
    void I3DstiEval(int ep);
    void I3D2stEval(int ep);
    void I4D3sEval(int ep);
    void I4D2sijEval(int ep);

    void I2DstuEval(int idx, int ep, int s, int t, int u, int m, int n, double qsq);
    void I3D2stiEval(int ep);
    void I4D3siEval(int ep);
    void I4D3sijEval(int ep);

    void I2DstuiEval(int ep, int s, int t, int u, int i, int ip, double qsq);
    void I3D2stijEval(int ep);
    void I4D3sijkEval(int ep);

    void I4D4sEval(int ep);
    void I4D4siEval(int ep);
    void I3D3stiEval(int ep);
    void I4D4sijEval(int ep);
    void I2D2stuiEval(int ep, int s, int t, int u, int i, int ip, double qsq);
    void I3D3stijEval(int ep);
    void I4D4sijkEval(int ep);
    void I2D2stuijEval(int ep, int s, int t, int u, int i, int ip, double qsq);
    void I3D3stijkEval(int ep);
    void I4D4sijklEval(int ep);

    // Gram4
    void I2D2stuEval(int idx, int ep, int s, int t, int u, int m, int n, double qsq);
    void I3D3stEval(int ep);
    void I2D3stuEval(int idx, int ep, int s, int t, int u, int m, int n, double qsq);
    void I3D4stEval(int ep);
    void I2D4stuEval(int idx, int ep, int s, int t, int u, int m, int n, double qsq);
    void I3D5stEval(int ep);
    void I2D5stuEval(int idx, int ep, int s, int t, int u, int m, int n, double qsq);
    void I3D6stEval(int ep);
    void I2D6stuEval(int idx, int ep, int s, int t, int u, int m, int n, double qsq);
    void I3D7stEval(int ep);
};

class Minor4 : public Minor<4>
{
  public:
    friend class SPtr<Minor4>;
    typedef SPtr<Minor4> Ptr;
    static Ptr create(const Kinem4 &k, Minor5::Ptr mptr5, int s)
    {
      return Ptr(new Minor4(k, mptr5, s));
    };

    ncomplex evalD(int ep);
    ncomplex evalD(int ep, int i);
    ncomplex evalD(int ep, int i, int j);
    ncomplex evalD(int ep, int i, int j, int k);
    ncomplex evalD(int ep, int i, int j, int k, int l);

  private:
    Minor4(const Kinem4 &k, Minor5::Ptr mptr, int s);

    Kinem4 kinem;

    Minor5::Ptr pm5;
    int ps;
};

class Minor3 : public Minor<3>
{
  public:
    friend class SPtr<Minor3>;
    typedef SPtr<Minor3> Ptr;
    static Ptr create(const Kinem3 &k, Minor5::Ptr mptr5, int s, int t)
    {
      return Ptr(new Minor3(k, mptr5, s, t));
    };

    ncomplex evalC(int ep);
    ncomplex evalC(int ep, int i);
    ncomplex evalC(int ep, int i, int j);
    ncomplex evalC(int ep, int i, int j, int k);

  private:
    Minor3(const Kinem3 &k);
    Minor3(const Kinem3 &k, Minor5::Ptr mptr5, int s, int t);

    Kinem3 kinem;

    Minor5::Ptr pm5;
    int ps, pt;
};

class Minor2 : public Minor<2>
{
  public:
    friend class SPtr<Minor2>;
    typedef SPtr<Minor2> Ptr;
    static Ptr create(const Kinem2 &k, Minor5::Ptr mptr5, int s, int t, int u)
    {
      return Ptr(new Minor2(k, mptr5, s, t, u));
    };

    ncomplex evalB(int ep);
    ncomplex evalB(int ep, int i);
    ncomplex evalB(int ep, int i, int j);

  private:
    Minor2(const Kinem2 &k);
    Minor2(const Kinem2 &k, Minor5::Ptr mptr5, int s, int t, int u);

    Kinem2 kinem;

    Minor5::Ptr pm5;
    int ps, pt, pu;
};


/* ===============================================
 *
 *               Utility functions
 *
 * ===============================================
 */

// Completely antysymmetric i,j,k size 6 matrix index
inline
int MinorBase::im3(int i, int j, int k)
{
  return idxtbl[(1<<i)|(1<<j)|(1<<k)];
}

// Completely antysymmetric i,j size 6 matrix index
inline
int MinorBase::im2(int i, int j)
{
  return idxtbl[(1<<i)|(1<<j)];
}

// Signature[{i,j,k}]*Signature[{l,m,n}]
inline
int MinorBase::signM3ud(int i, int j, int k, int l, int m, int n)
{
  int t=(j-i)*(k-j)*(k-i)*(m-l)*(n-m)*(n-l);
  return t==0 ? 0 : 2*(t>>(sizeof(int)*8-1))+1;
}

// Signature[{i,j}]*Signature[{l,m}]
inline
int MinorBase::signM2ud(int i, int j, int l, int m)
{
  int t=(j-i)*(m-l);
  return t==0 ? 0 : 2*(t>>(sizeof(int)*8-1))+1;
}

#endif /* QUL_MINOR_H */
