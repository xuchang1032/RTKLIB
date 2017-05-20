/* PSINS c++ hearder file PSINS.h */
/*
	By     : Yan Gongmin @ NWPU
	Date   : 2015-02-17, 2015-04-25, 2017-04-29
	From   : College of Automation, 
	         Northwestern Polytechnical University, 
			 Xi'an 710072, China
*/

#ifndef _PSINS_H
#define _PSINS_H

#ifdef _MSC_VER
#ifdef PSINS_DLL
#define PSINS_API __declspec(dllexport)
#else
#define PSINS_API __declspec(dllimport)
#endif
#else
#define PSINS_API
#endif


#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <time.h>

/* type re-define */
#ifndef BOOL
typedef int		BOOL;
#endif

/*constant define*/
#ifndef TRUE
#define TRUE	1
#define FALSE	0
#endif

#ifndef NULL
#define NULL	((void *)0)
#endif

#ifndef PI
#define PI		3.14159265358979
#endif

#define PI_2	(PI/2.0)
#define PI_4	(PI/4.0)
#define _2PI	(2.0*PI)

#define EPS		2.22044604925031e-16
#define INF		123456789.0e100

BOOL	assert(BOOL b);
int		sign(double val, double eps = EPS);
double	range(double val, double minVal, double maxVal);
double	atan2Ex(double y, double x);

#define asinEx(x)		asin(range(x, -1.0, 1.0))

#define acosEx(x)		acos(range(x, -1.0, 1.0))

#ifndef max
#define max(x,y)        { (x)>=(y)?(x):(y); }
#endif

#ifndef min
#define min(x,y)        { (x)<=(y)?(x):(y); }
#endif

// class define
class CGLV;
class CVect3;	class CMat3;	class CQuat;
class CEarth;	class CIMU;		class CSINS;	class CAligni0;
class CVect;	class CMat;		class CKalman;	class CTDKF;

// Max Matrix Dimension define
#define MMD		16
#define MMD2	(MMD*MMD)


// global variables and functions, can not be changed in any way
extern const CVect3	O31;
extern const CQuat	qI;
extern const CMat3	I33;
extern const CGLV	glv;


class PSINS_API CGLV
{
public:
	double Re, f, g0, wie;											// the Earth's parameters
	double e, e2;
	double mg, ug, deg, min, sec, hur, ppm, ppmpsh;					// commonly used units
	double dps, dph, dpsh, dphpsh, ugpsh, ugpsHz, mpsh, mpspsh, secpsh;

	CGLV(double Re=6378137.0, double f=(1.0/298.257), double wie0=7.2921151467e-5, double g0=9.7803267714);
};

class PSINS_API CVect3
{
public:
	double i, j, k;

	CVect3(void);
	CVect3(double xx, double yy=0.0, double zz=0.0);
	CVect3(const double *pdata);
	CVect3(const CMat3 &Cnb);								// DCM to Euler angles
	CVect3(const CQuat &qnb);								// quaterion to  Euler angles

	CVect3 operator+(const CVect3 &v);						// vector addition
	CVect3 operator-(const CVect3 &v);						// vector subtraction
	CVect3 operator*(const CVect3 &v);						// vector cross multiplication
	CVect3 operator*(double f);								// vector multiply scale
	CVect3 operator/(double f);								// vector divide scale
	CVect3& operator+=(const CVect3 &v);					// vector addition
	CVect3& operator-=(const CVect3 &v);					// vector subtraction
	CVect3& operator*=(double f);							// vector multiply scale
	CVect3& operator/=(double f);							// vector divide scale
	BOOL IsZero(double eps=EPS);							// assert if all elements are zeros
	BOOL IsZeroXY(double eps=EPS);							// assert if x&&y-elements are zeros
	friend CVect3 operator*(double f, const CVect3 &v);		// scale multiply vector
	friend CVect3 operator-(const CVect3 &v);				// minus
	friend double norm(const CVect3 &v);					// vector norm
	friend double normXY(const CVect3 &v);					// vector norm or X & Y components
	friend double dot(const CVect3 &v1, const CVect3 &v2);	// vector dot multiplication
	friend CQuat rv2q(const CVect3 &rv);					// rotation vector to quaternion
	friend CMat3 askew(const CVect3 &v);					// askew matrix;
	friend CMat3 pos2Cen(const CVect3 &pos);				// to geographical position matrix
	friend CVect3 pp2vn(const CVect3 &pos1, const CVect3 &pos0, double ts=1.0, CEarth *pEth=NULL);  // position difference to velocity
};

class PSINS_API CQuat
{
public:
	double q0, q1, q2, q3;

	CQuat(void);
	CQuat(double qq0, double qq1=0.0, double qq2=0.0, double qq3=0.0);
	CQuat(const double *pdata);
	CQuat(const CVect3 &att);					// Euler angles to quaternion
	CQuat(const CMat3 &Cnb);					// DCM to quaternion

	CQuat operator+(const CVect3 &phi);			// true quaternion add misalign angles
	CQuat operator-(const CVect3 &phi);			// calculated quaternion delete misalign angles
	CVect3 operator-(CQuat &quat);				// get misalign angles from calculated quaternion & true quaternion
	CQuat operator*(const CQuat &q);			// quaternion multiplication
	CVect3 operator*(const CVect3 &v);			// quaternion multiply vector
	CQuat& operator*=(const CQuat &q);			// quaternion multiplication
	CQuat& operator-=(const CVect3 &phi);		// calculated quaternion delete misalign angles
	void normlize(CQuat *q);					// quaternion norm
	friend CQuat operator~(const CQuat &q);		// quaternion conjugate
	friend CVect3 q2rv(const CQuat &q);			// quaternion to rotation vector
};

class PSINS_API CMat3
{
public:
	double e00, e01, e02, e10, e11, e12, e20, e21, e22;

	CMat3(void);
	CMat3(double xx, double xy, double xz,
		  double yx, double yy, double yz,
		  double zx, double zy, double zz );
	CMat3(const CVect3 &v0, const CVect3 &v1, const CVect3 &v2);  // M = [v0; v1; v2]
	CMat3(const CVect3 &att);								// Euler angles to DCM
	CMat3(const CQuat &qnb);								// quaternion to DCM

	CMat3 operator+(const CMat3 &m);						// matirx addition
	CMat3 operator-(const CMat3 &m);						// matirx subtraction
	CMat3 operator*(const CMat3 &m);						// matirx multiplication
	CMat3 operator*(double f);								// matirx multiply scale
	CVect3 operator*(const CVect3 &v);						// matirx multiply vector
	friend CMat3 operator~(const CMat3 &m);					// matirx transposition
	friend CMat3 operator*(double f, const CMat3 &m);		// scale multiply matirx
	friend double det(const CMat3 &m);						// matirx determination
	friend CMat3 inv(const CMat3 &m);						// matirx inverse
	friend CVect3 diag(const CMat3 &m);						// diagonal of a matrix
	friend CMat3 diag(const CVect3 &v);						// diagonal matrix
	friend CMat3 dv2att(CVect3 &va1, const CVect3 &va2, CVect3 &vb1, const CVect3 &vb2);  // attitude determination using double-vector
};

class PSINS_API CVect
{
public:
	int row, clm;
	double dd[MMD];

	CVect(void);
	CVect(int row0, int clm0=1);
	CVect(int row0, double f);
	CVect(const CVect3 &v);
	CVect(const CVect3 &v1, const CVect3 v2);

	void Set(double f, ...);
	void Set2(double f, ...);
	CVect operator+(const CVect &v);		// vector addition
	CVect operator-(const CVect &v);		// vector subtraction
	CVect operator*(double f);				// vector multiply scale
	CVect& operator+=(const CVect &v);		// vector addition
	CVect& operator-=(const CVect &v);		// vector subtraction
	CVect& operator*=(double f);			// vector multiply scale
	CVect operator*(const CMat &m);			// row-vector multiply matrix
	CMat operator*(const CVect &v);			// 1xn vector multiply nx1 vector, or nx1 vector multiply 1xn vector
	double& operator()(int r);				// vector element
	friend CVect operator~(const CVect &v);	// vector transposition
	friend double norm(const CVect &v);		// vector norm
};

class PSINS_API CMat
{
public:
	int row, clm, rc;
	double dd[MMD2];

	CMat(void);
	CMat(int row0, int clm0);
	CMat(int row0, int clm0, double f);

	void SetDiag(double f, ...);
	void SetDiag2(double f, ...);
	CMat operator+(const CMat &m);				// matirx addition
	CMat operator-(const CMat &m);				// matirx subtraction
	CMat operator*(double f);					// matirx multiply scale
	CVect operator*(const CVect &v);			// matirx multiply vector
	CMat operator*(const CMat &m);				// matirx multiplication
	CMat& operator+=(const CMat &m0);			// matirx addition
	CMat& operator+=(const CVect &v);			// matirx + diag(vector)
	CMat& operator-=(const CMat &m0);			// matirx subtraction
	CMat& operator*=(double f);					// matirx multiply scale
	CMat& operator++();							// 1.0 + diagonal
	double& operator()(int r, int c);	// get element m(r,c)
	void SetRow(int i, const CVect &v);			// set i-row from vector
	void SetClm(int j, const CVect &v);			// set j-column from vector
	CVect GetRow(int i);						// get i-row from matrix
	CVect GetClm(int j);						// get j-column from matrix
	void SetRowVect3(int i, int j, const CVect3 &v); // set i...(i+2)-row&j-column from CVect3
	void SetClmVect3(int i, int j, const CVect3 &v); // set i-row&j...(j+2)-column from CVect3
	void SetMat3(int i, int j, const CMat3 &m);	// set i...(i+2)-row&j...(j+2)-comumn from CMat3
	void ZeroRow(int i);						// set i-row to 0
	void ZeroClm(int j);						// set j-column to 0
	friend CMat operator~(const CMat &m);		// matirx transposition
	friend void symmetry(CMat &m);				// matirx symmetrization
	friend double MaxAbs(CMat &m);				// max abs element
	friend CVect diag(const CMat &m);			// diagonal of a matrix
	friend CMat diag(const CVect &v);			// diagonal matrix
};

class PSINS_API CRAvar
{
public:
	#define RAMAX 10
	int nR0, Rmaxflag[RAMAX];
	double ts, R0[RAMAX], Rmax[RAMAX], Rmin[RAMAX], tstau[RAMAX], r0[RAMAX];

	CRAvar(int nR0, double ts);
	void setR0(double r0, ...);
	void setTau(double tau, ...);
	void setRmax(double rmax, ...);
	void setRmin(double rmin, ...);
	void Update(double r, ...);
	double operator()(int k);			// get element sqrt(R0(k))
};

class PSINS_API CEarth
{
public:
	double a, b;
	double f, e, e2, ep, ep2;
	double wie;

	double sl, sl2, sl4, cl, tl, RMh, RNh, clRNh, f_RMh, f_RNh, f_clRNh;
	CVect3 pos, vn, wnie, wnen, wnin, gn, gcc;

	CEarth(double a0=glv.Re, double f0=glv.f, double g0=glv.g0);
	void Update(const CVect3 &pos, const CVect3 &vn=CVect3(0.0));
	CVect3 vn2dpos(const CVect3 &vn, double ts=1.0);
};

class PSINS_API CIMU
{
public:
	int nSamples, prefirst;
	CVect3 phim, dvbm, wm_1, vm_1;

	CIMU(void);
	void Update(const CVect3 *wm, const CVect3 *vm, int nSamples);
};

class PSINS_API CSINS
{
public:
	double nts, tk;
	CEarth eth;
	CIMU imu;
	CQuat qnb;
	CMat3 Cnb, Cnb0, Cbn, Kg, Ka;
	CVect3 wib, fb, fn, an, web, wnb, att, vn, vb, pos, eb, db, ebMax, dbMax;

	CSINS(const CQuat &qnb0=qI, const CVect3 &vn0=O31, const CVect3 &pos0=O31);    // initialization using quat attitude, velocity & position
	void Update(const CVect3 *wm, const CVect3 *vm, int nSamples, double ts);		// SINS update using Gyro&Acc samples
	void etm(CMat3 &Maa, CMat3 &Mav, CMat3 &Map, 
		CMat3 &Mva, CMat3 &Mvv, CMat3 &Mvp, CMat3 &Mpv, CMat3 &Mpp);	// SINS error transform matrix coefficients
};

class PSINS_API CAligni0
{
public:
	double tk;
	int t0, t1, t2;
	CVect3 wmm, vmm, vib0, vi0, Pib01, Pib02, Pi01, Pi02, tmpPib0, tmpPi0;
	CQuat qib0b;
	CEarth eth;
	CIMU imu;

	CAligni0(const CVect3 &pos=O31);
	CQuat Update(const CVect3 *wm, const CVect3 *vm, int nSamples, double ts);
};

class PSINS_API CKalman
{
public:
	double tk;
	int nq, nr, measflag;
	CMat Ft, Pk, Hk;
	CVect Xk, Zk, Qt, Rk, Pmax, Pmin;

	CKalman(int nq0, int nr0);
	void SetFt(CSINS &sins);			// process matrix set
	void SetHk(void);					// measurement matrix set
	void TimeUpdate(double ts);			// time update
	void MeasUpdate(double fading=1.0);	// measurement update
	void SetMeasFlag(int flag);			// measurement flag set
	void PkConstrain(void);				// Pk constrain: Pmin<diag(Pk)<Pmax
	void Feedback(CSINS &sins, double tauphi=INF, double taudvn=INF, double taudpos=INF, double taueb=INF, double taudb=INF);
};

class PSINS_API CTDKF:public CKalman  // Time-Distributed KF
{
public:
	double nts;
	int iter, ifn;
	CMat Fk, Pk1; 
	CVect Pxz, Qk, Kk, Hi, tmeas;
	CVect3 meanfn;

	CTDKF(int nq0, int nr0);
	int TDUpdate(CSINS &sins, double ts, int nStep=1);
	virtual void MeasRearrange(CSINS &sins) {};
};

class PSINS_API CIMUFile
{
	FILE *f;
public:
	CVect3 att0, vn0, pos0, gf, af;
	double t0, ts, t, g0;

	CIMUFile(char *fname);
	int skip(int n);
	int load(CVect3 *wm, CVect3 *vm, int n);
	~CIMUFile();
};

class PSINS_API CFileRdWt
{
	FILE *f;
public:
	double buff[128];
	int columns;

	CFileRdWt(char *fname, int columns);
	int load(int lines=1);
	CFileRdWt& operator<<(double d);
	CFileRdWt& operator<<(const CVect3 &v);
	CFileRdWt& operator<<(const CVect &v);
	CFileRdWt& operator<<(const CMat &m);
	CFileRdWt& operator<<(const CRAvar &R);
	CFileRdWt& operator<<(const CSINS &sins);
	CFileRdWt& operator<<(const CKalman &kf);
	CFileRdWt& operator>>(double &d);
	CFileRdWt& operator>>(CVect3 &v);
	CFileRdWt& operator>>(CVect &v);
	CFileRdWt& operator>>(CMat &m);
	~CFileRdWt();
};

CVect3 AlignCoarse(CVect3 wmm, CVect3 vmm, double latitude);

#endif
