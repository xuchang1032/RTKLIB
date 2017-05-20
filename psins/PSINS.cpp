#include "PSINS.h"

const CVect3 O31(0,0,0);
const CQuat qI(1.0,0,0,0);
const CMat3 I33(1,0,0, 0,1,0, 0,0,1);
const CGLV glv;

/***************************  class CGLV  *********************************/
CGLV::CGLV(double Re, double f, double wie0, double g0)
{
	this->Re = Re; this->f = f; this->wie = wie0; this->g0 = g0;
	e = sqrt(2*f-f*f); e2 = e*e;
    mg = 1.0e-3*g0;
    ug = 1.0e-6*glv.g0;
    deg = PI/180.0;
    min = deg/60.0;
    sec = min/60.0;
    ppm = 1.0e-6;
    hur = 3600.0;
	dps = deg/1.0;
    dph = deg/hur;
    dpsh = deg/sqrt(hur);
    dphpsh = dph/sqrt(hur);
    ugpsHz = ug/sqrt(1.0);
    ugpsh = ug/sqrt(hur);
    mpsh = 1/sqrt(hur); 
    mpspsh = 1/1/sqrt(hur);
    ppmpsh = ppm/sqrt(hur);
    secpsh = sec/sqrt(hur);
}

/***************************  class CVect3  *********************************/
CVect3::CVect3(void)
{
}

CVect3::CVect3(double xx, double yy, double zz)
{
	i=xx, j=yy, k=zz;
}

CVect3::CVect3(const double *pdata)
{
	i=*pdata++, j=*pdata++, k=*pdata++;
}

CVect3::CVect3(const CMat3 &Cnb)
{
	i=asinEx(Cnb.e21), j=atan2Ex(-Cnb.e20, Cnb.e22), k=atan2Ex(-Cnb.e01, Cnb.e11);
}

CVect3::CVect3(const CQuat &qnb)
{
	CMat3 Cnb(qnb);
	*this = CVect3(Cnb);
}

BOOL CVect3::IsZero(double eps)
{
	return (i<eps&&i>-eps && j<eps&&j>-eps && k<eps&&k>-eps);
}

BOOL CVect3::IsZeroXY(double eps)
{
	return (i<eps&&i>-eps && j<eps&&j>-eps);
}

CVect3 CVect3::operator+(const CVect3 &v)
{
	return CVect3(this->i+v.i, this->j+v.j, this->k+v.k);
}

CVect3 CVect3::operator-(const CVect3 &v)
{
	return CVect3(this->i-v.i, this->j-v.j, this->k-v.k);
}

CVect3 CVect3::operator*(const CVect3 &v)
{
	return CVect3(this->j*v.k-this->k*v.j, this->k*v.i-this->i*v.k, this->i*v.j-this->j*v.i);
}
	
CVect3 CVect3::operator*(double f)
{
	return CVect3(i*f, j*f, k*f);
}
	
CVect3 CVect3::operator/(double f)
{
	return CVect3(i/f, j/f, k/f);
}

CVect3& CVect3::operator+=(const CVect3 &v)
{ 
	i += v.i, j += v.j, k += v.k;
	return *this;
}

CVect3& CVect3::operator-=(const CVect3 &v)
{ 
	i -= v.i, j -= v.j, k -= v.k;
	return *this;
}

CVect3& CVect3::operator*=(double f)
{ 
	i *= f, j *= f, k *= f;
	return *this;
}

CVect3& CVect3::operator/=(double f)
{ 
	i /= f, j /= f, k /= f;
	return *this;
}

CVect3 operator*(double f, const CVect3 &v)
{
	return CVect3(v.i*f, v.j*f, v.k*f);
}
	
CVect3 operator-(const CVect3 &v)
{
	return CVect3(-v.i, -v.j, -v.k);
}

double norm(const CVect3 &v)
{
	return sqrt(v.i*v.i + v.j*v.j + v.k*v.k);
}

double normXY(const CVect3 &v)
{
	return sqrt(v.i*v.i + v.j*v.j);
}

double dot(const CVect3 &v1, const CVect3 &v2)
{
	return (v1.i*v2.i + v1.j*v2.j + v1.k*v2.k);
}

CQuat rv2q(const CVect3 &rv)
{
#define F1	(   2 * 1)		// define: Fk=2^k*k! 
#define F2	(F1*2 * 2)
#define F3	(F2*2 * 3)
#define F4	(F3*2 * 4)
#define F5	(F3*2 * 5)
	double n2 = rv.i*rv.i+rv.j*rv.j+rv.k*rv.k, c, f;
	if(n2<(PI/180.0*PI/180.0))	// 0.017^2 
	{
		double n4=n2*n2;
		c = 1.0 - n2*(1.0/F2) + n4*(1.0/F4);
		f = 0.5 - n2*(1.0/F3) + n4*(1.0/F5);
	}
	else
	{
		double n_2 = sqrt(n2)/2.0;
		c = cos(n_2);
		f = sin(n_2)/n_2*0.5;
	}
	return CQuat(c, f*rv.i, f*rv.j, f*rv.k);
}

CMat3 askew(CVect3 &v)
{
	return CMat3(0,  -v.k, v.j, 
				 v.k, 0.0,  -v.i,
				-v.j, v.i, 0);
}

CMat3 pos2Cen(const CVect3 &pos)
{
	double si = sin(pos.i), ci = cos(pos.i), sj = sin(pos.j), cj = cos(pos.j);
	return CMat3(	-sj, -si*cj,  ci*cj,  
					 cj, -si*sj,  ci*sj,  
					 0,   ci,     si      );	//Cen
}

CVect3 pp2vn(CVect3 &pos1, CVect3 &pos0, double ts, CEarth *pEth)
{
	double sl, cl, sl2, sq, sq2, RMh, RNh, clRNh;
	if(pEth)
	{
		RMh = pEth->RMh; clRNh = pEth->clRNh;
	}
	else
	{
		sl=sin(pos0.i); cl=cos(pos0.i); sl2=sl*sl;
		sq = 1-glv.e2*sl2; sq2 = sqrt(sq);
		RMh = glv.Re*(1-glv.e2)/sq/sq2+pos0.k;
		RNh = glv.Re/sq2+pos0.k;    clRNh = cl*RNh;
	}
    CVect3 vn = pos1 - pos0;
    return CVect3(vn.j*clRNh/ts, vn.i*RMh/ts, vn.k/ts);
}

/***************************  class CQuat  *********************************/
CQuat::CQuat(void)
{
}

CQuat::CQuat(double qq0, double qq1, double qq2, double qq3)
{
	q0=qq0, q1=qq1, q2=qq2, q3=qq3;
}

CQuat::CQuat(const double *pdata)
{
	q0=*pdata++, q1=*pdata++, q2=*pdata++, q3=*pdata++;
}

CQuat::CQuat(const CVect3 &att)
{
	CMat3 Cnb(att);
	*this = CQuat(Cnb);
}

CQuat::CQuat(const CMat3 &Cnb)
{
	double tmp;

	tmp = 1.0 + Cnb.e00 - Cnb.e11 - Cnb.e22;
	q1 = sqrt(fabs(tmp))/2.0;
	tmp = 1.0 - Cnb.e00 + Cnb.e11 - Cnb.e22;
	q2 = sqrt(fabs(tmp))/2.0;
	tmp = 1.0 - Cnb.e00 - Cnb.e11 + Cnb.e22;
	q3 = sqrt(fabs(tmp))/2.0;
	tmp = 1.0 - q1*q1 - q2*q2 - q3*q3;
	q0 = sqrt(fabs(tmp));

	if(Cnb.e21 - Cnb.e12 < 0)	/* sign decision */
	{
		q1 = -q1;
	}
	if(Cnb.e02 - Cnb.e20 < 0)
	{
		q2 = -q2;
	}
	if(Cnb.e10 - Cnb.e01 < 0)
	{
		q3 = -q3;
	}

	double nq = sqrt(q0*q0+q1*q1+q2*q2+q3*q3);
	q0 /= nq; q1 /= nq; q2 /= nq; q3 /= nq; 
}

CQuat CQuat::operator+(const CVect3 &phi)
{
	CQuat qtmp = rv2q(-phi);
	return qtmp*(*this);
}

CQuat CQuat::operator-(const CVect3 &phi)
{
	CQuat qtmp = rv2q(phi);
	return qtmp*(*this);
}

CVect3 CQuat::operator-(CQuat &quat)
{
	CQuat dq;
	
	dq = quat*(~(*this));
	if(dq.q0<0)
	{
		dq.q0=-dq.q0, dq.q1=-dq.q1, dq.q2=-dq.q2, dq.q3=-dq.q3;
	}
	double n2 = acos(dq.q0), f;
	if( sign(n2)!=0 )
	{
		f = 2.0/(sin(n2)/n2);
	}
	else
	{
		f = 2.0;
	}
	return CVect3(dq.q1,dq.q2,dq.q3)*f;
}

CQuat CQuat::operator*(const CQuat &quat)
{
	CQuat qtmp;
	qtmp.q0 = q0*quat.q0 - q1*quat.q1 - q2*quat.q2 - q3*quat.q3;
	qtmp.q1 = q0*quat.q1 + q1*quat.q0 + q2*quat.q3 - q3*quat.q2;
	qtmp.q2 = q0*quat.q2 + q2*quat.q0 + q3*quat.q1 - q1*quat.q3;
	qtmp.q3 = q0*quat.q3 + q3*quat.q0 + q1*quat.q2 - q2*quat.q1;
	return qtmp;
}

CQuat& CQuat::operator*=(const CQuat &quat)
{
	return (*this=*this*quat);
}

CQuat& CQuat::operator-=(const CVect3 &phi)
{
	CQuat qtmp = rv2q(phi);
	return (*this=qtmp*(*this));
}

CQuat operator~(const CQuat &q)
{
	return CQuat(q.q0,-q.q1,-q.q2,-q.q3);
}

CVect3 CQuat::operator*(const CVect3 &v)
{
	CQuat qtmp;
	CVect3 vtmp;
	qtmp.q0 =         - q1*v.i - q2*v.j - q3*v.k;
	qtmp.q1 = q0*v.i           + q2*v.k - q3*v.j;
	qtmp.q2 = q0*v.j           + q3*v.i - q1*v.k;
	qtmp.q3 = q0*v.k           + q1*v.j - q2*v.i;
	vtmp.i = -qtmp.q0*q1 + qtmp.q1*q0 - qtmp.q2*q3 + qtmp.q3*q2;
	vtmp.j = -qtmp.q0*q2 + qtmp.q2*q0 - qtmp.q3*q1 + qtmp.q1*q3;
	vtmp.k = -qtmp.q0*q3 + qtmp.q3*q0 - qtmp.q1*q2 + qtmp.q2*q1;
	return vtmp;
}

void normlize(CQuat *q)
{
	double nq=sqrt(q->q0*q->q0+q->q1*q->q1+q->q2*q->q2+q->q3*q->q3);
	q->q0 /= nq, q->q1 /= nq, q->q2 /= nq, q->q3 /= nq;
}

CVect3 q2rv(const CQuat &q)
{
	CQuat dq;
	dq = q;
	if(dq.q0<0)  { dq.q0=-dq.q0, dq.q1=-dq.q1, dq.q2=-dq.q2, dq.q3=-dq.q3; }
	if(dq.q0>1.0) dq.q0=1.0;
	double n2 = acos(dq.q0), f;
	if(n2>1.0e-20)
	{
		f = 2.0/(sin(n2)/n2);
	}
	else
	{
		f = 2.0;
	}
	return CVect3(dq.q1,dq.q2,dq.q3)*f;
}

/***************************  class CMat3  *********************************/
CMat3::CMat3(void)
{
}

CMat3::CMat3(double xx, double xy, double xz, 
		  double yx, double yy, double yz,
		  double zx, double zy, double zz )
{
	e00=xx,e01=xy,e02=xz; e10=yx,e11=yy,e12=yz; e20=zx,e21=zy,e22=zz;
}

CMat3::CMat3(const CVect3 &att)
{
	double
		si = sin(att.i), ci = cos(att.i),
		sj = sin(att.j), cj = cos(att.j),
		sk = sin(att.k), ck = cos(att.k);
	e00 =  cj*ck - si*sj*sk;	e01 =  -ci*sk;	e02 = sj*ck + si*cj*sk;
	e10 =  cj*sk + si*sj*ck;	e11 =  ci*ck;	e12 = sj*sk - si*cj*ck;
	e20 = -ci*sj;				e21 =  si;		e22 = ci*cj;
}

CMat3::CMat3(const CQuat &qnb)
{
	double 
		q11 = qnb.q0*qnb.q0, q12 = qnb.q0*qnb.q1, q13 = qnb.q0*qnb.q2, q14 = qnb.q0*qnb.q3, 
		q22 = qnb.q1*qnb.q1, q23 = qnb.q1*qnb.q2, q24 = qnb.q1*qnb.q3,     
		q33 = qnb.q2*qnb.q2, q34 = qnb.q2*qnb.q3,  
		q44 = qnb.q3*qnb.q3;
    e00 = q11+q22-q33-q44,  e01 = 2*(q23-q14),     e02 = 2*(q24+q13),
	e10 = 2*(q23+q14),      e11 = q11-q22+q33-q44, e12 = 2*(q34-q12),
	e20 = 2*(q24-q13),      e21 = 2*(q34+q12),     e22 = q11-q22-q33+q44 ;
}

CMat3::CMat3(const CVect3 &v0, const CVect3 &v1, const CVect3 &v2)
{
	e00 = v0.i, e01 = v0.j, e02 = v0.k;
	e10 = v1.i, e11 = v1.j, e12 = v1.k;
	e20 = v2.i, e21 = v2.j, e22 = v2.k;
}

CMat3 dv2att(CVect3 &va1, const CVect3 &va2, CVect3 &vb1, const CVect3 &vb2)
{
	CVect3 a=va1*va2, b=vb1*vb2, aa=a*va1, bb=b*vb1;
	CMat3 Ma(va1/norm(va1),a/norm(a),aa/norm(aa)), Mb(vb1/norm(vb1),b/norm(b),bb/norm(bb));
	return (~Ma)*(Mb);  //Cab
}

CMat3 operator~(const CMat3 &m)
{
	return CMat3(m.e00,m.e10,m.e20, m.e01,m.e11,m.e21, m.e02,m.e12,m.e22);
}

CMat3 CMat3::operator*(const CMat3 &mat)
{
	CMat3 mtmp;
	mtmp.e00 = e00*mat.e00 + e01*mat.e10 + e02*mat.e20;
	mtmp.e01 = e00*mat.e01 + e01*mat.e11 + e02*mat.e21;
	mtmp.e02 = e00*mat.e02 + e01*mat.e12 + e02*mat.e22;
	mtmp.e10 = e10*mat.e00 + e11*mat.e10 + e12*mat.e20;
	mtmp.e11 = e10*mat.e01 + e11*mat.e11 + e12*mat.e21;
	mtmp.e12 = e10*mat.e02 + e11*mat.e12 + e12*mat.e22;
	mtmp.e20 = e20*mat.e00 + e21*mat.e10 + e22*mat.e20;
	mtmp.e21 = e20*mat.e01 + e21*mat.e11 + e22*mat.e21;
	mtmp.e22 = e20*mat.e02 + e21*mat.e12 + e22*mat.e22;
	return mtmp;
}

CMat3 CMat3::operator+(const CMat3 &mat)
{
	CMat3 mtmp;
	mtmp.e00 = e00 + mat.e00;  mtmp.e01 = e01 + mat.e01;  mtmp.e02 = e02 + mat.e02;  
	mtmp.e10 = e10 + mat.e10;  mtmp.e11 = e11 + mat.e11;  mtmp.e12 = e12 + mat.e12;  
	mtmp.e20 = e20 + mat.e20;  mtmp.e21 = e21 + mat.e21;  mtmp.e22 = e22 + mat.e22;  
	return mtmp;
}

CMat3 CMat3::operator-(const CMat3 &mat)
{
	CMat3 mtmp;
	mtmp.e00 = e00 - mat.e00;  mtmp.e01 = e01 - mat.e01;  mtmp.e02 = e02 - mat.e02;  
	mtmp.e10 = e10 - mat.e10;  mtmp.e11 = e11 - mat.e11;  mtmp.e12 = e12 - mat.e12;  
	mtmp.e20 = e20 - mat.e20;  mtmp.e21 = e21 - mat.e21;  mtmp.e22 = e22 - mat.e22;  
	return mtmp;
}

CMat3 CMat3::operator*(double f)
{
	return CMat3(e00*f,e01*f,e02*f, e10*f,e11*f,e12*f, e21*f,e20*f,e22*f);
}

CMat3 operator*(double f, const CMat3 &m)
{
	return CMat3(m.e00*f,m.e01*f,m.e02*f, m.e10*f,m.e11*f,m.e12*f, m.e20*f,m.e21*f,m.e22*f);
}

CVect3 CMat3::operator*(const CVect3 &v)
{
	return CVect3(e00*v.i+e01*v.j+e02*v.k,e10*v.i+e11*v.j+e12*v.k,e20*v.i+e21*v.j+e22*v.k);
}

double det(const CMat3 &m)
{
	return m.e00*(m.e11*m.e22-m.e12*m.e21) - m.e01*(m.e10*m.e22-m.e12*m.e20) + m.e02*(m.e10*m.e21-m.e11*m.e20);
}

CMat3 inv(const CMat3 &m)
{
	double nm;
	nm = m.e00*(m.e11*m.e22-m.e12*m.e21) - m.e01*(m.e10*m.e22-m.e12*m.e20) + m.e02*(m.e10*m.e21-m.e11*m.e20);
	CMat3 mtmp;
	mtmp.e00 =  (m.e11*m.e22-m.e12*m.e21)/nm;
	mtmp.e10 = -(m.e10*m.e22-m.e12*m.e20)/nm;
	mtmp.e20 =  (m.e10*m.e21-m.e11*m.e20)/nm;
	mtmp.e01 = -(m.e01*m.e22-m.e02*m.e21)/nm;
	mtmp.e11 =  (m.e00*m.e22-m.e02*m.e20)/nm;
	mtmp.e21 = -(m.e00*m.e21-m.e01*m.e20)/nm;
	mtmp.e02 =  (m.e01*m.e12-m.e02*m.e11)/nm;
	mtmp.e12 = -(m.e00*m.e12-m.e02*m.e10)/nm;
	mtmp.e22 =  (m.e00*m.e11-m.e01*m.e10)/nm;
	return mtmp;
}

CVect3 diag(const CMat3 &m)
{
	return CVect3(m.e00, m.e11, m.e22);
}

CMat3 diag(const CVect3 &v)
{
	return CMat3(v.i,0,0, 0,v.j,0, 0,0,v.k);
}

/***************************  class CMat  *********************************/
CMat::CMat(void)
{
}
	
CMat::CMat(int row0, int clm0)
{
	row=row0; clm=clm0; rc=row*clm;
}

CMat::CMat(int row0, int clm0, double f)
{
	row=row0; clm=clm0; rc=row*clm;
	for(double *pd=dd, *pEnd=&dd[rc]; pd<pEnd; pd++)  *pd = f;
}

CMat CMat::operator*(const CMat &m0)
{
	assert(this->clm==m0.row);
	CMat mtmp(this->row,m0.clm);
	int m=this->row, k=this->clm, n=m0.clm;
	double *p=mtmp.dd, *p1i=this->dd; const double *p2=m0.dd;
	for(int i=0; i<m; i++,p1i+=k)
	{
		for(int j=0; j<n; j++)
		{
			double f=0.0, *p1is=p1i, *p1isEnd=&p1i[k]; const double *p2sj=&p2[j];
			for(; p1is<p1isEnd; p1is++,p2sj+=n)	f += (*p1is) * (*p2sj);
			*p++ = f;
		}
	}
	return mtmp;
}

CVect CMat::operator*(const CVect &v)
{
	assert(this->clm==v.row);
	CVect vtmp(this->row);
	double *p=vtmp.dd, *pEnd=&vtmp.dd[vtmp.row], *p1ij=this->dd; const double *p2End=&v.dd[v.row];
	for(; p<pEnd; p++)
	{
		double f=0.0; const double *p2j=v.dd;
		for(; p2j<p2End; p1ij++,p2j++)	f += (*p1ij) * (*p2j);
		*p = f;
	}
	return vtmp;
}

CMat CMat::operator+(const CMat &m0)
{
	assert(row==m0.row&&clm==m0.clm);
	CMat mtmp(row,clm);
	double *p=mtmp.dd, *pEnd=&mtmp.dd[rc], *p1=this->dd; const double *p2=m0.dd;
	while(p<pEnd)
	{ *p++ = (*p1++) + (*p2++); } 
	return mtmp;
}

CMat& CMat::operator+=(const CVect &v)
{
	assert(row==v.row||clm==v.clm);
	int row1 = row+1;
	double *p=dd, *pEnd=&dd[rc];
	for(const double *p1=v.dd; p<pEnd; p+=row1, p1++)	*p += *p1;
	return *this;
}

CMat CMat::operator-(const CMat &m0)
{
	assert(row==m0.row&&clm==m0.clm);
	CMat mtmp(row,clm);
	double *p=mtmp.dd, *pEnd=&mtmp.dd[rc], *p1=this->dd; const double *p2=m0.dd;
	while(p<pEnd)
	{ *p++ = (*p1++) - (*p2++); } 
	return mtmp;
}

CMat CMat::operator*(double f)
{
	CMat mtmp(row,clm);
	double *p=mtmp.dd, *pEnd=&mtmp.dd[rc], *p1=this->dd;
	while(p<pEnd)
	{ *p++ = (*p1++) * f; } 
	return mtmp;
}

CMat& CMat::operator+=(const CMat &m0)
{
	assert(row==m0.row&&clm==m0.clm);
	double *p=dd, *pEnd=&dd[rc]; const double *p1=m0.dd;
	while(p<pEnd)
	{ *p++ += *p1++; } 
	return *this;
}

CMat& CMat::operator-=(const CMat &m0)
{
	assert(row==m0.row&&clm==m0.clm);
	double *p=dd, *pEnd=&dd[rc]; const double *p1=m0.dd;
	while(p<pEnd)
	{ *p++ -= *p1++; } 
	return *this;
}

CMat& CMat::operator*=(double f)
{
	double *p=dd, *pEnd=&dd[rc];
	while(p<pEnd)
	{ *p++ *= f; } 
	return *this;
}

CMat& CMat::operator++()
{
	int row1=row+1;
	for(double *p=dd, *pEnd=&dd[rc]; p<pEnd; p+=row1)	*p += 1.0;
	return *this;
}

CMat operator~(const CMat &m0)
{
	CMat mtmp(m0.clm,m0.row);
	const double *pm=m0.dd;
	for(int i=0; i<m0.row; i++)
	{ for(int j=i; j<m0.rc; j+=m0.row) mtmp.dd[j] = *pm++; }
	return mtmp;
}

void symmetry(CMat &m)
{
	assert(m.row==m.clm);
	for(int i=0; i<m.clm; i++)
	{
		double *prow=&m.dd[i*m.clm+i+1], *prowEnd=&m.dd[i*m.clm+m.clm], *pclm=&m.dd[i*m.clm+i+m.clm];
		for(; prow<prowEnd; prow++,pclm+=m.clm)  *prow=*pclm=(*prow+*pclm)*0.5;
	}
}

double& CMat::operator()(int r, int c)
{
	return this->dd[r*this->clm+c];
}

void CMat::SetRow(int i, const CVect &v)
{
	assert(clm==v.clm);
	const double *p=v.dd;
	for(double *p1=&dd[i*clm],*pEnd=p1+clm; p1<pEnd; p++,p1++) *p1 = *p;
	return;
}

void CMat::SetClm(int j, const CVect &v)
{
	assert(row==v.row);
	const double *p=v.dd;
	for(double *p1=&dd[j],*pEnd=&dd[rc]; p1<pEnd; p++,p1+=clm) *p1 = *p;
	return;
}

void CMat::SetRowVect3(int i, int j, const CVect3 &v)
{
	*(CVect3*)&dd[i*clm+j] = v;
}

void CMat::SetClmVect3(int i, int j, const CVect3 &v)
{
	double *p=&dd[i*clm+j];
	*p = v.i; p += clm;
	*p = v.j; p += clm;
	*p = v.k;
}

void CMat::SetMat3(int i, int j, const CMat3 &m)
{
	double *p=&dd[i*clm+j];
	*(CVect3*)p = *(CVect3*)&m.e00;  p += clm;
	*(CVect3*)p = *(CVect3*)&m.e10;  p += clm;
	*(CVect3*)p = *(CVect3*)&m.e20;
}

CVect CMat::GetRow(int i)
{
	CVect v;
	v.row=1; v.clm=clm;
	for(double *p=v.dd,*p1=&dd[i*clm],*pEnd=p1+clm; p1<pEnd; p++,p1++) *p = *p1;
	return v;
}

CVect CMat::GetClm(int j)
{
	CVect v;
	v.row=row; v.clm=1;
	for(double *p=v.dd,*p1=&dd[j],*pEnd=&dd[rc]; p1<pEnd; p++,p1+=clm) *p = *p1;
	return v;
}

void CMat::ZeroRow(int i)
{
	for(double *p=&dd[i*clm],*pEnd=p+clm; p<pEnd; p++) *p = 0.0;
	return;
}

void CMat::ZeroClm(int j)
{
	for(double *p=&dd[j],*pEnd=&dd[rc]; p<pEnd; p+=clm) *p = 0.0;
	return;
}

void CMat::SetDiag(double f, ...)
{
	va_list vl;
	va_start(vl, f);
	double *p=dd, *pEnd=&dd[rc];
	for(int row1=row+1; p<pEnd; p+=row1)
	{ *p = f;  f = va_arg(vl, double);	}
	va_end(vl);
}

void CMat::SetDiag2(double f, ...)
{
	va_list vl;
	va_start(vl, f);
	double *p=dd, *pEnd=&dd[rc];
	for(int row1=row+1; p<pEnd; p+=row1)
	{ *p = f*f;  f = va_arg(vl, double);	}
	va_end(vl);
}

double MaxAbs(CMat &m)
{
	double positive=0.0, negative=0.0;
	for(double *p=m.dd,*pEnd=&m.dd[m.rc]; p<pEnd; p++)
	{
		if(*p>positive)	 positive = *p;
		else if(*p<negative)  negative = *p;
	}
	return positive>-negative ? positive : -negative;
}

CVect diag(const CMat &m)
{
	int row1 = m.row+1;
	CVect vtmp(m.row,1);
	double *p=vtmp.dd, *pEnd=&vtmp.dd[vtmp.row];
	for(const double *p1=m.dd; p<pEnd; p++, p1+=row1)	*p = *p1;
	return vtmp;
}

CMat diag(const CVect &v)
{
	int rc = v.row>v.clm ? v.row : v.clm, rc1=rc+1;
	CMat mtmp(rc,rc,0.0);
	double *p=mtmp.dd;
	for(const double *p1=v.dd, *p1End=&v.dd[rc]; p1<p1End; p+=rc1, p1++)	*p = *p1;
	return mtmp;
}

/***************************  class CVect  *********************************/
CVect::CVect(void)
{
}

CVect::CVect(int row0, int clm0)
{
	if(clm0==1) { row=row0; clm=1;   }
	else		{ row=1;    clm=clm0;}
 }

CVect::CVect(int row0, double f)
{
	row=row0; clm=1;
	for(int i=0;i<row;i++) dd[i]=f;
 }

CVect::CVect(const CVect3 &v)
{
	row=3; clm=1;
	dd[0]=v.i; dd[1]=v.j; dd[2]=v.k;
 }

CVect::CVect(const CVect3 &v1, const CVect3 v2)
{
	row=6; clm=1;
	dd[0]=v1.i; dd[1]=v1.j; dd[2]=v1.k;
	dd[3]=v2.i; dd[4]=v2.j; dd[5]=v2.k;
}

CVect operator~(const CVect &v)
{
	CVect vtmp=v;
	vtmp.row=v.clm; vtmp.clm=v.row;
	return vtmp;
}

CVect CVect::operator*(const CMat &m)
{
	assert(clm==m.row);
	CVect vtmp(row,clm);
	double *p=vtmp.dd, *p1End=&dd[clm];
	for(int j=0; j<clm; p++,j++)
	{
		double f=0.0, *p1j=dd; const double *p2jk=&m.dd[j];
		for(; p1j<p1End; p1j++,p2jk+=m.clm)	 f += (*p1j) * (*p2jk);
		*p = f;
	}
	return vtmp;
}

CMat CVect::operator*(const CVect &v)
{
	assert(clm==v.row);
	CMat mtmp(row,v.clm);
	if(row==1 && v.clm==1)  // (1x1) = (1xn)*(nx1)
	{
		double f = 0.0;
		for(int i=0; i<clm; i++)  f += dd[i]*v.dd[i];
		mtmp.dd[0] = f;
	}
	else    // (nxn) = (nx1)*(1xn)
	{
		double *p = mtmp.dd;
		for(int i=0; i<row; i++)
		{
			for(int j=0; j<v.clm; j++)  *p++ = dd[i]*v.dd[j];
		}
	}
	return mtmp;
}

CVect CVect::operator+(const CVect &v)
{
	assert(row==v.row&&clm==v.clm);
	const double *p2=v.dd;
	CVect vtmp(row,clm);
	for(double *p=vtmp.dd,*p1=dd,*p1End=&dd[row>clm?row:clm]; p1<p1End; p++,p1++,p2++)  { *p=*p1+*p2; }
	return vtmp;
}

CVect CVect::operator-(const CVect &v)
{
	assert(row==v.row&&clm==v.clm);
	const double *p2=v.dd;
	CVect vtmp(row,clm);
	for(double *p=vtmp.dd,*p1=dd,*p1End=&dd[row>clm?row:clm]; p1<p1End; p++,p1++,p2++)  { *p=*p1-*p2; }
	return vtmp;
}
	
CVect CVect::operator*(double f)
{
	CVect vtmp(row,clm);
	for(double *p=vtmp.dd,*p1=dd,*p1End=&dd[row>clm?row:clm]; p1<p1End; p++,p1++)  { *p=*p1*f; }
	return vtmp;
}

CVect& CVect::operator+=(const CVect &v)
{
	assert(row==v.row&&clm==v.clm);
	const double *p1 = v.dd;
	for(double *p=dd, *pEnd=&dd[row>clm?row:clm]; p<pEnd; p++,p1++)  { *p += *p1; }
	return *this;
}

CVect& CVect::operator-=(const CVect &v)
{
	assert(row==v.row&&clm==v.clm);
	const double *p1 = v.dd;
	for(double *p=dd, *pEnd=&dd[row>clm?row:clm]; p<pEnd; p++,p1++)  { *p -= *p1; }
	return *this;
}

CVect& CVect::operator*=(double f)
{
	for(double *p=dd, *pEnd=&dd[row>clm?row:clm]; p<pEnd; p++)  { *p *= f; }
	return *this;
}

double norm(const CVect &v)
{
	const double *p=v.dd, *pEnd=&v.dd[v.row>v.clm?v.row:v.clm];
	double f=0.0;
	for(; p<pEnd; p++)  { f += (*p)*(*p); }
	return sqrt(f);
}

double& CVect::operator()(int r)
{
	return this->dd[r];
}

void CVect::Set(double f, ...)
{
	assert(row<=MMD&&clm<=MMD);
	va_list vl;
	va_start(vl, f);
	for(int i=0, rc=row>clm?row:clm; i<rc; i++)
	{ dd[i] = f;  f = va_arg(vl, double);	}
	va_end(vl);
}

void CVect::Set2(double f, ...)
{
	assert(row<=MMD&&clm<=MMD);
	va_list vl;
	va_start(vl, f);
	for(int i=0, rc=row>clm?row:clm; i<rc; i++)
	{ dd[i] = f*f;  f = va_arg(vl, double);	}
	va_end(vl);
}

/***************************  class CRAvar  *********************************/
CRAvar::CRAvar(int nR0, double ts)
{
	assert(nR0<RAMAX);
	this->nR0 = nR0;
	this->ts = ts;
}

void CRAvar::setR0(double r0, ...)
{
	va_list vl;
	va_start(vl, r0);
	for(int i=0; i<nR0; i++)
	{
		R0[i] = r0*r0;  Rmax[i] = 100.0*R0[i];  Rmin[i] = 0.01*R0[i];
		r0 = va_arg(vl, double);
		this->r0[i] = 0.0;  Rmaxflag[i] = 1;
	}
	va_end(vl);
}

void CRAvar::setTau(double tau, ...)
{
	va_list vl;
	va_start(vl, tau);
	for(int i=0; i<nR0; i++)
	{
		tstau[i] = ts>tau ? 1.0 : ts/tau;
		tau = va_arg(vl, double);
	}
	va_end(vl);
}

void CRAvar::setRmax(double rmax, ...)
{
	va_list vl;
	va_start(vl, rmax);
	for(int i=0; i<nR0; i++)
	{
		double rmax2 = rmax*rmax;
		if(rmax2>R0[i]) Rmax[i] = rmax2;
		rmax = va_arg(vl, double);
	}
	va_end(vl);
}

void CRAvar::setRmin(double rmin, ...)
{
	va_list vl;
	va_start(vl, rmin);
	for(int i=0; i<nR0; i++)
	{
		double rmin2 = rmin*rmin;
		if(rmin2>EPS&&rmin2<R0[i]) Rmin[i] = rmin2;
		rmin = va_arg(vl, double);
	}
	va_end(vl);
}

void CRAvar::Update(double r, ...)
{
	va_list vl;
	va_start(vl, r);
	for(int i=0; i<nR0; i++)
	{
		double dr2=r-r0[i]; dr2=dr2*dr2; r0[i]=r;
		if(dr2>R0[i]) R0[i]=dr2; else R0[i]=(1.0-tstau[i])*R0[i]+tstau[i]*dr2;
		if(R0[i]<Rmin[i]) R0[i]=Rmin[i];
		if(R0[i]>Rmax[i]) {R0[i]=Rmax[i];Rmaxflag[i]=1;} else {Rmaxflag[i]=0;}
		r = va_arg(vl, double);
	}
	va_end(vl);
}

double CRAvar::operator()(int k)
{
	return Rmaxflag[k] ? INF : sqrt(R0[k]);
}

/***************************  class CKalman  *********************************/
CKalman::CKalman(int nq0, int nr0)
{
	assert(nq0<=MMD&&nr0<=MMD);
	tk = 0.0;
	nq = nq0; nr = nr0;
	Ft = Pk = CMat(nq,nq,0.0);
	Hk = CMat(nr,nq,0.0);
	Qt = Pmin = Xk = CVect(nq,0.0);  Pmax = CVect(nq,INF);
	Rk = Zk = CVect(nr,0.0);
	measflag = 0;
}

void CKalman::TimeUpdate(double ts)
{
	tk += ts;
	CMat Fk=++(Ft*ts);  // Fk = I+Ft*ts
	Xk = Fk * Xk;
	Pk = Fk*Pk*(~Fk), Pk += Qt*ts;
}

void CKalman::SetMeasFlag(int flag)
{
	measflag = (flag==0) ? 0 : (measflag|flag);
}

void CKalman::MeasUpdate(double fading)
{
	CVect Pxz, Kk, Hi;
	for(int i=0; i<nr; i++)
	{
		if(measflag&(0x01<<i))
		{
			Hi = Hk.GetRow(i);
			Pxz = Pk*(~Hi);
			double Pzz = (Hi*Pxz)(0,0) + Rk(i);
			Kk = Pxz*(1.0/Pzz);
			Xk += Kk*(Zk(i)-(Hi*Xk)(0,0));
			Pk -= Kk*(~Pxz);
		}
	}
	if(fading>1.0) Pk *= fading;
	symmetry(Pk);
}

void CKalman::PkConstrain(void)
{
	int i=0, nq1=nq+1;
	for(double *p=Pk.dd,*pmin=Pmin.dd,*pminEnd=&Pmin.dd[nq],*pmax=Pmax.dd; pmin<pminEnd; p+=nq1,pmin++,pmax++)
	{
		if(*p<*pmin && *p>EPS)
		{
			*p = *pmin;
		}
		else if(*p>*pmax)
		{
			double sqf=sqrt(*pmax/(*p));
			for(double *prow=&Pk.dd[i*Pk.clm],*prowEnd=prow+nq,*pclm=&Pk.dd[i]; prow<prowEnd; prow++,pclm+=nq)
			{
				*prow *= sqf;
				*pclm *= sqf;
			}
		}
		i++;
	}
}

void CKalman::Feedback(CSINS &sins, double tauphi, double taudvn, double taudpos, double taueb, double taudb)
{
	double afa, hur10=36000;
	if(tauphi<hur10)
	{
		afa = sins.nts<tauphi ? sins.nts/tauphi : 1.0;
		CVect3 *phi = (CVect3*)&Xk.dd[0];
//		sins.qnb = sins.qnb - (*phi)*afa;  *phi = (*phi)*(1-afa);
		sins.qnb -= CVect3(phi->i*afa,phi->j*afa,phi->k*0.1*afa); *phi = CVect3(phi->i*(1-afa),phi->j*(1-afa),phi->k*(1-0.1*afa));
	}
	if(taudvn<hur10)
	{
		afa = sins.nts<taudvn ? sins.nts/taudvn : 1.0;
		CVect3 *dvn = (CVect3*)&Xk.dd[3];
		sins.vn -= (*dvn)*afa;    *dvn = (*dvn)*(1-afa);
	}
	if(taudpos<hur10)
	{
		afa = sins.nts<taudpos ? sins.nts/taudpos : 1.0;
		CVect3 *dpos = (CVect3*)&Xk.dd[6];
		sins.pos -= (*dpos)*afa;    *dpos = (*dpos)*(1-afa);
	}
	if(taueb<hur10)
	{
		afa = sins.nts<taueb ? sins.nts/taueb : 1.0;
		int i=0;
		for(double *peb=&sins.eb.i, *pMax=&sins.ebMax.i, *p=&Xk.dd[9]; i<3; i++, peb++, pMax++, p++)
			if((*peb>*pMax&&*p<0) || (*peb<-*pMax&&*p>0) || (*peb>=-*pMax&&*peb<=*pMax)) 
				*peb += (*p)*afa,  *p = (*p)*(1-afa);
	}
	if(taudb<hur10)
	{
		afa = sins.nts<taudb ? sins.nts/taudb : 1.0;
		int i=0;
		for(double *pdb=&sins.db.i, *pMax=&sins.dbMax.i, *p=&Xk.dd[12]; i<3; i++, pdb++, pMax++, p++)
			if((*pdb>*pMax&&*p<0) || (*pdb<-*pMax&&*p>0) || (*pdb>=-*pMax&&*pdb<=*pMax)) 
				*pdb += (*p)*afa,  *p = (*p)*(1-afa);
	}
}

void CKalman::SetFt(CSINS &sins)
{
	CMat3 Maa, Mav, Map, Mva, Mvv, Mvp, Mpv, Mpp, Cnb;
	sins.etm(Maa, Mav, Map, Mva, Mvv, Mvp, Mpv, Mpp);
	Cnb = sins.Cnb;
//	Ft = [ Maa    Mav    Map    -ins.Cnb  O33 
//         Mva    Mvv    Mvp     O33      ins.Cnb 
//         O33    Mpv    Mpp     O33      O33
//         zeros(6,9)  diag(-1./[ins.tauG;ins.tauA]) ];
	// phi
	                 ; Ft(0,1) = Maa.e01; Ft(0,2) = Maa.e02;                     ;  Ft(0,4) = Mav.e01;
	Ft(1,0) = Maa.e10;                  ; Ft(1,2) = Maa.e12;    Ft(1,3) = Mav.e10;
	Ft(2,0) = Maa.e20; Ft(2,1) = Maa.e21;                  ;    Ft(2,3) = Mav.e20;
	                 ;                  ; Ft(0,8) = Map.e02;    Ft(0,9) =-Cnb.e00;  Ft(0,10) =-Cnb.e01; Ft(0,11) =-Cnb.e02; 
	Ft(1,6) = Map.e10;                  ; Ft(1,8) = Map.e12;    Ft(1,9) =-Cnb.e10;  Ft(1,10) =-Cnb.e11; Ft(1,11) =-Cnb.e12; 
	Ft(2,6) = Map.e20;                  ; Ft(2,8) = Map.e22;    Ft(2,9) =-Cnb.e20;  Ft(2,10) =-Cnb.e21; Ft(2,11) =-Cnb.e22; 
	// dv
	                 ; Ft(3,1) = Mva.e01; Ft(3,2) = Mva.e02;    Ft(3,3) = Mvv.e00;  Ft(3,4) = Mvv.e01;  Ft(3,5) = Mvv.e02; 
	Ft(4,0) = Mva.e10;                  ; Ft(4,2) = Mva.e12;    Ft(4,3) = Mvv.e10;  Ft(4,4) = Mvv.e11;  Ft(4,5) = Mvv.e12; 
	Ft(5,0) = Mva.e20; Ft(5,1) = Mva.e21;                  ;    Ft(5,3) = Mvv.e20;  Ft(5,4) = Mvv.e21; 
	Ft(3,6) = Mvp.e00;                  ; Ft(3,8) = Mvp.e02;    Ft(3,12) = Cnb.e00; Ft(3,13) = Cnb.e01; Ft(3,14) = Cnb.e02; 
	Ft(4,6) = Mvp.e10;                  ; Ft(4,8) = Mvp.e12;    Ft(4,12) = Cnb.e10; Ft(4,13) = Cnb.e11; Ft(4,14) = Cnb.e12; 
	Ft(5,6) = Mvp.e20;                  ; Ft(5,8) = Mvp.e22;    Ft(5,12) = Cnb.e20; Ft(5,13) = Cnb.e21; Ft(5,14) = Cnb.e22; 
	// dpos
	                 ; Ft(6,4) = Mpv.e01;                  ;                     ;                    ; Ft(6,8) = Mpp.e02; 
	Ft(7,3) = Mpv.e10;                  ;                  ;    Ft(7,6) = Mpp.e10;                    ; Ft(7,8) = Mpp.e12; 
	                 ;                  ; Ft(8,5) = Mpv.e22;                     ;                    ;                  ;
	// dKgz
	Ft(2,15) = -sins.wib.k*sins.Cnb.e22;
}

void CKalman::SetHk(void)
{
	Hk(0,3) = Hk(1,4) = Hk(2,5) = 1.0; Hk(3,6) = Hk(4,7) = Hk(5,8) = 1.0;
}

/***************************  class CEarth  *********************************/
CEarth::CEarth(double a0, double f0, double g0)
{
	a = a0;	f = f0; wie = glv.wie; 
	b = (1-f)*a;
	e = sqrt(a*a-b*b)/a;	e2 = e*e;
	gn = CVect3(0, 0, -g0);
}

void CEarth::Update(const CVect3 &pos, const CVect3 &vn)
{
	this->pos = pos;  this->vn = vn;
	sl = sin(pos.i), cl = cos(pos.i), tl = sl/cl;
	double sq = 1-e2*sl*sl, sq2 = sqrt(sq);
	RMh = a*(1-e2)/sq/sq2+pos.k;	f_RMh = 1.0/RMh;
	RNh = a/sq2+pos.k;    clRNh = cl*RNh;  f_RNh = 1.0/RNh; f_clRNh = 1.0/clRNh;
	wnie.i = 0,				wnie.j = wie*cl,		wnie.k = wie*sl;
	wnen.i = -vn.j*f_RMh,	wnen.j = vn.i*f_RNh,	wnen.k = wnen.j*tl;
	wnin = wnie + wnen;
	sl2 = sl*sl, sl4 = sl2*sl2;
	gn.k = -( glv.g0*(1+5.27094e-3*sl2+2.32718e-5*sl4)-3.086e-6*pos.k );
	gcc = gn - (wnie+wnin)*vn;
}

CVect3 CEarth::vn2dpos(const CVect3 &vn, double ts)
{
	return CVect3(vn.j*f_RMh, vn.i*f_clRNh, vn.k)*ts;
}

/***************************  class CIMU  *********************************/
CIMU::CIMU(void)
{
	prefirst = 1;
}

void CIMU::Update(const CVect3 *wm, const CVect3 *vm, int nSamples)
{
	static double conefactors[5][4] = {				// coning coefficients
		{2./3},										// 2
		{9./20, 27./20},							// 3
		{54./105, 92./105, 214./105},				// 4
		{250./504, 525./504, 650./504, 1375./504}	// 5
		};
	int i;
	double *pcf = conefactors[nSamples-2];
	CVect3 cm(0.0), sm(0.0), wmm(0.0), vmm(0.0);

	this->nSamples = nSamples;
	if(nSamples==1)  // one-plus-previous sample
	{
		if(prefirst==1) {wm_1=wm[0]; vm_1=vm[0]; prefirst=0;}
		cm = 1.0/12*wm_1; wm_1=wm[0]; 
		sm = 1.0/12*vm_1; vm_1=vm[0];
	}
	if(nSamples>1) prefirst=1;
	for(i=0; i<nSamples-1; i++)
	{
		cm += pcf[i]*wm[i];
		sm += pcf[i]*vm[i];
		wmm += wm[i];
		vmm += vm[i];
	}
	wmm += wm[i];
	vmm += vm[i];
	phim = wmm + cm*wm[i];
	dvbm = vmm + 1.0/2*wmm*vmm + (cm*vm[i]+sm*wm[i]);
}

/***************************  class CSINS  *********************************/
CSINS::CSINS(const CQuat &qnb0, const CVect3 &vn0, const CVect3 &pos0)
{
	tk = 0.0;
	qnb = qnb0;	vn = vn0, pos = pos0;
	eth.Update(pos0, vn0);
	Cnb = CMat3(qnb); att = CVect3(Cnb); Cnb0 = Cnb; Cbn = ~Cnb;
	Kg = Ka = I33; eb = db = O31; ebMax = dbMax = CVect3(1,1,1)*INF;
	wib = fb = fn = an = wnb = web = O31;
}

void CSINS::Update(const CVect3 *wm, const CVect3 *vm, int nSamples, double ts)
{
	nts = nSamples*ts;	tk += nts;
	double nts2 = nts/2;
	imu.Update(wm, vm, nSamples);
	imu.phim = Kg*imu.phim - eb*nts; imu.dvbm = Ka*imu.dvbm - db*nts;  // IMU calibration
	CVect3 vn01 = vn+an*nts2, pos01 = pos+eth.vn2dpos(vn01,nts2);
	eth.Update(pos01, vn01);
	wib = imu.phim/nts; fb = imu.dvbm/nts;
	web = wib - Cbn*eth.wnie;
	wnb = wib - (qnb*rv2q(imu.phim/2))*eth.wnin;
	fn = qnb*fb;
	an = rv2q(-eth.wnin*nts2)*fn+eth.gcc;
	CVect3 vn1 = vn + an*nts;
	pos = pos + eth.vn2dpos(vn+vn1, nts2);	vn = vn1;
	Cnb0 = Cnb;
	qnb = rv2q(-eth.wnin*nts)*qnb*rv2q(imu.phim);
	Cnb = CMat3(qnb); att = CVect3(Cnb); Cbn = ~Cnb; vb = Cbn*vn;
}

void CSINS::etm(CMat3 &Maa, CMat3 &Mav, CMat3 &Map, CMat3 &Mva, CMat3 &Mvv, CMat3 &Mvp, CMat3 &Mpv, CMat3 &Mpp)
{
	double tl=eth.tl, secl=1.0/eth.cl, secl2=secl*secl, 
		wN=eth.wnie.j, wU=eth.wnie.k, vE=vn.i, vN=vn.j;
	double f_RMh=eth.f_RMh, f_RNh=eth.f_RNh, f_clRNh=eth.f_clRNh, 
		f_RMh2=f_RMh*f_RMh, f_RNh2=f_RNh*f_RNh;
	CMat3 Avn=askew(vn),
		Mp1(0,0,0, -wU,0,0, wN,0,0),
		Mp2(0,0,vN*f_RMh2, 0,0,-vE*f_RNh2, vE*secl2*f_RNh,0,-vE*tl*f_RNh2);
	CVect3 _wnin = -eth.wnin; 	Maa = askew(_wnin);  // for Keil/VS2017 ???
//	Maa = askew(-eth.wnin);
	Mav = CMat3(0,-f_RMh,0, f_RNh,0,0, tl*f_RNh,0,0);
	Map = Mp1+Mp2;
	Mva = askew(fn);
	CVect3 wnien = eth.wnie+eth.wnin;
	Mvv = Avn*Mav - askew(wnien);
	Mvp = Avn*(Mp1+Map);
	double scl = eth.sl*eth.cl;
    Mvp.e20 = Mvp.e20-glv.g0*(5.27094e-3*2*scl+2.32718e-5*4*eth.sl2*scl); Mvp.e22 = Mvp.e22+3.086e-6;
	Mpv = CMat3(0,f_RMh,0, f_clRNh,0,0, 0,0,1);
	Mpp = CMat3(0,0,-vN*f_RMh2, vE*tl*f_clRNh,0,-vE*secl*f_RNh2, 0,0,0);
}

/*********************  class CIMUFile (PSINS Format File)  ************************/
CIMUFile::CIMUFile(char *fname)
{
	f = fopen(fname, "rt");
	// skip notation
	char prechar=0, curchar=0;
	while(1)
	{
		fscanf(f, "%c", &curchar);
		if(prechar=='\n'&&curchar!='%')	{ fseek(f,-1L,SEEK_CUR); break; }
		prechar = curchar;
	}
	// read data info
	fscanf(f, "%lf %lf %lf %lf %lf %lf", &att0.i, &att0.j, &att0.k, &vn0.i, &vn0.j, &vn0.k);
	fscanf(f, "%lf %lf %lf %lf %lf %lf", &pos0.i, &pos0.j, &pos0.k, &t0, &ts, &g0);
	fscanf(f, "%lf %lf %lf %lf %lf %lf", &gf.i, &gf.j, &gf.k, &af.i, &af.j, &af.k);
	att0 = att0*glv.deg; pos0.i *= glv.deg; pos0.j *= glv.deg; ts /= 1000.0;
	gf = gf*glv.sec; af = af*glv.ug;
	t = t0;
}

CIMUFile::~CIMUFile()
{
	if(f) fclose(f);
}

int CIMUFile::skip(int n)
{
	int gx, gy, gz, ax, ay, az;
	for(int i=0; i<n; i++)
	{
		fscanf(f, "%d %d %d %d %d %d", &gx, &gy, &gz, &ax, &ay, &az);
	}
	if(feof(f)) return 0;
	t = t + n*ts;
	return n;
}

int CIMUFile::load(CVect3 *wm, CVect3 *vm, int n)
{
	int gx, gy, gz, ax, ay, az;
	for(int i=0; i<n; i++)
	{
		fscanf(f, "%d %d %d %d %d %d", &gx, &gy, &gz, &ax, &ay, &az);
		wm[i].i=gx*gf.i; wm[i].j=gy*gf.j; wm[i].k=gz*gf.k; vm[i].i=ax*af.i; vm[i].j=ay*af.j; vm[i].k=az*af.k;
	}
	if(feof(f)) return 0;
	t = t + n*ts;
	if(fmod(t+0.01*ts,50.0)<n*ts)
		printf("%.4lf\n", t);
	return n;
}

/******************************  File Read or Write *********************************/
CFileRdWt::CFileRdWt(char *fname, int columns)
{
	if(columns>0)  // txt file read
	{
		f = fopen(fname, "rt");
		unsigned char cpre='\n', c;
		while(1)  // skip comments
		{
			fscanf(f, "%c", &c);
			if(c==' ') continue;
			if(cpre=='\n' && c!='%') break;
			cpre = c;
		}
		int res = fseek(f, -2L, SEEK_CUR);
		this->columns = columns;
	}
	else if(columns==0)  // bin file write
	{
		f = fopen(fname, "wb");
	}
	else if(columns==-1)  // bin file read
	{
		f = fopen(fname, "rb");
	}
}

int CFileRdWt::load(int lines)
{
	for(int i=0; i<lines; i++)
	{
		for(int j=0; j<columns; j++)	{ fscanf(f, "%lf,", &buff[j]); }		
//		for(int j=0; j<columns-1; j++)	{ fscanf(f, "%lf,", &buff[j]); }	fscanf(f, "%lf ", &buff[j]);	
		if(feof(f))  return 0;
	}
	return 1;
}

CFileRdWt& CFileRdWt::operator<<(double d)
{
	fwrite(&d, 1, sizeof(double), f);
	return *this;
}

CFileRdWt& CFileRdWt::operator<<(const CVect3 &v)
{
	fwrite(&v, 1, sizeof(v), f);
	return *this;
}

CFileRdWt& CFileRdWt::operator<<(const CVect &v)
{
	fwrite(v.dd, v.clm*v.row, sizeof(double), f);
	return *this;
}

CFileRdWt& CFileRdWt::operator<<(const CMat &m)
{
	fwrite(m.dd, m.clm*m.row, sizeof(double), f);
	return *this;
}

CFileRdWt& CFileRdWt::operator<<(const CRAvar &R)
{
	fwrite(R.R0, R.nR0, sizeof(double), f);
	return *this;
}

CFileRdWt& CFileRdWt::operator<<(const CSINS &sins)
{
	return *this<<sins.att<<sins.vn<<sins.pos<<sins.eb<<sins.db<<sins.tk;
}

CFileRdWt& CFileRdWt::operator<<(const CKalman &kf)
{
	return *this<<kf.Xk<<diag(kf.Pk)<<kf.tk;
}

CFileRdWt::~CFileRdWt()
{
	if(f) fclose(f); 
}

CFileRdWt& CFileRdWt::operator>>(double &d)
{
	fread(&d, 1, sizeof(double), f);
	return *this;
}

CFileRdWt& CFileRdWt::operator>>(CVect3 &v)
{
	fread(&v, 1, sizeof(v), f);
	return *this;
}

CFileRdWt& CFileRdWt::operator>>(CVect &v)
{
	fread(v.dd, v.clm*v.row, sizeof(double), f);
	return *this;
}

CFileRdWt& CFileRdWt::operator>>(CMat &m)
{
	fread(m.dd, m.clm*m.row, sizeof(double), f);
	return *this;
}

/***************************  function AlignCoarse  *********************************/
CVect3 AlignCoarse(CVect3 wmm, CVect3 vmm, double latitude)
{
	double T11, T12, T13, T21, T22, T23, T31, T32, T33;
	double cl = cos(latitude), tl = tan(latitude), nn;
	CVect3 wbib = wmm / norm(wmm),  fb = vmm / norm(vmm);
	T31 = fb.i,				 T32 = fb.j,			 	T33 = fb.k;
	T21 = wbib.i/cl-T31*tl,	 T22 = wbib.j/cl-T32*tl,	T23 = wbib.k/cl-T33*tl;		nn = sqrt(T21*T21+T22*T22+T23*T23);  T21 /= nn, T22 /= nn, T23 /= nn;
	T11 = T22*T33-T23*T32,	 T12 = T23*T31-T21*T33,		T13 = T21*T32-T22*T31;		nn = sqrt(T11*T11+T12*T12+T13*T13);  T11 /= nn, T12 /= nn, T13 /= nn;
	CMat3 Cnb(T11, T12, T13, T21, T22, T23, T31, T32, T33);
	return CVect3(Cnb);
}

CAligni0::CAligni0(const CVect3 &pos)
{
	eth.Update(pos);
	tk = 0;
	t0 = t1 = 10, t2 = 0; 
	wmm = vmm = vib0 = vi0 = Pib01 = Pib02 = Pi01 = Pi02 = O31;
	qib0b = CQuat(1.0);
}

CQuat CAligni0::Update(const CVect3 *wm, const CVect3 *vm, int nSamples, double ts)
{
	double nts = nSamples*ts;
	imu.Update(wm, vm, nSamples);
	wmm = wmm + imu.phim;  vmm = vmm + imu.dvbm;
	// vtmp = qib0b * (vm + 1/2 * wm X vm)
	CVect3 vtmp = qib0b*imu.dvbm;
	// vtmp1 = qni0' * [dvn+(wnin+wnie)Xvn-gn] * ts;
	tk += nts;
	CMat3 Ci0n = pos2Cen(CVect3(eth.pos.i,eth.wie*tk,0.0));
	CVect3 vtmp1 = Ci0n*(-eth.gn*nts);
	// Pib02 = Pib02 + vib0*ts, Pi02 = Pi02 + vi0*ts
	vib0 = vib0 + vtmp,		 vi0 = vi0 + vtmp1;
	Pib02 = Pib02 + vib0*nts, Pi02 = Pi02 + vi0*nts;
	//
	if(++t2>3*t0)
	{
		t0 = t1, Pib01 = tmpPib0, Pi01 = tmpPi0;
	}
	else if(t2>2*t0 && t1==t0)
	{
		t1 = t2, tmpPib0 = Pib02, tmpPi0 = Pi02;
	}
	//
	qib0b = qib0b*rv2q(imu.phim);
	// qnb=qni0*qiib0*qib0b
	CQuat qnb;
	if(t2<100)
	{
		qnb = CQuat(1.0);
	}
	else if(t2<1000)
	{
		qnb = CQuat(AlignCoarse(wmm, vmm, eth.pos.i));
	}
	else
	{
		qnb = ~CQuat(Ci0n)*CQuat(dv2att(Pi01, Pi02, Pib01, Pib02))*qib0b;
	}
	return qnb;
}

//#define assert(b)
BOOL assert(BOOL b)
{
	int res;

	if(b)
	{
		res = 1;
	}
	else
	{
		res = 0;
	}
	return res;
}

// determine the sign of 'val' with the sensitivity of 'eps'
int sign(double val, double eps)
{
	int s;

	if(val<-eps)
	{
		s = -1;
	}
	else if(val>eps)
	{
		s = 1;
	}
	else
	{
		s = 0; 
	}
	return s;
}

// set double value 'val' between range 'minVal' and 'maxVal'
double range(double val, double minVal, double maxVal)
{
	double res;

	if(val<minVal)
	{ 
		res = minVal; 
	}
	else if(val>maxVal)	
	{ 
		res = maxVal; 
	}
	else
	{ 
		res = val;
	}
	return res;
}

double atan2Ex(double y, double x)
{
	double res;

	if((sign(y)==0) && (sign(x)==0))
	{
		res = 0.0;
	}
	else
	{
		res = atan2(y, x);
	}
	return res;
}

/***************************  class CTDKF  *********************************/
CTDKF::CTDKF(int nq0, int nr0):CKalman(nq0,nr0)
{
	iter = -2;  ifn = 0;  nts = 0.0;
	Fk = Pk1 = CMat(nq,nq, 0.0);
	Pxz = Qk = Kk = tmeas = CVect(nr, 0.0);
	meanfn = O31;
}

int CTDKF::TDUpdate(CSINS &sins, double ts, int nStep)
{
	int measRes=0;

	if(nStep<=0) { nStep=2*(nq+nr)+3; }

	nts += ts; tk += ts;
	meanfn = meanfn+sins.fn; ifn++;
	for(int i=0; i<nStep; i++)
	{
		if(iter==-2)			// -2
		{
			CVect3 vtmp=meanfn*(1.0/ifn); meanfn = O31; ifn = 0;
			sins.fn=vtmp; SetFt(sins); sins.fn = vtmp;			
			MeasRearrange(sins);
		}
		else if(iter==-1)			// -1
		{
			Fk = ++(Ft*nts);  // Fk = I+Ft*ts
			Qk = Qt*nts;
			Xk = Fk*Xk;
			nts = 0.0;
		}
		else if(iter<nq)		// 0 -> (nq-1)
		{
			int row=iter;
			Pk1.SetRow(row, Fk.GetRow(row)*Pk);
		}
		else if(iter<2*nq)		// nq -> (2*nq-1)
		{
			int row=iter-nq;
			if(row==0) Fk = ~Fk;
			Pk.SetRow(row, Pk1.GetRow(row)*Fk);
			if(row==nq-1) {	Pk += Qk; }
		}
		else if(iter<2*nq+2*nr)	// (2*nq) -> (2*nq+2*nr-1)
		{
			int row=(iter-2*Ft.row)/2;
			int flag = measflag&(0x01<<row);
			if(flag)
			{
				measRes |= flag;
				if((iter-2*Ft.row)%2==0)
				{
					Hi = Hk.GetRow(row);
					Pxz = Pk*(~Hi);
					double Pzz = (Hi*Pxz)(0,0) + Rk(row);
					Kk = Pxz*(1.0/Pzz);
				}
				else
				{
					Xk += Kk*(Zk(row)-(Hi*Xk)(0,0));
					Pk -= Kk*(~Pxz);
				}
			}
		}
		else if(iter>=2*(nq+nr))	// 2*(nq+nr)
		{
			PkConstrain();
			symmetry(Pk);
			SetMeasFlag(0);
			iter = -3;
		}
		iter++;
	}
	Feedback(sins, 10.0, 1.0, 1.0, 100.0, 100.0);
	return measRes;
}