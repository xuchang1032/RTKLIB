/*------------------------------------------------------------------------------
* pntpos.c : dynamic standard point position
*
*          Copyright (C) 2007-2015 by chang.xu, All rights reserved.
*
* version : $Revision:$ $Date:$
*-----------------------------------------------------------------------------*/
#include "rtklib.h"

/* constants/macros ----------------------------------------------------------*/

#define SQR(x)      ((x)*(x))
#define SQRT(x)     ((x)<=0.0?0.0:sqrt(x))
#define MIN(x,y)    ((x)<=(y)?(x):(y))
#define ROUND(x)    (int)floor((x)+0.5)

#define VAR_POS		SQR(30.0)
#define VAR_VEL		SQR(10.0)
#define VAR_DT      SQR(30.0)
#define VAR_DDT		SQR(10.0)

#define NOISE_DT	SQR(5)
#define NOISE_DDT	SQR(0.5)

#define ERR_SAAS    0.3             /* saastamoinen model error std (m) */
#define ERR_BRDCI   0.5             /* broadcast iono model error factor */
#define ERR_CBIAS   0.3             /* code bias error std (m) */
#define REL_HUMI    0.7             /* relative humidity for saastamoinen model */


#define NP(opt)     ((opt)->dynamics?6:3) /* number of pos solution */
#define IC(s,opt)   (NP(opt)+(s))      /* state index of clocks (s=0:gps,1:glo,2:bds) */
#define ICD(opt)    (IC(0,opt)+NSYS)   /* state index of clock drift */
/* number of solutions */
#define NX_S(opt)     (ICD(opt)+1) /* number of estimated states */

int getnx_single(const prcopt_t *opt)
{
	return NX_S(opt);
}

/* pseudorange measurement error variance ------------------------------------*/
static double varerr(const prcopt_t *opt, double el, int sys)
{
	double fact,varr;
	fact=sys==SYS_GLO?EFACT_GLO:(sys==SYS_SBS?EFACT_SBS:EFACT_GPS);
	varr=SQR(opt->err[0])*(SQR(opt->err[1])+SQR(opt->err[2])/sin(el));
	if (opt->ionoopt==IONOOPT_IFLC) varr*=SQR(3.0); /* iono-free */
	return SQR(fact)*varr;
}

/* initialize state and covariance -------------------------------------------*/
static void initx(rtk_t *rtk, double xi, double var, int i)
{
	int j;
	rtk->x[i]=xi;
	for (j=0;j<rtk->nx;j++) {
		rtk->P[i+j*rtk->nx]=rtk->P[j+i*rtk->nx]=i==j?var:0.0;
	}
}

/* pseudorange and doppler residual */
static int res_dynspp(rtk_t *rtk, const obsd_t *obs, int n, const double *rs,
					  const double *dts, const double *vare, const int *svh,
					  const nav_t *nav, double *x, double *v,
					  double *H, double *R, double *azel)
{
	prcopt_t *opt=&(rtk->opt);
	double r,dion,dtrp,vion,vtrp,rr[3],pos[3],dtr,e[3],lam_L1;
	int i,j,nv=0,sys,nx=rtk->nx;
	double lam,rate,E[9],a[3],vs[3],cosel;
	double var[MAXOBS*2];
	double clockjump=0;

	for (i=0;i<3;i++) rr[i]=x[i];

	ecef2pos(rr,pos); xyz2enu(pos,E);

	for (i=0;i<n&&i<MAXOBS;i++) {
		azel[i*2]=azel[1+i*2]=0.0;
		if (!(sys=satsys(obs[i].sat,NULL))) continue;

		/* reject duplicated observation data */
		if (i<n-1&&i<MAXOBS-1&&obs[i].sat==obs[i+1].sat) {
			trace(2,"duplicated observation data %s sat=%2d\n",
				time_str(obs[i].time,3),obs[i].sat);
			i++;
			continue;
		}
		/* geometric distance/azimuth/elevation angle */
		if ((r=geodist(rs+i*6,rr,e))<=0.0||
			satazel(pos,e,azel+i*2)<opt->elmin) continue;
		
		/* excluded satellite? */
		if (satexclude(obs[i].sat,svh[i],opt)) continue;

		/* ionospheric corrections */
		if (!ionocorr(obs[i].time,nav,obs[i].sat,pos,azel+i*2,
			IONOOPT_BRDC,&dion,&vion)) continue;

		/* GPS-L1 -> L1/B1 */
		if ((lam_L1=nav->lam[obs[i].sat-1][0])>0.0) {
			dion*=SQR(lam_L1/lam_carr[0]);
		}

		/* tropospheric corrections */
		if (!tropcorr(obs[i].time,nav,pos,azel+i*2,
			TROPOPT_SAAS,&dtrp,&vtrp)) {
				continue;
		}

		if (sys == SYS_GPS) dtr=x[IC(0,opt)];
		else if (sys == SYS_GLO) dtr=x[IC(1,opt)];
		else if (sys == SYS_CMP) dtr=x[IC(2,opt)];

		/* psudorange residual */		
		v[nv]=obs[i].P[0]-(r+dtr-CLIGHT*dts[i*2]+dion+dtrp);
		
		trace(6,"dynspp psrres sat %02d v %7.3f P %15.3f r %15.3f dtr %12.3f dts %10.3f dion %5.2f dtrp %5.2f\n",
			obs[i].sat,v[nv],obs[i].P[0],r,dtr,CLIGHT*dts[i*2],dion,dtrp);

		if (fabs(v[nv]) > 1000) continue;

		/* design matrix */
		for (j=0;j<rtk->nx;j++) 
		{
			if (j<3) H[j+nv*rtk->nx] = -e[j]; //x y z
			else if (j<NP(opt)) H[j+nv*rtk->nx] = 0; //vx vy vz
			else if (j<ICD(opt)) //dtg dtr dtc
			{
				H[IC(0,opt)+nv*rtk->nx] = (sys==SYS_GPS?1:0);
				H[IC(1,opt)+nv*rtk->nx] = (sys==SYS_GLO?1:0);
				H[IC(2,opt)+nv*rtk->nx] = (sys==SYS_CMP?1:0);
			}
			else H[IC(2,opt)+nv*rtk->nx] = 0; //ddt
		}
		
		/* psudorange variance */
		var[nv]=varerr(opt,azel[1+i*2],sys)+vare[i]+vion+vtrp;
		var[nv]*=50;
		nv++;
		
		lam=nav->lam[obs[i].sat-1][0];
		if (obs[i].D[0]==0.0||lam==0.0||norm(rs+3+i*6,3)<=0.0) {
			continue;
		}
		/* line-of-sight vector in ecef */
		cosel=cos(azel[1+i*2]);
		a[0]=sin(azel[i*2])*cosel;
		a[1]=cos(azel[i*2])*cosel;
		a[2]=sin(azel[1+i*2]);
		matmul("TN",3,1,3,1.0,E,a,0.0,e);

		/* satellite velocity relative to receiver in ecef */
		for (j=0;j<3;j++) vs[j]=rs[j+3+i*6]-x[j+3];

		/* range rate with earth rotation correction */
		rate=dot(vs,e,3)+OMGE/CLIGHT*(rs[4+i*6]*rr[0]+rs[1+i*6]*x[3]-
			rs[3+i*6]*rr[1]-rs[  i*6]*x[4]);

		/* doppler residual */
		v[nv]=-lam*obs[i].D[0]-(rate+x[ICD(opt)]-CLIGHT*dts[1+i*2]);

		trace(6, "dynspp dopres sat %2d v %6.3f D %7.3f rate %7.3f clkdrift %6.3f satdrift %6.3f\n",
			obs[i].sat,v[nv],-lam*obs[i].D[0],rate,x[ICD(opt)],CLIGHT*dts[1+i*2]);

		if (fabs(v[nv]) > 100) continue;

		/* design matrix */
		for (j=0;j<rtk->nx;j++) 
		{
			if (j<3) H[j+nv*rtk->nx] = 0; //x y z
			else if (j<NP(opt)) H[j+nv*rtk->nx] = -e[j-3]; //vx vy vz
			else if (j<ICD(opt)) H[j+nv*rtk->nx] = 0; //dtg dtr dtc
			else H[j+nv*rtk->nx] = 1; //ddt
		}

		/* doppler variance */
		var[nv++]=(varerr(opt,azel[1+i*2],sys)+vare[i]+vion+vtrp)/200.;
	}

	for (i=0;i<nv;i++) 
	{
		for (j=0;j<nv;j++)
			R[i+j*nv]=(i==j?var[i]:0.0);
	}
	
	trace(6,"x=\n"); tracemat(6,x, 1,nx,8,3);
	trace(6,"v=\n"); tracemat(6,v, 1,nv,8,3);
	trace(6,"H=\n"); tracemat(6,H,nx,nv,8,3);
	trace(6,"R=\n"); tracemat(6,R,nv,nv,8,5);
	return nv;
}


/*
    X = [x y z vx vy vz dtg dtr dtc ddt]

		[1 0 0 t 0 0 0 0 0 0]
		[0 1 0 0 t 0 0 0 0 0]
		[0 0 1 0 0 t 0 0 0 0]
		[0 0 0 1 0 0 0 0 0 0]
		[0 0 0 0 1 0 0 0 0 0]
	F =	[0 0 0 0 0 1 0 0 0 0]
		[0 0 0 0 0 0 1 0 0 t]
		[0 0 0 0 0 0 0 1 0 t]
		[0 0 0 0 0 0 0 0 1 t]
		[0 0 0 0 0 0 0 0 0 1]
*/
static int udstate_spp(rtk_t *rtk)
{
	prcopt_t *opt = &(rtk->opt);
	double tt = fabs(rtk->tt);
	double *F,*FP,*xp,pos[3],*Q, *Qv;
	int i,j;
	double Qx=10., Qy=10., Qz=10., Qf=10, Qg=1.;
	double coef1, coef2, coef3;

	/* initialize state for first epoch */
	if (norm(rtk->x,3)<=0.0) {
		for (i=0;i<3;i++) initx(rtk,rtk->sol.rr[i],VAR_POS,i); //x y z
		for (i=3;i<NP(opt);i++) initx(rtk,rtk->sol.rr[i],VAR_VEL,i); //vx vy vz
		initx(rtk,rtk->sol.dtr[0],VAR_DT,IC(0,opt)); //dt gps
		initx(rtk,(rtk->sol.dtr[0]+rtk->sol.dtr[1]),VAR_DT,IC(1,opt)); //dt glo
		initx(rtk,(rtk->sol.dtr[0]+rtk->sol.dtr[3]),VAR_DT,IC(2,opt)); //dt bds
		initx(rtk,rtk->sol.ddt,VAR_DDT,ICD(opt)); //ddt

		trace(6,"x(0)=\n"); tracemat(6,rtk->x,1,rtk->nx,13,4);

		return 0;
	}

	/* state transition of position/velocity */
	F=eye(rtk->nx); FP=mat(rtk->nx,rtk->nx); xp=mat(rtk->nx,1);
	F[0+(3)*rtk->nx]=tt;
	F[1+(4)*rtk->nx]=tt;
	F[2+(5)*rtk->nx]=tt;
	F[6+(9)*rtk->nx]=tt;
	F[7+(9)*rtk->nx]=tt;
	F[8+(9)*rtk->nx]=tt;

	trace(6,"F=\n"); tracemat(6,F,rtk->nx,rtk->nx,8,1);
	trace(6,"X(t-1)=\n"); tracemat(6,rtk->x,1,rtk->nx,8,1);
	trace(6,"P(t-1)=\n"); tracemat(6,rtk->P,rtk->nx,rtk->nx,8,1);
	
	/* x=F*x, P=F*P*F+Q */
	matmul("NN",rtk->nx,1,rtk->nx,1.0,F,rtk->x,0.0,xp);
	matcpy(rtk->x,xp,rtk->nx,1);
	matmul("NN",rtk->nx,rtk->nx,rtk->nx,1.0,F,rtk->P,0.0,FP);
	matmul("NT",rtk->nx,rtk->nx,rtk->nx,1.0,FP,F,0.0,rtk->P);

	/* process noise */
	Q = mat(rtk->nx,rtk->nx); Qv = mat(rtk->nx,rtk->nx);
	coef1 = rtk->tt; coef2 = coef1*coef1/2; coef3 = coef1*coef1*coef1/3;
	Q[0+0*rtk->nx] = coef3*Qx;
	Q[1+1*rtk->nx] = coef3*Qy;
	Q[2+2*rtk->nx] = coef3*Qz;
	Q[3+3*rtk->nx] = coef1*Qx;
	Q[4+4*rtk->nx] = coef1*Qy;
	Q[5+5*rtk->nx] = coef1*Qz;
	Q[0+3*rtk->nx] = Q[3+0*rtk->nx] = coef2*Qx;
	Q[1+4*rtk->nx] = Q[4+1*rtk->nx] = coef2*Qy;
	Q[2+5*rtk->nx] = Q[5+2*rtk->nx] = coef2*Qz;
	Q[IC(0,opt)+IC(0,opt)*rtk->nx] = coef1*Qf + coef3*Qg;
	Q[IC(1,opt)+IC(1,opt)*rtk->nx] = coef1*Qf + coef3*Qg;
	Q[IC(2,opt)+IC(2,opt)*rtk->nx] = coef1*Qf + coef3*Qg;
	Q[ICD(opt)+ICD(opt)*rtk->nx] = coef1*Qg;
	Q[IC(0,opt)+ICD(opt)*rtk->nx] = Q[ICD(opt)+IC(0,opt)*rtk->nx] = coef2*Qg;
	Q[IC(1,opt)+ICD(opt)*rtk->nx] = Q[ICD(opt)+IC(1,opt)*rtk->nx] = coef2*Qg;
	Q[IC(2,opt)+ICD(opt)*rtk->nx] = Q[ICD(opt)+IC(2,opt)*rtk->nx] = coef2*Qg;
	ecef2pos(rtk->x,pos);
	//covecef(pos,Q,Qv);
	for (i=0;i<rtk->nx;i++) 
	{
		for (j=0;j<rtk->nx;j++)
			rtk->P[i+j*rtk->nx]+=Q[i+j*rtk->nx];
	}

	trace(6,"X'(t)=\n"); tracemat(6,rtk->x,1,rtk->nx,8,1);
	trace(6,"Qxyz=\n"); tracemat(6,Q,rtk->nx,rtk->nx,8,1);
	trace(6,"P'(t)=\n"); tracemat(6,rtk->P,rtk->nx,rtk->nx,8,1);
	free(F); free(FP); free(xp);
	return 1;
}

extern void dynspp(rtk_t *rtk, const obsd_t *obs, int n, const nav_t *nav)
{
	prcopt_t *opt = &(rtk->opt);
	double *rs,*dts,*vare,*xp,*Pp,*v,*H,*R,*azel_;
	int i,info,nv,svh[MAXOBS];
	rs=mat(6,n); dts=mat(2,n); vare=mat(1,n);
	
	xp=mat(rtk->nx,1); Pp=zeros(rtk->nx,rtk->nx);
	nv=n*2; v=mat(nv,1); H=mat(rtk->nx,nv); R=mat(nv,nv); azel_=zeros(2,n);

	/* satellite positons, velocities */
	satposs(obs[0].time,obs,n,nav,opt->sateph,rs,dts,vare,svh);

	/* adjust time ms jump */
	if (fabs(fabs(rtk->sol.dtr[0]-rtk->x[IC(0,opt)]) - 0.001*CLIGHT) < 1000)
	{
		rtk->x[IC(0,opt)] = rtk->sol.dtr[0];
		rtk->x[IC(1,opt)] = rtk->sol.dtr[0]+rtk->sol.dtr[1];
		rtk->x[IC(2,opt)] = rtk->sol.dtr[0]+rtk->sol.dtr[2];
	}

	/* state propogation */
	if (0==udstate_spp(rtk)) return;

	matcpy(xp,rtk->x,rtk->nx,1);
	/* code and doppler residual */
	if ((nv=res_dynspp(rtk,obs,n,rs,dts,vare,svh,nav,xp,v,H,R,azel_))<=0)
		return;

	/* measurement update */
	matcpy(Pp,rtk->P,rtk->nx,rtk->nx);
	if ((info=filter(xp,Pp,H,v,R,rtk->nx,nv))) {
		trace(2,"spp filter error %s info=%d\n",time_str(rtk->sol.time,0),
			info);
		return;
	}

	/* update state and covariance matrix */
	matcpy(rtk->x,xp,rtk->nx,1);
	matcpy(rtk->P,Pp,rtk->nx,rtk->nx);

	trace(6,"X(t)=\n"); tracemat(6,rtk->x,1,rtk->nx,8,1);
	trace(6,"P(t)=\n"); tracemat(6,rtk->P,rtk->nx,rtk->nx,8,1);

	/* update solution status */
	rtk->sol.ns=0;
	for (i=0;i<n&&i<MAXOBS;i++) {
		if (!rtk->ssat[obs[i].sat-1].vsat[0]) continue;
		rtk->sol.ns++;
	}
	rtk->sol.stat=SOLQ_DGPS;
	for (i=0;i<3;i++) {
		rtk->sol.rr[i]=rtk->x[i];
		rtk->sol.rr[i+3]=rtk->x[i+3];
		rtk->sol.qr[i]=(float)rtk->P[i+i*rtk->nx];
	}
	rtk->sol.qr[3]=(float)rtk->P[1];
	rtk->sol.qr[4]=(float)rtk->P[2+rtk->nx];
	rtk->sol.qr[5]=(float)rtk->P[2];

	free(rs); free(dts); free(vare);
	free(xp); free(Pp); free(v); free(H); free(R);

	trace(6, "END OF PROCESS\n\n");
}