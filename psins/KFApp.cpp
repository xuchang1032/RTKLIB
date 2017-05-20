#include "KFApp.h"

/***************************  class CCarAHRS  *********************************/
CCarAHRS::CCarAHRS(double ts):CTDKF(16,13)
{
	levelAlignOK = yawAlignOK = initPosOK = FALSE;
	iter = -2;  ifn = 0;
	measGPSVn = measGPSPos = measINSvn = measCMvb = measMag = O31;
	double sts = sqrt(ts);
	Pmax.Set2(10.0*glv.deg,10.0*glv.deg,30.0*glv.deg, 50.0,50.0,50.0, 1.0e4/glv.Re,1.0e4/glv.Re,1.0e4, 
		1000.0*glv.dph,1000.0*glv.dph,1000.0*glv.dph, 100.0*glv.mg,100.0*glv.mg,100.0*glv.mg, 10000.0*glv.ppm);
	Pmin.Set2(1.0*glv.min,1.0*glv.min,10.0*glv.min, 0.01,0.01,0.1, 1.0/glv.Re,1.0/glv.Re,0.1, 
		0.1*glv.dph,0.1*glv.dph,0.1*glv.dph, 0.1*glv.mg,0.1*glv.mg,0.1*glv.mg, 10*glv.ppm);
	Pk.SetDiag2(10.0*glv.deg,10.0*glv.deg,10.0*glv.deg, 1.0,1.0,1.0, 100.0/glv.Re,100.0/glv.Re,100.0, 
		100.0*glv.dph,100.0*glv.dph,100.0*glv.dph, 10.0*glv.mg,10.0*glv.mg,10.0*glv.mg, 1000*glv.ppm);
	Qt.Set2(1.1*glv.dpsh,1.1*glv.dpsh,1.1*glv.dpsh, 500.0*glv.ugpsHz,500.0*glv.ugpsHz,500.0*glv.ugpsHz, 0.0,0.0,0.0,
		0.0,0.0,0.0, 100.0*glv.ugpsh,100.0*glv.ugpsh,100.0*glv.ugpsh, 100.0*glv.ppmpsh);
	Rk.Set2(0.5, 0.5, 0.5, 10.0/glv.Re, 10.0/glv.Re, 10.0, 0.1/sts, 0.1/sts, 0.1/sts, 10.0/sts, 1000.0/sts, 1.0/sts, 1.0*glv.deg);
	SetHk(); 
	Hk(ZUPT,3) = Hk(ZUPT+1,4) = Hk(ZUPT+2,5) = 1.0;
	Hk(GPSYAW,2) = -1.0;  measGPSYaw = 0;
}

void CCarAHRS::SetMeasGPSVn(CVect3& vnGPS)
{
	measGPSVn = vnGPS;
	tmeas.dd[GPSVN] = tk;
}

void CCarAHRS::SetMeasGPSPos(CVect3& posGPS)
{
	measGPSPos = posGPS;
	tmeas.dd[GPSPOS] = tk;
	if(!initPosOK)
	{
		sins.pos = posGPS;  initPosOK = TRUE;
	}
}

void CCarAHRS::SetMeasZUPT(void)
{
	measINSvn = sins.vn;
	tmeas.dd[ZUPT] = tk;
}

void CCarAHRS::SetMeasMag(CVect3& mag)
{
	measMag = mag;
}

void CCarAHRS::SetMeasMC(void)
{
	Hk(CARMC+0,3) = sins.Cnb.e00; Hk(CARMC+0,4) = sins.Cnb.e10; Hk(CARMC+0,5) = sins.Cnb.e20;
	Hk(CARMC+1,3) = sins.Cnb.e01; Hk(CARMC+1,4) = sins.Cnb.e11; Hk(CARMC+1,5) = sins.Cnb.e21;
	Hk(CARMC+2,3) = sins.Cnb.e02; Hk(CARMC+2,4) = sins.Cnb.e12; Hk(CARMC+2,5) = sins.Cnb.e22;
	measCMvb = sins.vb;
	if(measCMvb.j>20) measCMvb.j-=20;
	else if(measCMvb.j<-1) measCMvb.j-=-1;
	tmeas.dd[CARMC] = tk;
}

void CCarAHRS::SetMeasGPSYaw(double gpsYaw)
{
	measGPSYaw = gpsYaw;
	tmeas.dd[GPSYAW] = tk;
}

void CCarAHRS::MeasRearrange(CSINS &sins)
{
	if(yawAlignOK && measGPSVn.k!=0 && tk-tmeas.dd[GPSVN]<0.5)   // GPSVn
	{
		*(CVect3*)&Zk.dd[GPSVN] = sins.vn-measGPSVn;
		measGPSVn = O31; SetMeasFlag(0x07);
	}
	if(yawAlignOK && measGPSPos.k!=0 && tk-tmeas.dd[GPSPOS]<0.5)   // GPSPos
	{
		*(CVect3*)&Zk.dd[GPSPOS] = sins.pos-measGPSPos-sins.eth.vn2dpos(sins.vn,tk-tmeas.dd[GPSPOS]-sins.nts);
		measGPSPos = O31; SetMeasFlag(0x38);
	}
	if(levelAlignOK && measINSvn.k!=0 && tk-tmeas.dd[ZUPT]<0.1)   // ZUPT
	{
		*(CVect3*)&Zk.dd[ZUPT] = measINSvn-sins.an*(tk-tmeas.dd[ZUPT]-sins.nts);
		measINSvn = O31; SetMeasFlag(0x01C0);
	}
	else if(levelAlignOK && measCMvb.k!=0 && tk-tmeas.dd[CARMC]<0.1 && tk-tmeas.dd[GPSVN]>10.0)	  // CARMC
	{
		*(CVect3*)&Zk.dd[CARMC] = measCMvb;
		measCMvb = O31; SetMeasFlag(0x0E00);
	}
	if(yawAlignOK && measGPSYaw!=0 && tk-tmeas.dd[GPSYAW]<0.5)   // GPSYaw
	{
		Zk.dd[GPSYAW] = sins.att.k - measGPSYaw + 2*glv.deg;
		measGPSYaw = 0;
//		if(fabs(Zk.dd[11])<5*glv.deg)
//			SetMeasFlag(0x1000);
	}
}

int CCarAHRS::Update(CVect3 &wm, CVect3 &vm, double ts)
{
	int res = 0;
	if(!levelAlignOK)
	{
		sins.qnb = align.Update(&wm, &vm, 1, ts);
		if(align.tk>5)
			levelAlignOK = TRUE;
	}
	if(levelAlignOK && !yawAlignOK)
	{
/*		double nm = norm(measMag);
		if(nm>520-100 && nm<520+250)  // 520
		{
			CVect3 Hn = sins.qnb*measMag;  // [XYZ, H, DEC, DIP, F] = wrldmagm(400, 34.21, 108.86, 2010)
			double nh = normXY(Hn);
			if(nh>320-100 && nh<320+200)
			{
				sins.qnb -= CVect3(0,0,atan2(-Hn.i, -Hn.j));
				yawAlignOK = TRUE;
			}
		}*/
		if(normXY(measGPSVn)>2.0)
		{
			CVect3 att(sins.qnb); att.k=atan2(-measGPSVn.i, measGPSVn.j);
			sins.qnb = CQuat(att);
			yawAlignOK = TRUE;
			sins.tk = tk;
		}
		if(measGPSYaw!=0)
		{
			CVect3 att(sins.qnb); att.k=measGPSYaw;
			sins.qnb = CQuat(att);
			yawAlignOK = TRUE;
			sins.tk = tk;
		}
	}

	sins.Update(&wm, &vm, 1, ts);
	res = TDUpdate(sins, ts, 20);

	return res;
}
