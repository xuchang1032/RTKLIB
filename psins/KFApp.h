/* KFApp c++ hearder file KFApp.h */
/*
	By     : Yan Gongmin @ NWPU
	Date   : 2017-04-29
	From   : College of Automation, 
	         Northwestern Polytechnical University, 
			 Xi'an 710072, China
*/

#ifndef _KFAPP_H
#define _KFAPP_H


#include "PSINS.h"

#define GPSVN	0
#define GPSPOS	3
#define ZUPT	6
#define CARMC	9		// Car Moving Constraint
#define GPSYAW	12

class PSINS_API CCarAHRS:public CTDKF
{
public:
	BOOL levelAlignOK, yawAlignOK, initPosOK;
	double measGPSYaw;
	CVect3 measGPSVn, measGPSPos, measINSvn, measCMvb, measMag;
	CAligni0 align;
	CSINS sins;

	CCarAHRS(double ts);
	void SetMeasGPSVn(CVect3& vnGPS);
	void SetMeasGPSPos(CVect3& posGPS);
	void SetMeasGPSYaw(double gpsYaw);
	void SetMeasZUPT(void);
	void SetMeasMC(void);
	void SetMeasMag(CVect3 &mag);
	virtual void MeasRearrange(CSINS &sins);
	int Update(CVect3 &wm, CVect3 &vm, double ts);
};

#endif

