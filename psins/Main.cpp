#include "kfapp.h"

int main(void)
{
	CFileRdWt fin("H:\\ygm2016\\精准\\PSINS-CPP\\data1.txt", 29);
	CFileRdWt fins("H:\\ygm2016\\精准\\PSINS-CPP\\ins.bin", 0);
	CFileRdWt fkf("H:\\ygm2016\\精准\\PSINS-CPP\\kf.bin", 0);

	CVect3 wm, vm, gpsVn, gpsPos, pos0(34.2120468*glv.deg,108.8323824*glv.deg,404.684), mag;
	int SatNum=0, res=0;
	double ts=0.01, PDOP, gpsYaw=0, wz;
		
	CRAvar ravar(5, ts);
	ravar.setR0(0.1, 0.01, 1*glv.deg, 0.1, 3.0);
	ravar.setTau(1.0, 1.0, 1.0, 1.0, 1.0);

	CCarAHRS kf(ts);
	for(int i=0; i<2500/ts; i++)
	{
		PDOP = 0.0;
		if(!fin.load(1)) // 读数据
            break;
		wm = CVect3(&fin.buff[1])*glv.deg*ts;  
        vm = CVect3(&fin.buff[4])*ts;
		gpsVn = CVect3(&fin.buff[10]); 
        PDOP = fin.buff[23];	
        SatNum = (int)fin.buff[24];
		gpsPos = CVect3(fin.buff[13]*glv.deg, fin.buff[14]*glv.deg, fin.buff[15]);
		gpsYaw = fin.buff[22]*glv.deg;	
        gpsYaw = gpsYaw>PI ? 2*PI-gpsYaw : -gpsYaw;
		mag = CVect3(&fin.buff[25]);
//		baro = fin.buff[28];
		
		ravar.Update(norm(kf.sins.an), norm(kf.sins.vn), norm(kf.sins.wnb), norm(gpsVn), normXY(gpsPos)*glv.Re);
		wz = ravar(2);
		if(PDOP>1.0 && PDOP<7.0 && SatNum>4 && wz<5*glv.dps)   // GPS速度组合
		{
//			kf.SetMeasGPSVn(gpsVn);
		}
		if(PDOP>1.0 && PDOP<7.0 && SatNum>4 && wz<10*glv.dps)   // GPS位置组合
		{
//			kf.SetMeasGPSPos(gpsPos);
		}
		if(ravar(0)<0.15 && wz<0.15*glv.dps)   // 静态检测
		{
			kf.SetMeasZUPT();
		}
		else if(wz<5*glv.dps && norm(kf.sins.vn)>10.0)  // 车载约束 
		{
			kf.SetMeasMC();
		}
		if(gpsYaw!=0 && wz<1*glv.dps)   // GPS双天线航向组合
		{
			kf.SetMeasGPSYaw(gpsYaw);
		}
		double nm = norm(mag);
		if(nm>200 && nm<800)
		{
			kf.SetMeasMag(mag);
		}
		res = kf.Update(wm, vm, ts);
		if(kf.yawAlignOK)
		{
			fkf<<kf;
			fins<<kf.sins<<ravar<<res<<gpsYaw;
		}
		if(i%10000==0)
			printf("%d\n", i/100);
	}
	return 0;
}


/* Use the following command to show the results in Matlab/PSINS Toolbox:
glvs
psinstypedef(153);
fid=fopen('ins.bin'); avp=fread(fid, [23,inf], 'double')'; fclose(fid); insplot(avp(:,1:16));
fid=fopen('kf.bin'); xkpk=fread(fid, [33,inf], 'double')'; fclose(fid); kfplot(xkpk(:,[1:15,17:31,end]));
*/
