// spp.cpp : 定义控制台应用程序的入口点。
//

#include <iostream>
#include <string>
#include <fstream>
using namespace std;

#include "..\..\src\rtklib.h"
#pragma comment(lib,"..\\..\\lib\\rtklib.lib")


#define REF_ANALYSIS    1

#if REF_ANALYSIS
#define POSFILE     "D:\\TestData\\rtklib_ref.pos"
#else
#define POSFILE     "D:\\TestData\\NewVer.pos"
#endif

#define REFPOS      "D:\\StoredData\\i-MOTORS\\Train Test\\ref.pos"
#define REFVEL      "D:\\StoredData\\i-MOTORS\\Train Test\\ref_vel.pos"

#define CONFFILE    "D:\\TestData\\m8t_process.conf"
//#define CONFFILE    "D:\\StoredData\\rtklibexplorer\\union_1030\\demo5_m8n.conf"

#if   0
//#define ROVERNX     "D:\\StoredData\\i-MOTORS\\Field Test\\Group 4\\GS10_09.16o"
//#define ROVERNX     "D:\\StoredData\\i-MOTORS\\Field Test\\Group 4\\GS10_10.16o"
//#define ROVERNX     "D:\\StoredData\\i-MOTORS\\Field Test\\Group 4\\M8T_Tracker.obs"
#define ROVERNX     "D:\\StoredData\\i-MOTORS\\Field Test\\Group 4\\M8T_EVK.obs"
#define BASERNX     "D:\\StoredData\\i-MOTORS\\Field Test\\Group 4\\NGB2_All.16o"
#define NAVFILE     "D:\\StoredData\\i-MOTORS\\Field Test\\Group 4\\NGB2_All.16n"
#elif 1
#define ROVERNX     "D:\\StoredData\\i-MOTORS\\Train Test\\Point7_151135.obs"
#define BASERNX     "D:\\StoredData\\i-MOTORS\\Train Test\\base.obs"
#define NAVFILE     "D:\\StoredData\\i-MOTORS\\Train Test\\base.nav"
#elif 0
//#define ROVERNX     "D:\\StoredData\\rtklibsample\\oemv_20090515c.obs"
#define ROVERNX     "D:\\StoredData\\rtklibsample\\ubx_20090515c.obs"
#define BASERNX     "D:\\StoredData\\rtklibsample\\0263_20090515c.obs"
#define NAVFILE     "D:\\StoredData\\rtklibsample\\oemv_20090515c.nav"
#elif 0
#define ROVERNX     "D:\\StoredData\\comnav\\836DL\\802443191e.14O"
#define BASERNX     "D:\\StoredData\\comnav\\836DL\\802182191e.14O"
#define NAVFILE     "D:\\StoredData\\comnav\\836DL\\802182191e.14N"
#elif 0
/* rover and base are all ublox m8n */
#define ROVERNX     "D:\\StoredData\\rtklibexplorer\\argeles_car\\rover3.obs"
#define BASERNX     "D:\\StoredData\\rtklibexplorer\\argeles_car\\base3.obs"
#define NAVFILE     "D:\\StoredData\\rtklibexplorer\\argeles_car\\rover3.nav"
#elif 0
/* rover and base are all ublox m8n. <have problem with modified configuration at the beginning> */
#define ROVERNX     "D:\\StoredData\\rtklibexplorer\\argeles2_car\\rover1.obs"
#define BASERNX     "D:\\StoredData\\rtklibexplorer\\argeles2_car\\base1.obs"
#define NAVFILE     "D:\\StoredData\\rtklibexplorer\\argeles2_car\\rover1.nav"
#elif 0
#define ROVERNX     "D:\\StoredData\\rtklibexplorer\\collioure_car\\rover_1Hz.obs"
#define BASERNX     "D:\\StoredData\\rtklibexplorer\\collioure_car\\base_1Hz.obs"
#define NAVFILE     "D:\\StoredData\\rtklibexplorer\\collioure_car\\base_1Hz.nav"
#elif 0
#define ROVERNX     "D:\\StoredData\\rtklibexplorer\\niwot2_car\\rov2_m8t.obs"
#define BASERNX     "D:\\StoredData\\rtklibexplorer\\niwot2_car\\base2_m8t.obs"
#define NAVFILE     "D:\\StoredData\\rtklibexplorer\\niwot2_car\\rov2_m8t.nav"
#elif 0
#define ROVERNX     "D:\\StoredData\\rtklibexplorer\\radio_0906\\rov1.obs"
#define BASERNX     "D:\\StoredData\\rtklibexplorer\\radio_0906\\base1.obs"
#define NAVFILE     "D:\\StoredData\\rtklibexplorer\\radio_0906\\base1.nav"
#elif 0
#define ROVERNX     "D:\\StoredData\\rtklibexplorer\\roger_0904\\raw_rover.obs"
#define BASERNX     "D:\\StoredData\\rtklibexplorer\\roger_0904\\raw_base.obs"
#define NAVFILE     "D:\\StoredData\\rtklibexplorer\\roger_0904\\raw_base.nav"
#elif 0
#define ROVERNX     "D:\\StoredData\\rtklibexplorer\\union_0705\\rover3_m8t.obs"
#define BASERNX     "D:\\StoredData\\rtklibexplorer\\union_0705\\base3_m8t.obs"
#define NAVFILE     "D:\\StoredData\\rtklibexplorer\\union_0705\\base3_m8t.nav"
#elif 0
#define ROVERNX     "D:\\StoredData\\rtklibexplorer\\union_0705\\rover3.obs"
#define BASERNX     "D:\\StoredData\\rtklibexplorer\\union_0705\\base3_m8t.obs"
#define NAVFILE     "D:\\StoredData\\rtklibexplorer\\union_0705\\base3_m8t.nav"
#elif 0
//#define ROVERNX     "D:\\StoredData\\rtklibexplorer\\union_1030\\rover3.obs"
#define ROVERNX     "D:\\StoredData\\rtklibexplorer\\union_1030\\rover7.obs"
#define BASERNX     "D:\\StoredData\\rtklibexplorer\\union_1030\\tmgo3040.16o"
#define NAVFILE     "D:\\StoredData\\rtklibexplorer\\union_1030\\rover7.nav"
#endif

void main(int argc, char* argv[])
{
    printf("Rover: %s\n", ROVERNX);
    printf("Base: %s\n", BASERNX);
    printf("Solution: %s\n", POSFILE);

    char *rnxfile[] = { ROVERNX, BASERNX, NAVFILE };

    /* modify configuration */
    prcopt_t prcopt = prcopt_default;
    solopt_t solopt = solopt_default;
    filopt_t filopt = { 0 };
    gtime_t gt = { 0 };

    resetsysopts();
    if (loadopts(CONFFILE, sysopts))
    {
        getsysopts(&prcopt, &solopt, &filopt);
    }

    solopt.timef = TIMES_GPST;
    solopt.posf = SOLF_XYZ;
    solopt.trace = 6; // trace level 6
    solopt.sstat = 2; // output stat with residual
    
#if REF_ANALYSIS
    setdynref(REFPOS, REFVEL); /* use reference trajectory */
#endif

    postpos(gt, gt, 0, 0, &prcopt, &solopt, &filopt, rnxfile, 3, POSFILE, "rove", "base");

    //system("pause");
    return;
}

