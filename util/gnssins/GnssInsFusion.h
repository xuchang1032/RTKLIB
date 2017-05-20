#ifndef __GNSS_INS_FUSION__
#define __GNSS_INS_FUSION__

#include "rtklib.h"
#include "KFApp.h"

class CGnssInsFusion
{
public:
    CGnssInsFusion();
    ~CGnssInsFusion();

    CCarAHRS m_CarAHRS;
};

#endif