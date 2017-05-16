#include "rtklib.h"

static int g_reftype = 0;   /* 0: none ref, 1: static ref, 2: dynamic ref */

/* buffer for dynamic scenario */
static solbuf_t g_refpos = { 0 };
static solbuf_t g_refvel = { 0 };
static int g_lastposindex = 0;
static int g_lastvelindex = 0;

/* buffer for static scenario */
static double g_stacoor[3] = { 0 };

extern void setstaref(double coor[])
{
    g_stacoor[0] = coor[0];
    g_stacoor[1] = coor[1];
    g_stacoor[2] = coor[2];
    g_reftype = 1;
}

extern void setdynref(char *posfile, char *velfile)
{
    /* read position file */
    if (NULL != posfile)
    {
        memset(&g_refpos, 0, sizeof(solbuf_t));
        if (!readsol(&posfile, 1, &g_refpos))
        {
            g_reftype = 0;
            trace(6, "ERROR: read reference position file fail\n");
            return;
        }
    }
    if (NULL != velfile)
    {
        memset(&g_refvel, 0, sizeof(solbuf_t));
        if (!readsol(&velfile, 1, &g_refvel))
        {
            g_reftype = 0;
            trace(6, "ERROR: read reference velocity file fail\n");
            return;
        }
    }
    g_reftype = 2;
}

extern int  getrefpos(gtime_t t, double *refpos, double *refvel)
{
    static double rb[3] = { 0 }, pos[3] = { 0 }, venu[3] = {0};
	int bFindPos = 0, bFindVel = 0, i;
    
    if (g_reftype == 0) 
        return 0;

    if (g_reftype == 1)
    {
        refpos[0] = g_stacoor[0];
        refpos[1] = g_stacoor[1];
        refpos[2] = g_stacoor[2];
        refvel[0] = refvel[1] = refvel[2] = 0;
        return 1;
    }

    /* find position */
	for (i = g_lastposindex; i < g_refpos.n; i++)
    {
        if (fabs(timediff(t, g_refpos.data[i].time)) < 0.01)
        {
            refpos[0] = g_refpos.data[i].rr[0];
            refpos[1] = g_refpos.data[i].rr[1];
            refpos[2] = g_refpos.data[i].rr[2];
            rb[0] = g_refpos.data[i].rr[0];
            rb[1] = g_refpos.data[i].rr[1];
            rb[2] = g_refpos.data[i].rr[2];
            ecef2pos(rb, pos);
            g_lastposindex = i;
            bFindPos = 1;
        }
    }
    if (bFindPos == 0)
    {
        trace(6, "Cannot find reference position\n");
        g_lastposindex = 0;
        return 0;
    }

    /* find velocity */
    for (i = g_lastvelindex; i < g_refvel.n; i++)
    {
        if (fabs(timediff(t, g_refvel.data[i].time)) < 0.01)
        {
            venu[0] = g_refvel.data[i].rr[0];
            venu[1] = g_refvel.data[i].rr[1];
            venu[2] = g_refvel.data[i].rr[2];
            enu2ecef(pos, venu, refvel);
            g_lastvelindex = i;
            bFindVel = 1;
            break;
        }
    }
    if (!bFindVel)
    {
        trace(6, "Cannot find reference velocity\n");
        g_lastvelindex = 0;
        return 0;
    }
    return 1;
}