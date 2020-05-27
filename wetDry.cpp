#ifdef WIN32
#include <LibsDP.h>
#endif

#include <Arrays_dpreal.h> 
#include <IsOs.h>
#include <DpString.h>
#include <time.h>
#include <string>
#include <direct.h>
#include <math.h>
#include <array>
#include <supportFunct.cpp>

// wetting-drying according to Liang (2010)
inline void wetDry(dpreal& etaL, dpreal& etaR, 
				   dpreal& etaLWD, dpreal& etaRWD,			
				   dpreal& hL, dpreal& hR,
				   dpreal& hLWD, dpreal& hRWD,
				   dpreal& qxL, dpreal& qxR, 
				   dpreal& qxLWD, dpreal& qxRWD, 
				   dpreal& qyL, dpreal& qyR, 
				   dpreal& qyLWD, dpreal& qyRWD, 
				   dpreal& uL, dpreal& uR,
				   dpreal& vL, dpreal& vR,
				   dpreal& zbL, dpreal& zbR,
				   dpreal& zbface, dpreal& dz
				   )
{
	zbface = max(zbL, zbR);
	hLWD   = max(0.0, (etaL - zbface));
	hRWD   = max(0.0, (etaR - zbface));
	etaLWD = hLWD + zbface;
	etaRWD = hRWD + zbface;
	qxLWD  = uL*hLWD;
	qxRWD  = uR*hRWD;
	qyLWD  = vL*hLWD;
	qyRWD  = vR*hRWD;

	dz     = max(0.0, (zbface - etaL));
	zbface = zbface - dz;
	etaLWD = etaLWD - dz;
	etaRWD = etaRWD - dz;
}

