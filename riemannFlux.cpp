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

// compute face fluxes in x-direction, HLLC approximate Riemann solver
inline void xFlux(Vec(dpreal)& Fface,																									
				  dpreal& uL, dpreal& uR, 
				  dpreal& vL, dpreal& vR, 
				  dpreal& qxL, dpreal& qxR,
				  dpreal& qyL, dpreal& qyR,
				  dpreal& hL, dpreal& hR, 
				  dpreal& etaL, dpreal& etaR,
				  dpreal& zbface,
				  Vec(dpreal)& UL, Vec(dpreal)& UR, 
				  Vec(dpreal)& FL, Vec(dpreal)& FR,
				  dpreal& ustar, dpreal& hstar, dpreal& SL, dpreal& SR, dpreal& Sstar,
				  Vec(dpreal)& FstarL, Vec(dpreal)& FstarR																			
				  )
{
	// compute face values of the vectors
	FL(1) = qxL;
	FL(2) = hL*pow(uL,2) + 0.5*9.81*(pow(etaL,2) - 2.0*etaL*zbface);
	FL(3) = uL*qyL;

	FR(1) = qxR;
	FR(2) = hR*pow(uR,2) + 0.5*9.81*(pow(etaR,2) - 2.0*etaR*zbface);
	FR(3) = uR*qyR;

	UL(1) = etaL;
	UL(2) = qxL;
	UL(3) = qyL;

	UR(1) = etaR;
	UR(2) = qxR;
	UR(3) = qyR;

	// compute wave velocities
	ustar = 0.5*(uL + uR) + sqrt(9.81*hL) - sqrt(9.81*hR);
	hstar = (1.0/9.81)*sqr(0.5*(sqrt(9.81*hL) + sqrt(9.81*hR)) + 0.25*(uL - uR));
	if(hL <= 1e-10) SL = uR - 2*sqrt(9.81*hR);
	else            SL = min(uL - sqrt(9.81*hL), ustar - sqrt(9.81*hstar));
	if(hR <= 1e-10) SR = uL + 2*sqrt(9.81*hL);
	else            SR = max(uR + sqrt(9.81*hR), ustar + sqrt(9.81*hstar));

	Sstar     = (SL*hR*(uR-SR) - SR*hL*(uL-SL))/(hR*(uR-SR) - hL*(uL-SL));
	FstarL(1) = (SR*FL(1) - SL*FR(1) + SL*SR*(UR(1) - UL(1)))/(SR - SL);
	FstarL(2) = (SR*FL(2) - SL*FR(2) + SL*SR*(UR(2) - UL(2)))/(SR - SL);

	if(SL<=0 && 0<=Sstar)     {FstarL(3) = vL*FstarL(1);}
	else if(Sstar<0 && 0<=SR) {FstarL(3) = vR*FstarL(1);}

	FstarR(1) = FstarL(1);
	FstarR(2) = FstarL(2);
	if(SL<=0 && 0<=Sstar)     {FstarR(3) = vL*FstarR(1);}
	else if(Sstar<0 && 0<=SR) {FstarR(3) = vR*FstarR(1);}

	// assign Fe according to Riemann solution structure
	if (0<=SL)
	{
		Fface(1) = FL(1);
		Fface(2) = FL(2);
		Fface(3) = FL(3);
	} 
	else if (SL<0 && 0<=Sstar)
	{
		Fface(1) = FstarL(1);
		Fface(2) = FstarL(2);
		Fface(3) = FstarL(3);
	}
	else if (Sstar<0 && 0<=SR)
	{
		Fface(1) = FstarR(1);
		Fface(2) = FstarR(2);
		Fface(3) = FstarR(3);
	}
	else
	{
		Fface(1) = FR(1);
		Fface(2) = FR(2);
		Fface(3) = FR(3);
	}	
}

// compute face fluxes in y-direction, HLLC approximate Riemann solver
inline void yFlux(Vec(dpreal)& Gface,
				  dpreal& uD, dpreal& uU,
				  dpreal& vD, dpreal& vU,
				  dpreal& qxD, dpreal& qxU,
				  dpreal& qyD, dpreal& qyU,
				  dpreal& hD, dpreal& hU,
				  dpreal& etaD, dpreal& etaU,
				  dpreal& zbface,
				  Vec(dpreal)& UD, Vec(dpreal)& UU, 
				  Vec(dpreal)& GD, Vec(dpreal)& GU,
				  dpreal& vstar, dpreal& hstar, dpreal& SD, dpreal& SU, dpreal& Sstar,
				  Vec(dpreal)& GstarD, Vec(dpreal)& GstarU
				  )
{	
	// compute face values of the vectors
	GD(1) = qyD;
	GD(2) = vD*qxD;
	GD(3) = hD*pow(vD,2) + 0.5*9.81*(pow(etaD,2) - 2.0*etaD*zbface);

	GU(1) = qyU;
	GU(2) = vU*qxU;
	GU(3) = hU*pow(vU,2) + 0.5*9.81*(pow(etaU,2) - 2.0*etaU*zbface);

	UD(1) = etaD;
	UD(2) = qxD;
	UD(3) = qyD;

	UU(1) = etaU;
	UU(2) = qxU;
	UU(3) = qyU;

	// compute wave velocites
	vstar = 0.5*(vD + vU) + sqrt(9.81*hD) - sqrt(9.81*hU);
	hstar = 1.0/9.81*sqr(0.5*(sqrt(9.81*hD) + sqrt(9.81*hU)) + 0.25*(vD - vU));

	if(hD <= 1e-10) SD = vU - 2*sqrt(9.81*hU);
	else            SD = min(vD - sqrt(9.81*hD), vstar - sqrt(9.81*hstar));
	if(hU <= 1e-10) SU = vD + 2*sqrt(9.81*hD);
	else            SU = max(vU + sqrt(9.81*hU), vstar + sqrt(9.81*hstar));

	Sstar     = (SD*hU*(vU-SU) - SU*hD*(vD-SD))/(hU*(vU-SU) - hD*(vD-SD));
	GstarD(1) = (SU*GD(1) - SD*GU(1) + SD*SU*(UU(1) - UD(1)))/(SU - SD);
	GstarD(3) = (SU*GD(3) - SD*GU(3) + SD*SU*(UU(3) - UD(3)))/(SU - SD);

	if(SD<=0 && 0<=Sstar)     {GstarD(2) = uD*GstarD(1);}
	else if(Sstar<0 && 0<=SU) {GstarD(2) = uU*GstarD(1);}

	GstarU(1) = GstarD(1);
	GstarU(3) = GstarD(3);
	if(SD<=0 && 0<=Sstar)     {GstarU(2) = uD*GstarU(1);}
	else if(Sstar<0 && 0<=SU) {GstarU(2) = uU*GstarU(1);}

	// assign Gs according to Riemann solution structure
	if (0<=SD)
	{
		Gface(1) = GD(1);
		Gface(2) = GD(2);
		Gface(3) = GD(3);
	} 
	else if (SD<0 && 0<=Sstar)
	{
		Gface(1) = GstarD(1);
		Gface(2) = GstarD(2);
		Gface(3) = GstarD(3);
	}
	else if (Sstar<0 && 0<=SU)
	{
		Gface(1) = GstarU(1);
		Gface(2) = GstarU(2);
		Gface(3) = GstarU(3);
	}
	else
	{
		Gface(1) = GU(1);
		Gface(2) = GU(2);
		Gface(3) = GU(3);
	}
}

