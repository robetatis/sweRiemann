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

// min-mod slope limiter
inline void minMod(Vec(dpreal)& inVect, dpreal& limslope)
{
	if      ((inVect(1)*inVect(2)) < 0.0)       limslope = 0.0;
	else if (fabs(inVect(1)) < fabs(inVect(2))) limslope = inVect(1);
	else                                        limslope = inVect(2);
}


// reconstruct 2nd order external face values
inline void recon2out(ArrayGen(dpreal)& arr, 
					  dpreal& arrRe, dpreal& arrLw, dpreal& arrUn, dpreal& arrDs, 
					  Vec(dpreal)& inVect, dpreal& limslope, 
					  dpreal& dx, dpreal& dy, int& i, int& j, int& nx, int& ny)
{
	// east
	if(i==(nx+1))
	{
		arrRe = arr(j,i+1);
	}
	else
	{
		inVect(1) = (arr(j,i+2) - arr(j,i+1))/dx;
		inVect(2) = (arr(j,i+1) - arr(j,i))/dx;
		minMod(inVect, limslope);
		arrRe = arr(j,i+1) - 0.5*dx*limslope;
	}

	// west
	if(i==2)
	{
		arrLw = arr(j,i-1);
	}
	else
	{
		inVect(1) = (arr(j,i-1) - arr(j,i-2))/dx;
		inVect(2) = (arr(j,i)   - arr(j,i-1))/dx;
		minMod(inVect, limslope);
		arrLw = arr(j,i-1) + 0.5*dx*limslope;
	}

	// north
	if(j==2)
	{
		arrUn = arr(j-1,i);
	}
	else
	{
		inVect(1) = (arr(j-2,i) - arr(j-1,i))/dx;
		inVect(2) = (arr(j-1,i) - arr(j,i))/dx;
		minMod(inVect, limslope);
		arrUn = arr(j-1,i) - 0.5*dy*limslope;
	}

	// south
	if(j==(ny+1))
	{
		arrDs = arr(j+1,i);
	}
	else
	{
		inVect(1) = (arr(j+1,i) - arr(j+2,i))/dx;
		inVect(2) = (arr(j,i)   - arr(j+1,i))/dx;
		minMod(inVect, limslope);
		arrDs = arr(j+1,i) + 0.5*dy*limslope;
	}
}


// reconstruct 2nd order inner face values
inline void recon2inn(ArrayGen(dpreal)& arr, 
					  dpreal& arrLe, dpreal& arrRw, dpreal& arrDn, dpreal& arrUs, 
					  Vec(dpreal)& inVect, dpreal& limslope, 
					  dpreal& dx, dpreal& dy, int& i, int& j)
{
	// x-direction
	inVect(1) = (arr(j,i)   - arr(j,i-1))/dx;
	inVect(2) = (arr(j,i+1) - arr(j,i))/dx;
	minMod(inVect, limslope);
	arrLe = arr(j,i) + 0.5*dx*limslope;
	arrRw = arr(j,i) - 0.5*dx*limslope;

	// y-direction
	inVect(1) = (arr(j,i)   - arr(j+1,i))/dy;
	inVect(2) = (arr(j-1,i) - arr(j,i))/dy;
	minMod(inVect, limslope);	
	arrDn = arr(j,i) + 0.5*dy*limslope; 
	arrUs = arr(j,i) - 0.5*dy*limslope;				

}


// reconstruct face values for dry cell
inline void recon1dry(ArrayGen(dpreal)& eta, ArrayGen(dpreal)& h, 
					  ArrayGen(dpreal)& qx, ArrayGen(dpreal)& qy, 
					  ArrayGen(dpreal)& u, ArrayGen(dpreal)& v, 
					  ArrayGen(dpreal)& zb,
					  dpreal& etaLe, dpreal& etaRe, dpreal& etaLw, dpreal& etaRw, 
					  dpreal& etaDn, dpreal& etaUn, dpreal& etaDs, dpreal& etaUs,
					  dpreal& hLe, dpreal& hRe, dpreal& hLw, dpreal& hRw, 
					  dpreal& hDn, dpreal& hUn, dpreal& hDs, dpreal& hUs,
					  dpreal& qxLe, dpreal& qxRe, dpreal& qxLw, dpreal& qxRw, 
					  dpreal& qxDn, dpreal& qxUn, dpreal& qxDs, dpreal& qxUs,
					  dpreal& qyLe, dpreal& qyRe, dpreal& qyLw, dpreal& qyRw, 
					  dpreal& qyDn, dpreal& qyUn, dpreal& qyDs, dpreal& qyUs,
					  dpreal& uLe, dpreal& uRe, dpreal& uLw, dpreal& uRw, 
					  dpreal& uDn, dpreal& uUn, dpreal& uDs, dpreal& uUs,
					  dpreal& vLe, dpreal& vRe, dpreal& vLw, dpreal& vRw, 
					  dpreal& vDn, dpreal& vUn, dpreal& vDs, dpreal& vUs,
					  dpreal& zbLe, dpreal& zbRe, dpreal& zbLw, dpreal& zbRw, 
					  dpreal& zbDn, dpreal& zbUn, dpreal& zbDs, dpreal& zbUs,
					  int& i, int& j
					  )
{
	etaLe = eta(j,i);          etaLw = eta(j,i-1);     etaDn = eta(j,i);         etaDs = eta(j+1,i);
	etaRe = eta(j,i+1);        etaRw = eta(j,i);       etaUn = eta(j-1,i);       etaUs = eta(j,i);

	hLe   = h(j,i);            hLw   = h(j,i-1);       hDn   = h(j,i);           hDs   = h(j+1,i);
	hRe   = h(j,i+1);          hRw   = h(j,i);         hUn   = h(j-1,i);         hUs   = h(j,i);

	qxLe  = 0.0;		       qxRw = 0.0;               qxDn = 0.0;               qxUs = 0.0;
	qxRe  = qx(j,i+1);         qxLw = qx(j,i-1);         qxUn = qx(j-1,i);         qxDs = qx(j+1,i);

	qyLe  = 0.0;		       qyRw  = 0.0;              qyDn = 0.0;               qyUs = 0.0;
	qyRe  = qy(j,i+1);         qyLw = qy(j,i-1);         qyUn = qy(j-1,i);         qyDs = qy(j+1,i);

	uLe   = 0.0;               uRw   = 0.0;              uDn   = 0.0;              uUs   = 0.0;
	uRe   = u(j,i+1);          uLw   = u(j,i-1);         uUn   = u(j-1,i);         uDs   = u(j+1,i);

	vLe   = 0.0;               vRw   = 0.0;              vDn   = 0.0;              vUs   = 0.0;
	vRe   = v(j,i+1);          vLw   = v(j,i-1);         vUn   = v(j-1,i);         vDs   = v(j+1,i);

	zbLe  = zb(j,i);           zbRw  = zb(j,i);          zbDn  = zb(j,i);          zbUs  = zb(j,i);
	zbRe  = zb(j,i+1);         zbLw  = zb(j,i-1);        zbUn  = zb(j-1,i);        zbDs  = zb(j+1,i);
}


// reconstruct face values for wet cell with at least 1 wet neighbor
inline void recon1dryWet(ArrayGen(dpreal)& eta, ArrayGen(dpreal)& h, 
						 ArrayGen(dpreal)& qx, ArrayGen(dpreal)& qy, 
						 ArrayGen(dpreal)& u, ArrayGen(dpreal)& v, 
						 ArrayGen(dpreal)& zb,
						 dpreal& etaLe, dpreal& etaRe, dpreal& etaLw, dpreal& etaRw, 
						 dpreal& etaDn, dpreal& etaUn, dpreal& etaDs, dpreal& etaUs,
						 dpreal& hLe, dpreal& hRe, dpreal& hLw, dpreal& hRw, 
						 dpreal& hDn, dpreal& hUn, dpreal& hDs, dpreal& hUs,
						 dpreal& qxLe, dpreal& qxRe, dpreal& qxLw, dpreal& qxRw, 
						 dpreal& qxDn, dpreal& qxUn, dpreal& qxDs, dpreal& qxUs,
						 dpreal& qyLe, dpreal& qyRe, dpreal& qyLw, dpreal& qyRw, 
						 dpreal& qyDn, dpreal& qyUn, dpreal& qyDs, dpreal& qyUs,
						 dpreal& uLe, dpreal& uRe, dpreal& uLw, dpreal& uRw, 
						 dpreal& uDn, dpreal& uUn, dpreal& uDs, dpreal& uUs,
						 dpreal& vLe, dpreal& vRe, dpreal& vLw, dpreal& vRw, 
						 dpreal& vDn, dpreal& vUn, dpreal& vDs, dpreal& vUs,
						 dpreal& zbLe, dpreal& zbRe, dpreal& zbLw, dpreal& zbRw, 
						 dpreal& zbDn, dpreal& zbUn, dpreal& zbDs, dpreal& zbUs,
						 int& i, int& j
						 )
{
	etaLe = eta(j,i);          etaRw = eta(j,i);         etaDn = eta(j,i);         etaUs = eta(j,i);
	etaRe = eta(j,i+1);        etaLw = eta(j,i-1);       etaUn = eta(j-1,i);       etaDs = eta(j+1,i);

	hLe   = h(j,i);            hRw   = h(j,i);           hDn   = h(j,i);           hUs   = h(j,i);
	hRe   = h(j,i+1);          hLw   = h(j,i-1);         hUn   = h(j-1,i);         hDs   = h(j+1,i);

	qxLe  = qx(j,i);           qxRw = qx(j,i);           qxDn = qx(j,i);           qxUs = qx(j,i);
	qxRe  = qx(j,i+1);         qxLw = qx(j,i-1);         qxUn = qx(j-1,i);         qxDs = qx(j+1,i);

	qyLe  = qy(j,i);           qyRw = qy(j,i);           qyDn = qy(j,i);           qyUs = qy(j,i);
	qyRe  = qy(j,i+1);         qyLw = qy(j,i-1);         qyUn = qy(j-1,i);         qyDs = qy(j+1,i);

	uLe   = u(j,i);            uRw   = u(j,i);           uDn   = u(j,i);           uUs   = u(j,i);
	uRe   = u(j,i+1);          uLw   = u(j,i-1);         uUn   = u(j-1,i);         uDs   = u(j+1,i);

	vLe   = v(j,i);            vRw   = v(j,i);           vDn   = v(j,i);           vUs   = v(j,i);
	vRe   = v(j,i+1);          vLw   = v(j,i-1);         vUn   = v(j-1,i);         vDs   = v(j+1,i);

	zbLe  = zb(j,i);           zbRw  = zb(j,i);          zbDn  = zb(j,i);          zbUs  = zb(j,i);
	zbRe  = zb(j,i+1);         zbLw  = zb(j,i-1);        zbUn  = zb(j-1,i);        zbDs  = zb(j+1,i);
}

