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

// get maximum of integer vector
int intVectMax(VecSimple(int)& inVect)
{
	int ll = inVect.size();
	int max = inVect(1);      
	for(int i=2; i<=ll; i++)
	{
		if(inVect(i) > max)
			max = inVect(i);
	}
	return max;
}


// get maximum of double vector
dpreal doubVectMax(Vec(dpreal)& inVect)
{
	int ll = inVect.size();
	dpreal max = inVect(1);      
	for(int i=2; i<=ll; i++)
	{
		if(inVect(i) >= max)
		{
			max = inVect(i);
		}
	}
	return max;
}


// get minimum of double vector
dpreal doubVectMin(Vec(dpreal)& inVect)
{
	int ll = inVect.size();
	dpreal min = inVect(1);      
	for(int i=2; i<=ll; i++)
	{
		if(inVect(i) <= min)
		{
			min = inVect(i);
		}
	}
	return min;
}


// get maximum in double array
dpreal doubArrMax(ArrayGen(dpreal)& inArr, int& nx, int& ny)
{
	dpreal max = inArr(1,1);

	for(int i=1; i<=nx; i++)
	{
		for(int j=1; j<=ny; j++)
		{
			if(inArr(j,i) >= max)
			{
				max = inArr(j,i);
			}
		}
	}
	return max;
}


// get minimum in double array
dpreal doubArrMin(ArrayGen(dpreal)& inArr, int& nx, int& ny)
{
	dpreal min = inArr(1,1);

	for(int i=1; i<=nx; i++)
	{
		for(int j=1; j<=ny; j++)
		{
			if(inArr(j,i) <= min)
			{
				min = inArr(j,i);
			}
		}
	}
	return min;
}


// find adress of mimimum element in vector
int maxEl(Vec(dpreal)& inVect)
{
	int ll = inVect.size();
	dpreal max = inVect(1);
	int imax=1;
	for(int i=2; i<=ll; i++)
	{
		if(inVect(i) >= max)
		{
			max = inVect(i);
			imax=i;
		}
	}
	return imax;
}


// interpolate hydrograph BC to current time step
inline void hydrgInterp(ArrayGen(dpreal)& hydrgArr, int& nrowHydrgArr, int& ncolHydrgArr, VecSimple(int)& hydrgNoSteps,
						Vec(dpreal)& currHydr, dpreal& slopehydrg, dpreal& t)
{
	if(ncolHydrgArr > 0)
	{
		// loop through each hydrograph in hydrgArr (k = number of iterations, not column number)
		for(int k=1; k<=(ncolHydrgArr/2); k++)   
		{

			// loop through all time steps in current hydrograph
			for(int g=1; g<hydrgNoSteps(k); g++) 
			{

				// do interpolation using slope that corresponds to current t
				if(t < hydrgArr(g+1, 2*k-1))     
				{
					slopehydrg  = (hydrgArr(g+1, 2*k) - hydrgArr(g, 2*k)) / (hydrgArr(g+1, 2*k-1) - hydrgArr(g, 2*k-1));
					currHydr(k) = hydrgArr(g, 2*k) + slopehydrg*(t - hydrgArr(g, 2*k-1));
					break;
				}
			}
		}
	}
}


// update state variables taking dry cells and negative water depths into account (Brufau et al., 2004)
void updPrim(ArrayGen(dpreal)& uc, ArrayGen(dpreal)& vc, ArrayGen(dpreal)& hc, ArrayGen(dpreal)& etac, ArrayGen(dpreal)& qxc, ArrayGen(dpreal)& qyc,
			 ArrayGen(dpreal)& zb, ArrayGen(dpreal)& hcarry,
			 int& i, int& j, Vec(dpreal)& Up, Vec(dpreal)& neighh, int& nhmax
			 )
{
	// integers used in selection of appropriate neighbor for extracting water in case of h<0
	int jp=0, ip=0;

	// preliminary update of state variables
	etac(j,i) = Up(1);
	qxc(j,i)  = Up(2);
	qyc(j,i)  = Up(3);
	hc(j,i)   = etac(j,i) - zb(j,i);

	// if cell dry, velocities and fluxes are zero
	if(hc(j,i) <= 1e-10) 
	{
		uc(j,i)=0.0; vc(j,i)=0.0;
		qxc(j,i)=0.0; qyc(j,i)=0.0;		

		// if h was negative, find neighbor with the most water and subtract the negative h from it to ensure mass conservation.
		// this is done only if after subtraction that neighbor has a depth > 1e-6 m. otherwise, negative water depth carried onto next
		// iteration until suitable neighbor is identified or the cell is sufficiently recharged
		if(hc(j,i) < 0.0)
		{
			// update hcarry with the negative water depth
			hcarry(j,i) = hcarry(j,i) + hc(j,i);

			// find the right neighbor
			neighh(1)=hc(j,i+1); neighh(2)=hc(j,i-1); neighh(3)=hc(j-1,i); neighh(4)=hc(j+1,i);
			nhmax=maxEl(neighh);

			if(nhmax == 1) {jp=0;  ip=1;}
			if(nhmax == 2) {jp=0;  ip=-1;}
			if(nhmax == 3) {jp=-1; ip=0;}
			if(nhmax == 4) {jp=1;  ip=0;}

			// if neighbor has enough water, take the negaive water depth from it
			if((hc(j+jp,i+ip) + hcarry(j,i)) > 1e-6)
			{
				// adjust qx and qy to keep u and v unmodified
				qxc(j+jp,i+ip) = qxc(j+jp,i+ip) * (hc(j+jp,i+ip) + hcarry(j,i))/hc(j+jp,i+ip);
				qyc(j+jp,i+ip) = qyc(j+jp,i+ip) * (hc(j+jp,i+ip) + hcarry(j,i))/hc(j+jp,i+ip);

				// remove water and reset hcarry
				hc(j+jp,i+ip) = hc(j+jp,i+ip) + hcarry(j,i);
				hcarry(j,i)   = 0.0;
			}
		}

		// "fill" negative depth
		hc(j,i)=0.0; etac(j,i)=zb(j,i);
	}
	else
	{
		uc(j,i) = qxc(j,i)/hc(j,i); 
		vc(j,i) = qyc(j,i)/hc(j,i);
	}

}


