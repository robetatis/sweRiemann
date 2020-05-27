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

// set initial conditions on state variables: eta, u and v are read from ASCII raster, then h simply computed as eta - zb
void setIC(ArrayGen(dpreal)& h, ArrayGen(dpreal)& eta, ArrayGen(dpreal)& zb,
		   ArrayGen(dpreal)& u, ArrayGen(dpreal)& v, ArrayGen(dpreal)& qx, ArrayGen(dpreal)& qy, 
		   int& nx, int& ny, dpreal& dx, dpreal& dy, int& i, int& j, 
		   String& etaicname, String& uicname, String& vicname)

{
	s_o << "\nreading initial conditions\n-------------------------------------\n";

	// dummy variables. used to read the first 6 lines of the ASCII files
	String dummys;
	dpreal dummyd;
	dpreal unodata; // no data flag for u velocity
	dpreal vnodata; // no data flag for v velocity

	// create eta file stream
	ifstream etaicFile(etaicname.c_str());

	// read eta ASCII file header
	for (i=1; i<=6; i++)
		etaicFile >> dummys >> dummyd;

	// read initial eta data
	for(j=2; j<=(ny+1); j++)
		for(i=2; i<=(nx+1); i++)
			etaicFile >> eta(j,i);

	// compute h = eta - zb, handling eta < zb
	for(j=2; j<=(ny+1); j++)
	{
		for(i=2; i<=(nx+1); i++)
		{
			h(j,i) = eta(j,i) - zb(j,i);
			if(h(j,i)<=1e-10) 
			{
				h(j,i)   = 0.0;
				eta(j,i) = zb(j,i); 
			}
		}
	}

	// populate ghost cells with eta of adjacent internal cells
	i=1;
	for(j=2; j<=(ny+1); j++) eta(j,i) = eta(j,i+1);
	i=nx+2;
	for(j=2; j<=(ny+1); j++) eta(j,i) = eta(j,i-1);
	j=1;
	for(i=2; i<=(nx+1); i++) eta(j,i) = eta(j+1,i);
	j=ny+2;
	for(i=2; i<=(nx+1); i++) eta(j,i) = eta(j-1,i);


	// read initial velocity data; both components must exist, otherwise u and v are zero at simulation start
	if(uicname != "no" && vicname != "no") 
	{
		ifstream uicFile(uicname.c_str());
		ifstream vicFile(vicname.c_str());

		// read ASCII file header
		for (i=1; i<=5; i++)
		{
			uicFile >> dummys >> dummyd;
			vicFile >> dummys >> dummyd;
		}
		uicFile >> dummys >> unodata;
		vicFile >> dummys >> vnodata;

		// read u and v data
		for(j=2; j<=(ny+1); j++)
		{
			for(i=2; i<=(nx+1); i++)
			{
				uicFile >> u(j,i);
				vicFile >> v(j,i);

				// deal with no data flags by making u and v zero at those cells
				if(u(j,i) == unodata || v(j,i) == vnodata)
				{
					u(j,i) = 0; v(j,i) = 0;
				}

				// set initial qx and qy
				qx(j,i) = h(j,i)*u(j,i);
				qy(j,i) = h(j,i)*v(j,i);
			}
		}

		// populate ghost cells with u, v, qx and qy of adjacent internal cells
		i=1;
		for(j=2; j<=(ny+1); j++) {u(j,i) = u(j,i+1); v(j,i) = v(j,i+1); qx(j,i) = qx(j,i+1); qy(j,i) = qy(j,i+1);}
		i=nx+2;
		for(j=2; j<=(ny+1); j++) {u(j,i) = u(j,i-1); v(j,i) = v(j,i-1); qx(j,i) = qx(j,i-1); qy(j,i) = qy(j,i-1);}
		j=1;
		for(i=2; i<=(nx+1); i++) {u(j,i) = u(j+1,i); v(j,i) = v(j+1,i); qx(j,i) = qx(j+1,i); qy(j,i) = qy(j+1,i);}
		j=ny+2;
		for(i=2; i<=(nx+1); i++) {u(j,i) = u(j-1,i); v(j,i) = v(j-1,i); qx(j,i) = qx(j-1,i); qy(j,i) = qy(j-1,i);}
	}
}

