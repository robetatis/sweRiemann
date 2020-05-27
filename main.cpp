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
#include <inpoutp.cpp>
#include <setBC.cpp>
#include <setIC.cpp>
#include <supportFunct.cpp>
#include <timeLoop.cpp>

// 2D FV solver for inviscid shallow water equations on a Cartesian grid (Liang and Borthwick 2009, Liang 2010)

// to do:	
//          adaptive time step!!!
//			in readConFile: check file existence!!!
//							check file format dem, ic and man, make sure the three are the same form and that ncol and nrow match with actual array sizes
//			implement wq-beziehung in bcs
//			initialize with control file name from commandline/batch file
//			change control file reading --> loop till end of file and use switch to allocate control variables, such that ordering of keywords is not necessary
//			vassilis' sugestion -> do not compute fluxes twice!! Fe for cell j,i = Fw for cell j,i+1 ****** can reduce run time!!
// -------------------------------------------------------------------------------------------------------------------------------------------------------

// main program
int main(int argc, const char* argv[])
{
	initDiffpack(argc, argv);


	// control variables
	// --------------------------------------------------------------------------------------------------------------------------
	String dummy, demname, bcname, manname, etaicname,		// file names in control file
		uicname, vicname;
	String hout, etaout, vout, uout, vmagout, flowdout;		// switches for output files
	std::string rfpre, outdir; 
	dpreal tstop, dt;										// simulation time, time step
	dpreal outprf, outpc;									// output frequency (number of time steps) for results file and console, respectively
	dpreal Ctolup, Ctollow;									// max. and min. tolerable Courant values

	std::string cwd = getcwd(NULL,0);						// make directory for output files in current directory
	time_t sttime;											// variables for timing model run
	time(&sttime);		
	String sttimestr = ctime(&sttime);

	s_o << "\n---------------------------------------------------\nstarting time: " << sttimestr ;


	// declare and initialize variables to hold extent and resolution of computational domain
	// these data come from the DEM header
	// --------------------------------------------------------------------------------------------------------------------------
	int i=0;		 // column counter (x-coordinte, i grows eastwards, such that i=1 is the western edge of the domain)
	int j=0;         // row counter (y-coordinate, j grows southwards, such that j=1 is the northern edge of the domain)
	int nx=0, ny=0;  // no. columns and no. rows
	dpreal xll=0.0, yll=0.0, dx=0.0 ,dy=0.0, noData=0.0; // xll=coord. of lower left corner of domain, same for y, dx=dy=grid size


	// read control file
	// --------------------------------------------------------------------------------------------------------------------------
	readConFile(dummy, demname, bcname, etaicname, uicname, vicname, manname, 
		rfpre, outdir, tstop, dt, outprf, outpc, cwd,
		hout, etaout, vout, uout, vmagout, flowdout, Ctolup, Ctollow);


	// read DEM file header
	// --------------------------------------------------------------------------------------------------------------------------
	readDEMhead(demname, nx, ny, xll, yll, dx, dy, noData, dummy, i, j);


	// declare and initialize arrays for state variables
	// u         = x-velocity
	// v         = y-velocity
	// h         = water depth
	// qx        = x-unitary flow
	// qy        = y-unitary flow
	// eta       = water level
	// man       = manning coefficient
	// zb        = bottom elevation
	// vmag      = velocity magnitude = sqrt(vx^2 + vy^2)
	// flowd     = flow direction = atan2(vy/vx), atan2 considers sign ambiguity
	// "1,2,3"   = intermediate values for 4th-order Runge-Kutta
	// --------------------------------------------------------------------------------------------------------------------------
	ArrayGen(dpreal)
		u(ny+2,nx+2),  v(ny+2,nx+2),  h(ny+2,nx+2),  eta(ny+2,nx+2),  qx(ny+2,nx+2),  qy(ny+2,nx+2), zb(ny+2,nx+2), man(ny+2,nx+2), vmag(ny+2,nx+2), flowdir(ny+2,nx+2),
		u1(ny+2,nx+2), v1(ny+2,nx+2), h1(ny+2,nx+2), eta1(ny+2,nx+2), qx1(ny+2,nx+2), qy1(ny+2,nx+2),
		u2(ny+2,nx+2), v2(ny+2,nx+2), h2(ny+2,nx+2), eta2(ny+2,nx+2), qx2(ny+2,nx+2), qy2(ny+2,nx+2),
		u3(ny+2,nx+2), v3(ny+2,nx+2), h3(ny+2,nx+2), eta3(ny+2,nx+2), qx3(ny+2,nx+2), qy3(ny+2,nx+2);		

	u.fill(0.0); v.fill(0.0); h.fill(0.0); eta.fill(0.0); zb.fill(0.0); qx.fill(0.0); qy.fill(0.0); man.fill(0.0); vmag.fill(0.0); flowdir.fill(0.0);
	u1.fill(0.0); v1.fill(0.0); h1.fill(0.0); eta1.fill(0.0); qx1.fill(0.0); qy1.fill(0.0);
	u2.fill(0.0); v2.fill(0.0); h2.fill(0.0); eta2.fill(0.0); qx2.fill(0.0); qy2.fill(0.0);
	u3.fill(0.0); v3.fill(0.0); h3.fill(0.0); eta3.fill(0.0); qx3.fill(0.0); qy3.fill(0.0);	


	// declare and initialize variables to hold information on boundary conditions
	// --------------------------------------------------------------------------------------------------------------------------

	// flags for identifying boundary numbers along ghost cells
	Vec(dpreal) bcnE(ny+2), bcnW(ny+2), bcnS(nx+2), bcnN(nx+2);
	bcnE.fill(0.0); bcnW.fill(0.0); bcnS.fill(0.0); bcnN.fill(0.0);

	// flags for identifying boundary types along ghost cells
	Vec(dpreal) fube(ny+2), fubw(ny+2), fubs(nx+2), fubn(nx+2), 
		fope(ny+2), fopw(ny+2), fops(nx+2), fopn(nx+2), 
		fetae(ny+2), fetaw(ny+2), fetas(nx+2), fetan(nx+2), 
		fqe(ny+2), fqw(ny+2), fqs(nx+2), fqn(nx+2);

	fube.fill(0.0); fubw.fill(0.0); fubs.fill(0.0); fubn.fill(0.0);
	fope.fill(0.0); fopw.fill(0.0); fops.fill(0.0); fopn.fill(0.0);
	fetae.fill(0.0); fetaw.fill(0.0); fetas.fill(0.0); fetan.fill(0.0);
	fqe.fill(0.0); fqw.fill(0.0); fqs.fill(0.0); fqn.fill(0.0);

	// factors for flow distribution along boundaries
	Vec(dpreal) flowDistrE(ny+2), flowDistrW(ny+2), flowDistrS(nx+2), flowDistrN(nx+2); 

	flowDistrE.fill(0.0); flowDistrW.fill(0.0); flowDistrS.fill(0.0); flowDistrN.fill(0.0);

	// arrays for holding hydrographs (eta and Q)
	ArrayGen(dpreal) hydrgArr(0,0);
	int nrowHydrgArr=0;
	int ncolHydrgArr=0;

	// vector with no. of time steps in hydrographs
	VecSimple(int) hydrgNoSteps(1); hydrgNoSteps.fill(0);


	// read DEM file z-values
	// --------------------------------------------------------------------------------------------------------------------------
	readDEMz(demname, nx, ny, i, j, zb);


	// read manning coefficient
	// --------------------------------------------------------------------------------------------------------------------------
	readMan(manname, nx, ny, dummy, i, j, man);


	// set initial conditions
	// --------------------------------------------------------------------------------------------------------------------------
	setIC(h, eta, zb, u, v, qx, qy, nx, ny, dx, dy, i, j, etaicname, uicname, vicname);


	// set boundary conditions (default mode is closed boundary)
	// --------------------------------------------------------------------------------------------------------------------------
	setBC(dx, dy, xll, yll, nx, ny, tstop,
		u, v, h, eta, zb, man,
		i,j,
		bcname, dummy, 
		fube, fubw, fubs, fubn,
		fope, fopw, fops, fopn,
		fetae, fetaw, fetas, fetan,
		fqe, fqw, fqs, fqn,
		flowDistrE, flowDistrW, flowDistrS, flowDistrN,
		hydrgArr, nrowHydrgArr, ncolHydrgArr,
		hydrgNoSteps,
		bcnE, bcnW, bcnS, bcnN
		);


	// execute time loop
	// --------------------------------------------------------------------------------------------------------------------------
	timeLoop(u, v, h, eta, zb, man, qx, qy, vmag, flowdir,
		u1, v1, h1, eta1, qx1, qy1,
		u2, v2, h2, eta2, qx2, qy2,
		u3, v3, h3, eta3, qx3, qy3,
		dx, dy, i, j, nx, ny, dt, tstop, outprf, outpc,
		fube, fubw, fubs, fubn,
		fope, fopw, fops, fopn,
		fetae, fetaw, fetas, fetan,
		fqe, fqw, fqs, fqn,
		rfpre,
		xll, yll, noData, cwd, outdir,
		hydrgArr, nrowHydrgArr, ncolHydrgArr,
		hydrgNoSteps,
		bcnE, bcnW, bcnS, bcnN,
		flowDistrE, flowDistrW, flowDistrS, flowDistrN,
		hout, etaout, vout, uout, vmagout, flowdout, 
		Ctolup, Ctollow
		);


	// get and write run-time parameters to console
	// --------------------------------------------------------------------------------------------------------------------------
	runTime(sttime, sttimestr);

	return(0);			
}
