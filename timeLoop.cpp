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
#include <face.cpp> 
#include <inpoutp.cpp>
#include <main.cpp>
#include <riemannFlux.cpp>
#include <supportFunct.cpp>
#include <wetDry.cpp>

// time loop
void timeLoop(ArrayGen(dpreal)& u, ArrayGen(dpreal)& v, ArrayGen(dpreal)& h, ArrayGen(dpreal)& eta, ArrayGen(dpreal)& zb, ArrayGen(dpreal)& man,
			  ArrayGen(dpreal)& qx, ArrayGen(dpreal)& qy, ArrayGen(dpreal)& vmag, ArrayGen(dpreal)& flowdir,
			  ArrayGen(dpreal)& u1, ArrayGen(dpreal)& v1, ArrayGen(dpreal)& h1, ArrayGen(dpreal)& eta1, ArrayGen(dpreal)& qx1, ArrayGen(dpreal)& qy1,
			  ArrayGen(dpreal)& u2, ArrayGen(dpreal)& v2, ArrayGen(dpreal)& h2, ArrayGen(dpreal)& eta2, ArrayGen(dpreal)& qx2, ArrayGen(dpreal)& qy2,
			  ArrayGen(dpreal)& u3, ArrayGen(dpreal)& v3, ArrayGen(dpreal)& h3, ArrayGen(dpreal)& eta3, ArrayGen(dpreal)& qx3, ArrayGen(dpreal)& qy3,
			  dpreal& dx, dpreal& dy, int& i, int& j, int& nx, int& ny, dpreal& dt, dpreal& tstop, dpreal& outprf, dpreal& outpc,
			  Vec(dpreal)& fube, Vec(dpreal)& fubw, Vec(dpreal)& fubs, Vec(dpreal)& fubn,
			  Vec(dpreal)& fope, Vec(dpreal)& fopw, Vec(dpreal)& fops, Vec(dpreal)& fopn,
			  Vec(dpreal)& fetae, Vec(dpreal)& fetaw, Vec(dpreal)& fetas, Vec(dpreal)& fetan,
			  Vec(dpreal)& fqe, Vec(dpreal)& fqw, Vec(dpreal)& fqs, Vec(dpreal)& fqn,
			  std::string& rfpre,
			  dpreal& xll, dpreal& yll, dpreal& noData,
			  std::string& cwd, std::string& outdir,
			  ArrayGen(dpreal)& hydrgArr, int& nrowHydrgArr, int& ncolHydrgArr, VecSimple(int)& hydrgNoSteps,
			  Vec(dpreal)& bcnE, Vec(dpreal)& bcnW, Vec(dpreal)& bcnS, Vec(dpreal)& bcnN,
			  Vec(dpreal)& flowDistrE, Vec(dpreal)& flowDistrW, Vec(dpreal)& flowDistrS, Vec(dpreal)& flowDistrN,
			  String& hout,String& etaout,String& vout,String& uout, String& vmagout, String& flowdout,
			  dpreal& Ctolup, dpreal& Ctollow
			  )
{
	s_o << "\nstarting solver\n-------------------------------------\n";


	dpreal t = 0.0;		 // time counter
	dpreal toutc = 0.0;  // counter for console output
	dpreal toutf = 0.0;  // counter for file output


	// variables for computing courant number
	ArrayGen(dpreal) Cx(ny, nx);								// Courant number in x-direction
	ArrayGen(dpreal) Cy(ny, nx);								// Courant number in y-direction
	Cx.fill(0.0); Cy.fill(0.0);

	dpreal Cxmax=0.0;											// max courant number in domain, x-direction
	dpreal Cymax=0.0;											// max courant number in domain, y-direction
	dpreal Cmax=0.0;


	// variables for building output file names	
	int wwfile=0, wwcons=0, precfile=0, preccons=0;
	std::string etaname, hname, uname, vname, vmagname,
		flowdname, tString, fileext;


	// variables used in temporal hydrograph interpolation
	Vec(dpreal) currHydr(ncolHydrgArr/2); currHydr.fill(0.0);	 // vector for holding current hydrograph values
	int indCurrHydr=0;											 // index for selecting current hydrograph from currHydr	
	dpreal slopehydrg=0.0;                                       // used in hydrograph interpolation


	// variables for handling source terms
	dpreal Sfx=0.0, Sfy=0.0;
	dpreal Sfxp=0.0, Sfyp=0.0;
	dpreal Dx=0.0, Dy=0.0;
	dpreal Fx=0.0, Fy=0.0;
	dpreal Cf=0.0;												 
	dpreal dzw=0.0, dze=0.0, dzs=0.0, dzn=0.0;
	dpreal Sbw=0.0, Sbe=0.0, Sbn=0.0, Sbs=0.0;


	// variables and arrays to hold mass balance correction for individual cells for the case of h<0
	ArrayGen(dpreal) hcarry(ny+2, nx+2); hcarry.fill(0.0);
	int nhmax=0;


	// vectors to hold conserved quantities and Runge-Kutta coefficients (K1 etc.)
	Vec(dpreal) U(3), Up(3), Fe(3), Fw(3), Gs(3), Gn(3), S(3);
	ArrayGen(dpreal) K11(ny+2, nx+2), K12(ny+2, nx+2), K13(ny+2, nx+2);
	ArrayGen(dpreal) K21(ny+2, nx+2), K22(ny+2, nx+2), K23(ny+2, nx+2);
	ArrayGen(dpreal) K31(ny+2, nx+2), K32(ny+2, nx+2), K33(ny+2, nx+2);
	ArrayGen(dpreal) K41(ny+2, nx+2), K42(ny+2, nx+2), K43(ny+2, nx+2);
	K11.fill(0.0); K12.fill(0.0); K13.fill(0.0);
	K21.fill(0.0); K22.fill(0.0); K23.fill(0.0);
	K31.fill(0.0); K32.fill(0.0); K33.fill(0.0);
	K41.fill(0.0); K42.fill(0.0); K43.fill(0.0);


	// variables for computation of slope limiters and face values, inVect(1) = backward slope, inVect(2) = forward slope 	
	Vec(dpreal)         inVect(2); inVect.fill(0.0);            		 
	dpreal              limslope=0.0;


	// variables for Riemann solver: 
	// u, v, qx, qy, eta, h, zb = state variables
	// R, L, U, D = right, left, upper and lower sides of cell interface
	// e, w, s, n = eastern, western, northern and southern cell faces
	Vec(dpreal) neighh(4); neighh.fill(0.0); // vector used to store neighboring cell h values

	dpreal qxLe=0.0,  qxRe=0.0;
	dpreal qxLw=0.0,  qxRw=0.0;
	dpreal qxUs=0.0,  qxDs=0.0;
	dpreal qxUn=0.0,  qxDn=0.0;

	dpreal qyLe=0.0,  qyRe=0.0;
	dpreal qyLw=0.0,  qyRw=0.0;
	dpreal qyUs=0.0,  qyDs=0.0;
	dpreal qyUn=0.0,  qyDn=0.0;

	dpreal uLe=0.0,   uRe=0.0;									 
	dpreal uLw=0.0,   uRw=0.0;									 
	dpreal uDs=0.0,   uUs=0.0;
	dpreal uDn=0.0,   uUn=0.0;	

	dpreal vLe=0.0,   vRe=0.0;		 
	dpreal vLw=0.0,   vRw=0.0;
	dpreal vDn=0.0,   vUn=0.0;
	dpreal vDs=0.0,   vUs=0.0;

	dpreal hLe=0.0,   hRe=0.0;
	dpreal hLw=0.0,   hRw=0.0;
	dpreal hDs=0.0,   hUs=0.0;
	dpreal hDn=0.0,   hUn=0.0;

	dpreal etaLe=0.0, etaRe=0.0;
	dpreal etaLw=0.0, etaRw=0.0;
	dpreal etaDs=0.0, etaUs=0.0;
	dpreal etaDn=0.0, etaUn=0.0;

	dpreal zbLe=0.0,  zbRe=0.0;
	dpreal zbLw=0.0,  zbRw=0.0;
	dpreal zbDs=0.0,  zbUs=0.0;
	dpreal zbDn=0.0,  zbUn=0.0;

	dpreal zbe=0.0, zbw=0.0, zbs=0.0, zbn=0.0;

	dpreal ustar=0.0, vstar=0.0, hstar=0.0, SL=0.0, SR=0.0, SU=0.0, SD=0.0, Sstar=0.0;
	Vec(dpreal) FL(3), FR(3), GD(3), GU(3), UL(3), UR(3), UD(3), UU(3), FstarL(3), FstarR(3), GstarD(3), GstarU(3);


	// variables for wetting-drying
	dpreal etaLeWD=0.0, etaReWD=0.0;
	dpreal etaLwWD=0.0, etaRwWD=0.0;
	dpreal etaDsWD=0.0, etaUsWD=0.0;
	dpreal etaDnWD=0.0, etaUnWD=0.0;

	dpreal hLeWD=0.0, hReWD=0.0;
	dpreal hLwWD=0.0, hRwWD=0.0;
	dpreal hDsWD=0.0, hUsWD=0.0;
	dpreal hDnWD=0.0, hUnWD=0.0;

	dpreal qxLeWD=0.0, qxReWD=0.0;
	dpreal qxLwWD=0.0, qxRwWD=0.0;
	dpreal qxDsWD=0.0, qxUsWD=0.0;
	dpreal qxDnWD=0.0, qxUnWD=0.0;

	dpreal qyLeWD=0.0, qyReWD=0.0;
	dpreal qyLwWD=0.0, qyRwWD=0.0;
	dpreal qyDsWD=0.0, qyUsWD=0.0;
	dpreal qyDnWD=0.0, qyUnWD=0.0;

	dpreal dz=0.0;


	// fill arrays and vectors with zero
	U.fill(0.0); Up.fill(0.0); Fe.fill(0.0); Fw.fill(0.0); Gs.fill(0.0); Gn.fill(0.0); S.fill(0.0); 	
	FL.fill(0.0); FR.fill(0.0); GD.fill(0.0); GU.fill(0.0); UL.fill(0.0); UR.fill(0.0); UD.fill(0.0); UU.fill(0.0);	


	// determine width and precision for time label in result file names and console output
	timeLabel(tstop, dt, outprf, outpc, wwfile, wwcons, precfile, preccons);


	// if there are no hydrographs (only open and closed boundaries), currHydr is simply a vector holding a 1, such that it can be used 
	// when looping through ghost cells (otherwise an error is produced from trying to get element zero of this vector)
	if(ncolHydrgArr == 0)
	{
		currHydr.redim(1);
		currHydr.fill(1.0);
	}


	// time loop
	// ---------------------------------------------------------------------------------------
	// ---------------------------------------------------------------------------------------
	while(t <= (tstop + outprf))
	{
		// write current run information to console
		if(t >= toutc) 
		{
			printf("current time = % *.*f, Cmax = % *.*e, dt = % *.*e\n", 
				wwcons, preccons, t,
				wwcons, preccons, Cmax,
				wwcons, preccons, dt);
			toutc = toutc + outpc;
		}

		// write state variables to solution file		
		if (t >= toutf)
		{
			tString = oform("%0*.*f", wwfile, precfile, t); // make time stamp

			if(uout == "yes") 
				writeScal(uname, outdir, tString, rfpre, fileext="u", i, j, dx, nx, ny, xll, yll, noData, h, u);			

			if(vout == "yes")
				writeScal(vname, outdir, tString, rfpre, fileext="v", i, j, dx, nx, ny, xll, yll, noData, h, v);			

			if(etaout == "yes")
				writeScal(etaname, outdir, tString, rfpre, fileext="eta", i, j, dx, nx, ny, xll, yll, noData, h, eta);			

			if(hout == "yes")
				writeScal(hname, outdir, tString, rfpre, fileext="h", i, j, dx, nx, ny, xll, yll, noData, h, h);			

			if(vmagout == "yes")
				writeVmag(vmagname, outdir, tString, rfpre, fileext="vmag", i, j, dx, nx, ny, xll, yll, noData, h, u, v);

			if(flowdout == "yes")
				writeFlowDir(flowdname, outdir, tString, rfpre, fileext="fdir", i, j, dx, nx, ny, xll, yll, noData, h, v, u);

			toutf = toutf + outprf;
		}


		// interpolate BC hydrographs to current time step
		// --------------------------------------------------------------------------------------		
		hydrgInterp(hydrgArr, nrowHydrgArr, ncolHydrgArr, hydrgNoSteps, currHydr, slopehydrg, t);



		// RK step 1
		// ---------------------------------------------------------------------------------------
		// ---------------------------------------------------------------------------------------


		// loop through ghost cells
		// ---------------------------------------------------------------------------------------		

		// E edge	
		i=nx+2;
		for(j=2; j<=(ny+1); j++)
		{
			indCurrHydr = int(bcnE(j) + fope(j));
			v(j,i)   = v(j,i-1); // slip boundary condition
			u(j,i)   = fube(j)*(fope(j)*u(j,i-1) + fetae(j)*u(j,i-1) + fqe(j)*currHydr(indCurrHydr)*flowDistrE(j)) + (1.0-fube(j))*(-u(j,i-1));
			eta(j,i) = fube(j)*(fope(j)*eta(j,i-1) + fetae(j)*currHydr(indCurrHydr) + fqe(j)*eta(j,i-1)) + (1.0-fube(j))*(eta(j,i-1));
			h(j,i)   = eta(j,i) - zb(j,i);
			qx(j,i)  = u(j,i)*h(j,i);
			qy(j,i)  = v(j,i)*h(j,i);
		}

		// W edge
		i=1;
		for(j=2; j<=(ny+1); j++)
		{
			indCurrHydr = int(bcnW(j) + fopw(j));
			v(j,i)   = v(j,i+1); // slip boundary condition
			u(j,i)   = fubw(j)*(fopw(j)*u(j,i+1) + fetaw(j)*u(j,i+1) + fqw(j)*currHydr(indCurrHydr)*flowDistrW(j)) + (1.0-fubw(j))*(-u(j,i+1));
			eta(j,i) = fubw(j)*(fopw(j)*eta(j,i+1) + fetaw(j)*currHydr(indCurrHydr) + fqw(j)*eta(j,i+1)) + (1.0-fubw(j))*(eta(j,i+1));
			h(j,i)   = eta(j,i) - zb(j,i);
			qx(j,i)  = u(j,i)*h(j,i);
			qy(j,i)  = v(j,i)*h(j,i);
		}

		// S edge
		j=ny+2;
		for(i=2; i<=(nx+1); i++)
		{
			indCurrHydr = int(bcnS(i) + fops(i));
			v(j,i)   = fubs(i)*(fops(i)*v(j-1,i) + fetas(i)*v(j-1,i) + fqs(i)*currHydr(indCurrHydr)*flowDistrS(i)) + (1.0-fubs(i))*(-v(j-1,i));
			u(j,i)   = u(j-1,i); // slip boundary condition
			eta(j,i) = fubs(i)*(fops(i)*eta(j-1,i) + fetas(i)*currHydr(indCurrHydr) + fqs(i)*eta(j-1,i)) + (1.0-fubs(i))*(eta(j-1,i));
			h(j,i)   = eta(j,i) - zb(j,i);		
			qx(j,i)  = u(j,i)*h(j,i);
			qy(j,i)  = v(j,i)*h(j,i);
		}

		// N edge
		j=1;
		for(i=2; i<=(nx+1); i++)
		{
			indCurrHydr = int(bcnN(i) + fopn(i));
			v(j,i)   = fubn(i)*(fopn(i)*v(j+1,i) + fetan(i)*v(j+1,i) + fqn(i)*currHydr(indCurrHydr)*flowDistrN(i)) + (1.0-fubn(i))*(-v(j+1,i));
			u(j,i)   = u(j+1,i); // slip boundary condition
			eta(j,i) = fubn(i)*(fopn(i)*eta(j+1,i) + fetan(i)*currHydr(indCurrHydr) + fqn(i)*eta(j+1,i)) + (1.0-fubn(i))*(eta(j+1,i));
			h(j,i)   = eta(j,i) - zb(j,i);
			qx(j,i)  = u(j,i)*h(j,i);
			qy(j,i)  = v(j,i)*h(j,i);
		}


		// loop through internal cells
		// ---------------------------------------------------------------------------------------		

		for (j=2; j<=(ny+1); j++)
		{
			for (i=2; i<=(nx+1); i++)
			{

				// implicit friction term
				// ----------------------------------------------------------------------------------
				U(2) = qx(j,i);
				U(3) = qy(j,i);

				if(h(j,i) > 1e-10)
				{
					Cf  = (9.81*sqr(man(j,i)))/(pow(h(j,i),0.33));
					Sfx = Cf*u(j,i)*sqrt(sqr(u(j,i)) + sqr(v(j,i)));
					Sfy = Cf*v(j,i)*sqrt(sqr(u(j,i)) + sqr(v(j,i)));
					Dx  = 1.0 + 2.0*dt*Cf*fabs(U(2))/(sqr(h(j,i)));
					Dy  = 1.0 + 2.0*dt*Cf*fabs(U(3))/(sqr(h(j,i)));
					Fx  = Sfx/Dx;
					Fy  = Sfy/Dy;
				}
				else
				{
					Cf  = 0.0;
					Sfx = 0.0;
					Sfy = 0.0;
					Dx  = 1.0;
					Dy  = 1.0;
					Fx  = 0.0;
					Fy  = 0.0;
				}

				// update fluxes
				Up(2) = U(2) + 0.5*dt*Fx;
				Up(3) = U(3) + 0.5*dt*Fy;

				// limit friction force
				if(Up(2)*U(2) < 0.0) Up(2) = U(2);
				if(Up(3)*U(3) < 0.0) Up(3) = U(3);

				// compute Sf(t+1)
				Sfxp = (Up(2) - U(2))/(0.5*dt);
				Sfyp = (Up(3) - U(3))/(0.5*dt);


				// reconstruct face values
				// 2nd order for wet cells with all wet neighbors, 1st order otherwise
				// ----------------------------------------------------------------------------------

				if(h(j,i) <= 1e-10) // dry cell? -> 1st order reconstruction, zero inner velocities and fluxes
				{				
					recon1dry(eta, h, qx, qy, u, v, zb,
						etaLe, etaRe, etaLw, etaRw, 
						etaDn, etaUn, etaDs, etaUs,
						hLe, hRe, hLw, hRw, 
						hDn, hUn, hDs, hUs,
						qxLe, qxRe, qxLw, qxRw,
						qxDn, qxUn, qxDs, qxUs, 
						qyLe, qyRe, qyLw, qyRw, 
						qyDn, qyUn, qyDs, qyUs, 
						uLe, uRe, uLw, uRw,
						uDn, uUn, uDs, uUs,
						vLe, vRe, vLw, vRw, 
						vDn, vUn, vDs, vUs,
						zbLe, zbRe, zbLw, zbRw,
						zbDn, zbUn, zbDs, zbUs,
						i,j
						);
				}
				else // wet cell
				{
					// check h in neighbors
					neighh(1)=h(j,i+1); neighh(2)=h(j,i-1); neighh(3)=h(j-1,i); neighh(4)=h(j+1,i);

					if(doubVectMin(neighh) <= 1e-10) // wet cell with at least one dry neighbor? -> 1st order reconstruction
					{
						recon1dryWet(eta, h, qx, qy, u, v, zb,
							etaLe, etaRe, etaLw, etaRw, 
							etaDn, etaUn, etaDs, etaUs,
							hLe, hRe, hLw, hRw, 
							hDn, hUn, hDs, hUs,
							qxLe, qxRe, qxLw, qxRw,
							qxDn, qxUn, qxDs, qxUs, 
							qyLe, qyRe, qyLw, qyRw, 
							qyDn, qyUn, qyDs, qyUs, 
							uLe, uRe, uLw, uRw,
							uDn, uUn, uDs, uUs,
							vLe, vRe, vLw, vRw, 
							vDn, vUn, vDs, vUs,
							zbLe, zbRe, zbLw, zbRw,
							zbDn, zbUn, zbDs, zbUs,
							i,j);
					}
					else // wet cell with wet neighbors -> 2nd order
					{
						// inner faces
						recon2inn(eta, etaLe, etaRw, etaDn, etaUs,
							inVect, limslope,
							dx, dy, i, j);
						recon2inn(h, hLe, hRw, hDn, hUs,
							inVect, limslope,
							dx, dy, i, j);
						recon2inn(qx, qxLe, qxRw, qxDn, qxUs,
							inVect, limslope,
							dx, dy, i, j);
						recon2inn(qy, qyLe, qyRw, qyDn, qyUs,
							inVect, limslope,
							dx, dy, i, j);

						zbLe  = etaLe - hLe;
						zbRw  = etaRw - hRw;
						zbDn  = etaDn - hDn;
						zbUs  = etaUs - hUs;

						// outer faces
						recon2out(eta, etaRe, etaLw, etaUn, etaDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);
						recon2out(h, hRe, hLw, hUn, hDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);
						recon2out(qx, qxRe, qxLw, qxUn, qxDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);
						recon2out(qy, qyRe, qyLw, qyUn, qyDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);

						zbRe = etaRe - hRe;
						zbLw = etaLw - hLw;
						zbUn = etaUn - hUn;
						zbDs = etaDs - hDs;

						// velocities
						uLe = qxLe/hLe; uRw = qxRw/hRw; uDn = qxDn/hDn; uUs = qxUs/hUs;
						uRe = qxRe/hRe; uLw = qxLw/hLw; uUn = qxUn/hUn; uDs = qxDs/hDs;

						vLe = qyLe/hLe; vRw = qyRw/hRw; vDn = qyDn/hDn; vUs = qyUs/hUs;
						vRe = qyRe/hRe; vLw = qyLw/hLw; vUn = qyUn/hUn; vDs = qyDs/hDs;
					}
				}


				// wetting-drying (Liang 2010) & HLLC Riemann fluxes (Toro 2001)
				// -------------------------------------------------------------------------

				// east face			
				wetDry(etaLe, etaRe, etaLeWD, etaReWD,
					hLe, hRe, hLeWD, hReWD,
					qxLe, qxRe, qxLeWD, qxReWD, 
					qyLe, qyRe, qyLeWD, qyReWD, 
					uLe, uRe, vLe, vRe, 
					zbLe, zbRe, 
					zbe, dz
					);
				xFlux(Fe,
					uLe, uRe, 
					vLe, vRe,
					qxLeWD, qxReWD,
					qyLeWD, qyReWD,
					hLeWD, hReWD, 
					etaLeWD, etaReWD,
					zbe,
					UL, UR, 
					FL, FR, 
					ustar, hstar, 
					SL, SR, Sstar, 
					FstarL, FstarR
					);

				// west face
				wetDry(etaRw, etaLw, etaRwWD, etaLwWD,
					hRw, hLw, hRwWD, hLwWD,
					qxRw, qxLw, qxRwWD, qxLwWD, 
					qyRw, qyLw, qyRwWD, qyLwWD, 
					uRw, uLw, vRw, vLw, 
					zbRw, zbLw, 
					zbw, dz
					);
				xFlux(Fw,
					uLw, uRw, 
					vLw, vRw,
					qxLwWD, qxRwWD,
					qyLwWD, qyRwWD,
					hLwWD, hRwWD, 
					etaLwWD, etaRwWD,
					zbw,
					UL, UR, 
					FL, FR, 
					ustar, hstar, 
					SL, SR, Sstar, 
					FstarL, FstarR
					);

				// north face
				wetDry(etaDn, etaUn, etaDnWD, etaUnWD,
					hDn, hUn, hDnWD, hUnWD,
					qxDn, qxUn, qxDnWD, qxUnWD, 
					qyDn, qyUn, qyDnWD, qyUnWD, 
					uDn, uUn, vDn, vUn, 
					zbDn, zbUn, 
					zbn, dz
					);
				yFlux(Gn, 
					uDn, uUn, 
					vDn, vUn,
					qxDnWD, qxUnWD,
					qyDnWD, qyUnWD,
					hDnWD, hUnWD, 
					etaDnWD, etaUnWD,
					zbn,
					UD, UU, 
					GD, GU, 
					vstar, hstar, 
					SD, SU, Sstar, 
					GstarD, GstarU
					);

				// south face
				wetDry(etaUs, etaDs, etaUsWD, etaDsWD,
					hUs, hDs, hUsWD, hDsWD,
					qxUs, qxDs, qxUsWD, qxDsWD, 
					qyUs, qyDs, qyUsWD, qyDsWD, 
					uUs, uDs, vUs, vDs, 
					zbUs, zbDs, 
					zbs, dz
					);
				yFlux(Gs, 
					uDs, uUs,
					vDs, vUs,
					qxDsWD, qxUsWD,
					qyDsWD, qyUsWD,
					hDsWD, hUsWD, 
					etaDsWD, etaUsWD,
					zbs,
					UD, UU, 
					GD, GU, 
					vstar, hstar, 
					SD, SU, Sstar, 
					GstarD, GstarU
					);


				// Update cell
				// --------------------------------------------

				// get current values
				U(1) = eta(j,i);
				U(2) = qx(j,i);
				U(3) = qy(j,i);

				// compute source terms
				S(1) = 0.0;
				S(2) = -9.81*0.5*(etaRwWD + etaLeWD)*(zbe-zbw)/dx - Sfxp;
				S(3) = -9.81*0.5*(etaDnWD + etaUsWD)*(zbn-zbs)/dy - Sfyp;

				// Runge-Kutta coefficient
				K11(j,i) = -(Fe(1) - Fw(1))/dx - (Gn(1) - Gs(1))/dy + S(1);
				K12(j,i) = -(Fe(2) - Fw(2))/dx - (Gn(2) - Gs(2))/dy + S(2);
				K13(j,i) = -(Fe(3) - Fw(3))/dx - (Gn(3) - Gs(3))/dy + S(3);

				// compute intermediate quantities
				Up(1) = U(1) + 0.5*dt*K11(j,i);
				Up(2) = U(2) + 0.5*dt*K12(j,i);
				Up(3) = U(3) + 0.5*dt*K13(j,i);

				// compute intermediate quantities
				eta1(j,i) = Up(1);
				qx1(j,i)  = Up(2);
				qy1(j,i)  = Up(3);
				h1(j,i)   = eta1(j,i) - zb(j,i);
				if(h1(j,i) <= 1e-10) {u1(j,i) = 0.0; v1(j,i) = 0.0; qx1(j,i) = 0.0; qy1(j,i) = 0.0; h1(j,i)=0.0; eta1(j,i)=zb(j,i);}
				else {u1(j,i) = qx1(j,i)/h1(j,i); v1(j,i) = qy1(j,i)/h1(j,i);}
			}
		}


		// RK step 2
		// ---------------------------------------------------------------------------------------
		// ---------------------------------------------------------------------------------------

		// ghost cells
		// --------------------------------------------

		// E edge	
		i=nx+2;
		for(j=2; j<=(ny+1); j++)
		{
			indCurrHydr = int(bcnE(j) + fope(j));
			v1(j,i)   = v1(j,i-1); // slip boundary condition
			u1(j,i)   = fube(j)*(fope(j)*u1(j,i-1) + fetae(j)*u1(j,i-1) + fqe(j)*currHydr(indCurrHydr)*flowDistrE(j)) + (1.0-fube(j))*(-u1(j,i-1));
			eta1(j,i) = fube(j)*(fope(j)*eta1(j,i-1) + fetae(j)*currHydr(indCurrHydr) + fqe(j)*eta1(j,i-1)) + (1.0-fube(j))*(eta1(j,i-1));
			h1(j,i)   = eta1(j,i) - zb(j,i);
			qx1(j,i)  = u1(j,i)*h1(j,i);
			qy1(j,i)  = v1(j,i)*h1(j,i);
		}

		// W edge
		i=1;
		for(j=2; j<=(ny+1); j++)
		{
			indCurrHydr = int(bcnW(j) + fopw(j));
			v1(j,i)   = v1(j,i+1); // slip boundary condition
			u1(j,i)   = fubw(j)*(fopw(j)*u1(j,i+1) + fetaw(j)*u1(j,i+1) + fqw(j)*currHydr(indCurrHydr)*flowDistrW(j)) + (1.0-fubw(j))*(-u1(j,i+1));
			eta1(j,i) = fubw(j)*(fopw(j)*eta1(j,i+1) + fetaw(j)*currHydr(indCurrHydr) + fqw(j)*eta1(j,i+1)) + (1.0-fubw(j))*(eta1(j,i+1));
			h1(j,i)   = eta1(j,i) - zb(j,i);
			qx1(j,i)  = u1(j,i)*h1(j,i);
			qy1(j,i)  = v1(j,i)*h1(j,i);
		}

		// S edge
		// --------------------------------------------
		j=ny+2;
		for(i=2; i<=(nx+1); i++)
		{
			indCurrHydr = int(bcnS(i) + fops(i));
			v1(j,i)   = fubs(i)*(fops(i)*v1(j-1,i) + fetas(i)*v1(j-1,i) + fqs(i)*currHydr(indCurrHydr)*flowDistrS(i)) + (1.0-fubs(i))*(-v1(j-1,i));
			u1(j,i)   = u1(j-1,i); // slip boundary condition
			eta1(j,i) = fubs(i)*(fops(i)*eta1(j-1,i) + fetas(i)*currHydr(indCurrHydr) + fqs(i)*eta1(j-1,i)) + (1.0-fubs(i))*(eta1(j-1,i));
			h1(j,i)   = eta1(j,i) - zb(j,i);
			qx1(j,i)  = u1(j,i)*h1(j,i);
			qy1(j,i)  = v1(j,i)*h1(j,i);
		}

		// N edge
		// --------------------------------------------
		j=1;
		for(i=2; i<=(nx+1); i++)
		{
			indCurrHydr = int(bcnN(i) + fopn(i));
			v1(j,i)   = fubn(i)*(fopn(i)*v1(j+1,i) + fetan(i)*v1(j+1,i) + fqn(i)*currHydr(indCurrHydr)*flowDistrN(i)) + (1.0-fubn(i))*(-v1(j+1,i));
			u1(j,i)   = u1(j+1,i); // slip boundary condition
			eta1(j,i) = fubn(i)*(fopn(i)*eta1(j+1,i) + fetan(i)*currHydr(indCurrHydr) + fqn(i)*eta1(j+1,i)) + (1.0-fubn(i))*(eta1(j+1,i));
			h1(j,i)   = eta1(j,i) - zb(j,i);
			qx1(j,i)  = u1(j,i)*h1(j,i);
			qy1(j,i)  = v1(j,i)*h1(j,i);
		}


		// loop through internal cells
		// ---------------------------------------------------------------------------------------
		// ---------------------------------------------------------------------------------------

		for (j=2; j<=(ny+1); j++)
		{
			for (i=2; i<=(nx+1); i++)
			{

				// implicit friction term
				// ----------------------------------------------------------------------------------
				U(2) = qx1(j,i);
				U(3) = qy1(j,i);

				if(h1(j,i) > 1e-10)
				{
					Cf  = (9.81*sqr(man(j,i)))/(pow(h1(j,i),0.33));
					Sfx = Cf*u1(j,i)*sqrt(sqr(u1(j,i)) + sqr(v1(j,i)));
					Sfy = Cf*v1(j,i)*sqrt(sqr(u1(j,i)) + sqr(v1(j,i)));
					Dx  = 1.0 + 2.0*dt*Cf*fabs(U(2))/(sqr(h1(j,i)));
					Dy  = 1.0 + 2.0*dt*Cf*fabs(U(3))/(sqr(h1(j,i)));
					Fx  = Sfx/Dx;
					Fy  = Sfy/Dy;
				}
				else
				{
					Cf  = 0.0;
					Sfx = 0.0;
					Sfy = 0.0;
					Dx  = 1.0;
					Dy  = 1.0;
					Fx  = 0.0;
					Fy  = 0.0;
				}

				// update fluxes
				Up(2) = U(2) + 0.5*dt*Fx;
				Up(3) = U(3) + 0.5*dt*Fy;

				// limit friction force
				if(Up(2)*U(2) < 0.0) Up(2) = U(2);
				if(Up(3)*U(3) < 0.0) Up(3) = U(3);

				// compute Sf(t+1)
				Sfxp = (Up(2) - U(2))/(0.5*dt);
				Sfyp = (Up(3) - U(3))/(0.5*dt);


				// reconstruct face values
				// 2nd order for wet cells with all wet neighbors, 1st order otherwise
				// ----------------------------------------------------------------------------------

				if(h1(j,i) <= 1e-10) // dry cell? -> 1st order reconstruction, zero inner depth, velocities and fluxes
				{					
					recon1dry(eta1, h1, qx1, qy1, u1, v1, zb,
						etaLe, etaRe, etaLw, etaRw, 
						etaDn, etaUn, etaDs, etaUs,
						hLe, hRe, hLw, hRw, 
						hDn, hUn, hDs, hUs,
						qxLe, qxRe, qxLw, qxRw,
						qxDn, qxUn, qxDs, qxUs, 
						qyLe, qyRe, qyLw, qyRw, 
						qyDn, qyUn, qyDs, qyUs, 
						uLe, uRe, uLw, uRw,
						uDn, uUn, uDs, uUs,
						vLe, vRe, vLw, vRw, 
						vDn, vUn, vDs, vUs,
						zbLe, zbRe, zbLw, zbRw,
						zbDn, zbUn, zbDs, zbUs,
						i,j
						);
				}
				else // wet cell
				{
					// check h in neighbors
					neighh(1)=h1(j,i+1); neighh(2)=h1(j,i-1); neighh(3)=h1(j-1,i); neighh(4)=h1(j+1,i);

					if(doubVectMin(neighh) <= 1e-10) // wet cell with at least one dry neighbor? -> 1st order reconstruction
					{
						recon1dryWet(eta1, h1, qx1, qy1, u1, v1, zb,
							etaLe, etaRe, etaLw, etaRw, 
							etaDn, etaUn, etaDs, etaUs,
							hLe, hRe, hLw, hRw, 
							hDn, hUn, hDs, hUs,
							qxLe, qxRe, qxLw, qxRw,
							qxDn, qxUn, qxDs, qxUs, 
							qyLe, qyRe, qyLw, qyRw, 
							qyDn, qyUn, qyDs, qyUs, 
							uLe, uRe, uLw, uRw,
							uDn, uUn, uDs, uUs,
							vLe, vRe, vLw, vRw, 
							vDn, vUn, vDs, vUs,
							zbLe, zbRe, zbLw, zbRw,
							zbDn, zbUn, zbDs, zbUs,
							i,j);
					}
					else // wet cell with wet neighbors -> 2nd order
					{
						// inner faces
						recon2inn(eta1, etaLe, etaRw, etaDn, etaUs,
							inVect, limslope,
							dx, dy, i, j);
						recon2inn(h1, hLe, hRw, hDn, hUs,
							inVect, limslope,
							dx, dy, i, j);
						recon2inn(qx1, qxLe, qxRw, qxDn, qxUs,
							inVect, limslope,
							dx, dy, i, j);
						recon2inn(qy1, qyLe, qyRw, qyDn, qyUs,
							inVect, limslope,
							dx, dy, i, j);

						zbLe  = etaLe - hLe;
						zbRw  = etaRw - hRw;
						zbDn  = etaDn - hDn;
						zbUs  = etaUs - hUs;

						// outer faces
						recon2out(eta1, etaRe, etaLw, etaUn, etaDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);
						recon2out(h1, hRe, hLw, hUn, hDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);
						recon2out(qx1, qxRe, qxLw, qxUn, qxDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);
						recon2out(qy1, qyRe, qyLw, qyUn, qyDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);

						zbRe = etaRe - hRe;
						zbLw = etaLw - hLw;
						zbUn = etaUn - hUn;
						zbDs = etaDs - hDs;

						// velocities
						uLe = qxLe/hLe; uRw = qxRw/hRw; uDn = qxDn/hDn; uUs = qxUs/hUs;
						uRe = qxRe/hRe; uLw = qxLw/hLw; uUn = qxUn/hUn; uDs = qxDs/hDs;

						vLe = qyLe/hLe; vRw = qyRw/hRw; vDn = qyDn/hDn; vUs = qyUs/hUs;
						vRe = qyRe/hRe; vLw = qyLw/hLw; vUn = qyUn/hUn; vDs = qyDs/hDs;
					}
				}


				// wetting-drying (Liang 2010) & HLLC Riemann fluxes (Toro 2001)
				// -------------------------------------------------------------------------

				// east face							
				wetDry(etaLe, etaRe, etaLeWD, etaReWD,
					hLe, hRe, hLeWD, hReWD,
					qxLe, qxRe, qxLeWD, qxReWD, 
					qyLe, qyRe, qyLeWD, qyReWD, 
					uLe, uRe, vLe, vRe, 
					zbLe, zbRe, 
					zbe, dz
					);
				xFlux(Fe,
					uLe, uRe, 
					vLe, vRe,
					qxLeWD, qxReWD,
					qyLeWD, qyReWD,
					hLeWD, hReWD, 
					etaLeWD, etaReWD,
					zbe,
					UL, UR, 
					FL, FR, 
					ustar, hstar, 
					SL, SR, Sstar, 
					FstarL, FstarR
					);

				// west face
				wetDry(etaRw, etaLw, etaRwWD, etaLwWD,
					hRw, hLw, hRwWD, hLwWD,
					qxRw, qxLw, qxRwWD, qxLwWD, 
					qyRw, qyLw, qyRwWD, qyLwWD, 
					uRw, uLw, vRw, vLw, 
					zbRw, zbLw, 
					zbw, dz
					);
				xFlux(Fw,
					uLw, uRw, 
					vLw, vRw,
					qxLwWD, qxRwWD,
					qyLwWD, qyRwWD,
					hLwWD, hRwWD, 
					etaLwWD, etaRwWD,
					zbw,
					UL, UR, 
					FL, FR, 
					ustar, hstar, 
					SL, SR, Sstar, 
					FstarL, FstarR
					);

				// north face
				wetDry(etaDn, etaUn, etaDnWD, etaUnWD,
					hDn, hUn, hDnWD, hUnWD,
					qxDn, qxUn, qxDnWD, qxUnWD, 
					qyDn, qyUn, qyDnWD, qyUnWD, 
					uDn, uUn, vDn, vUn, 
					zbDn, zbUn, 
					zbn, dz
					);
				yFlux(Gn, 
					uDn, uUn, 
					vDn, vUn,
					qxDnWD, qxUnWD,
					qyDnWD, qyUnWD,
					hDnWD, hUnWD, 
					etaDnWD, etaUnWD,
					zbn,
					UD, UU, 
					GD, GU, 
					vstar, hstar, 
					SD, SU, Sstar, 
					GstarD, GstarU
					);

				// south face
				wetDry(etaUs, etaDs, etaUsWD, etaDsWD,
					hUs, hDs, hUsWD, hDsWD,
					qxUs, qxDs, qxUsWD, qxDsWD, 
					qyUs, qyDs, qyUsWD, qyDsWD, 
					uUs, uDs, vUs, vDs, 
					zbUs, zbDs, 
					zbs, dz
					);
				yFlux(Gs, 
					uDs, uUs,
					vDs, vUs,
					qxDsWD, qxUsWD,
					qyDsWD, qyUsWD,
					hDsWD, hUsWD, 
					etaDsWD, etaUsWD,
					zbs,
					UD, UU, 
					GD, GU, 
					vstar, hstar, 
					SD, SU, Sstar, 
					GstarD, GstarU
					);


				// Update cell
				// --------------------------------------------

				// get state variable values at t
				U(1) = eta(j,i);
				U(2) = qx(j,i);
				U(3) = qy(j,i);

				// compute source terms
				S(1) = 0.0;
				S(2) = -9.81*0.5*(etaRwWD + etaLeWD)*(zbe-zbw)/dx - Sfxp;
				S(3) = -9.81*0.5*(etaDnWD + etaUsWD)*(zbn-zbs)/dy - Sfyp;

				// Runge-Kutta coefficient
				K21(j,i) = -(Fe(1) - Fw(1))/dx - (Gn(1) - Gs(1))/dy + S(1);
				K22(j,i) = -(Fe(2) - Fw(2))/dx - (Gn(2) - Gs(2))/dy + S(2);
				K23(j,i) = -(Fe(3) - Fw(3))/dx - (Gn(3) - Gs(3))/dy + S(3);

				// compute intermediate quantities
				Up(1) = U(1) + 0.5*dt*K21(j,i);
				Up(2) = U(2) + 0.5*dt*K22(j,i);
				Up(3) = U(3) + 0.5*dt*K23(j,i);

				eta2(j,i) = Up(1);
				qx2(j,i)  = Up(2);
				qy2(j,i)  = Up(3);
				h2(j,i)   = eta2(j,i) - zb(j,i);
				if(h2(j,i) <= 1e-10) {u2(j,i) = 0.0; v2(j,i) = 0.0;  qx2(j,i) = 0.0; qy2(j,i) = 0.0; h2(j,i)=0.0; eta2(j,i)=zb(j,i);}
				else {u2(j,i) = qx2(j,i)/h2(j,i); v2(j,i) = qy2(j,i)/h2(j,i);}
			}
		}


		// RK step 3
		// ---------------------------------------------------------------------------------------
		// ---------------------------------------------------------------------------------------

		// ghost cells
		// --------------------------------------------

		// E edge	
		i=nx+2;
		for(j=2; j<=(ny+1); j++)
		{
			indCurrHydr = int(bcnE(j) + fope(j));
			v2(j,i)   = v2(j,i-1); // slip boundary condition
			u2(j,i)   = fube(j)*(fope(j)*u2(j,i-1) + fetae(j)*u2(j,i-1) + fqe(j)*currHydr(indCurrHydr)*flowDistrE(j)) + (1.0-fube(j))*(-u2(j,i-1));
			eta2(j,i) = fube(j)*(fope(j)*eta2(j,i-1) + fetae(j)*currHydr(indCurrHydr) + fqe(j)*eta2(j,i-1)) + (1.0-fube(j))*(eta2(j,i-1));
			h2(j,i)   = eta2(j,i) - zb(j,i);
			qx2(j,i)  = u2(j,i)*h2(j,i);
			qy2(j,i)  = v2(j,i)*h2(j,i);
		}

		// W edge
		i=1;
		for(j=2; j<=(ny+1); j++)
		{
			indCurrHydr = int(bcnW(j) + fopw(j));
			v2(j,i)   = v2(j,i+1); // slip boundary condition
			u2(j,i)   = fubw(j)*(fopw(j)*u2(j,i+1) + fetaw(j)*u2(j,i+1) + fqw(j)*currHydr(indCurrHydr)*flowDistrW(j)) + (1.0-fubw(j))*(-u2(j,i+1));
			eta2(j,i) = fubw(j)*(fopw(j)*eta2(j,i+1) + fetaw(j)*currHydr(indCurrHydr) + fqw(j)*eta2(j,i+1)) + (1.0-fubw(j))*(eta2(j,i+1));
			h2(j,i)   = eta2(j,i) - zb(j,i);
			qx2(j,i)  = u2(j,i)*h2(j,i);
			qy2(j,i)  = v2(j,i)*h2(j,i);
		}

		// S edge
		// --------------------------------------------
		j=ny+2;
		for(i=2; i<=(nx+1); i++)
		{
			indCurrHydr = int(bcnS(i) + fops(i));
			v2(j,i)   = fubs(i)*(fops(i)*v2(j-1,i) + fetas(i)*v2(j-1,i) + fqs(i)*currHydr(indCurrHydr)*flowDistrS(i)) + (1.0-fubs(i))*(-v2(j-1,i));
			u2(j,i)   = u2(j-1,i); // slip boundary condition
			eta2(j,i) = fubs(i)*(fops(i)*eta2(j-1,i) + fetas(i)*currHydr(indCurrHydr) + fqs(i)*eta2(j-1,i)) + (1.0-fubs(i))*(eta2(j-1,i));
			h2(j,i)   = eta2(j,i) - zb(j,i);
			qx2(j,i)  = u2(j,i)*h2(j,i);
			qy2(j,i)  = v2(j,i)*h2(j,i);
		}

		// N edge
		// --------------------------------------------
		j=1;
		for(i=2; i<=(nx+1); i++)
		{
			indCurrHydr = int(bcnN(i) + fopn(i));
			v2(j,i)   = fubn(i)*(fopn(i)*v2(j+1,i) + fetan(i)*v2(j+1,i) + fqn(i)*currHydr(indCurrHydr)*flowDistrN(i)) + (1.0-fubn(i))*(-v2(j+1,i));
			u2(j,i)   = u2(j+1,i); // slip boundary condition
			eta2(j,i) = fubn(i)*(fopn(i)*eta2(j+1,i) + fetan(i)*currHydr(indCurrHydr) + fqn(i)*eta2(j+1,i)) + (1.0-fubn(i))*(eta2(j+1,i));
			h2(j,i)   = eta2(j,i) - zb(j,i);
			qx2(j,i)  = u2(j,i)*h2(j,i);
			qy2(j,i)  = v2(j,i)*h2(j,i);
		}


		// loop through internal cells
		// ---------------------------------------------------------------------------------------
		// ---------------------------------------------------------------------------------------

		for (j=2; j<=(ny+1); j++)
		{
			for (i=2; i<=(nx+1); i++)
			{

				// implicit friction term
				// ----------------------------------------------------------------------------------
				U(2) = qx2(j,i);
				U(3) = qy2(j,i);

				if(h2(j,i) > 1e-10)
				{
					Cf  = (9.81*sqr(man(j,i)))/(pow(h2(j,i),0.33));
					Sfx = Cf*u2(j,i)*sqrt(sqr(u2(j,i)) + sqr(v2(j,i)));
					Sfy = Cf*v2(j,i)*sqrt(sqr(u2(j,i)) + sqr(v2(j,i)));
					Dx  = 1.0 + 2.0*dt*Cf*fabs(U(2))/(sqr(h2(j,i)));
					Dy  = 1.0 + 2.0*dt*Cf*fabs(U(3))/(sqr(h2(j,i)));
					Fx  = Sfx/Dx;
					Fy  = Sfy/Dy;
				}
				else
				{
					Cf  = 0.0;
					Sfx = 0.0;
					Sfy = 0.0;
					Dx  = 1.0;
					Dy  = 1.0;
					Fx  = 0.0;
					Fy  = 0.0;
				}

				// update fluxes
				Up(2) = U(2) + dt*Fx;
				Up(3) = U(3) + dt*Fy;

				// limit friction force
				if(Up(2)*U(2) < 0.0) Up(2) = U(2);
				if(Up(3)*U(3) < 0.0) Up(3) = U(3);

				// compute Sf(t+1)
				Sfxp = (Up(2) - U(2))/dt;
				Sfyp = (Up(3) - U(3))/dt;


				// reconstruct face values
				// 2nd order for wet cells with all wet neighbors, 1st order otherwise
				// ----------------------------------------------------------------------------------

				if(h2(j,i) <= 1e-10) // dry cell? -> 1st order reconstruction, zero inner depth, velocities and fluxes
				{					
					recon1dry(eta2, h2, qx2, qy2, u2, v2, zb,
						etaLe, etaRe, etaLw, etaRw, 
						etaDn, etaUn, etaDs, etaUs,
						hLe, hRe, hLw, hRw, 
						hDn, hUn, hDs, hUs,
						qxLe, qxRe, qxLw, qxRw,
						qxDn, qxUn, qxDs, qxUs, 
						qyLe, qyRe, qyLw, qyRw, 
						qyDn, qyUn, qyDs, qyUs, 
						uLe, uRe, uLw, uRw,
						uDn, uUn, uDs, uUs,
						vLe, vRe, vLw, vRw, 
						vDn, vUn, vDs, vUs,
						zbLe, zbRe, zbLw, zbRw,
						zbDn, zbUn, zbDs, zbUs,
						i,j
						);
				}
				else // wet cell
				{
					// check h in neighbors
					neighh(1)=h2(j,i+1); neighh(2)=h2(j,i-1); neighh(3)=h2(j-1,i); neighh(4)=h2(j+1,i);

					if(doubVectMin(neighh) <= 1e-10) // wet cell with at least one dry neighbor? -> 1st order reconstruction
					{
						recon1dryWet(eta2, h2, qx2, qy2, u2, v2, zb,
							etaLe, etaRe, etaLw, etaRw, 
							etaDn, etaUn, etaDs, etaUs,
							hLe, hRe, hLw, hRw, 
							hDn, hUn, hDs, hUs,
							qxLe, qxRe, qxLw, qxRw,
							qxDn, qxUn, qxDs, qxUs, 
							qyLe, qyRe, qyLw, qyRw, 
							qyDn, qyUn, qyDs, qyUs, 
							uLe, uRe, uLw, uRw,
							uDn, uUn, uDs, uUs,
							vLe, vRe, vLw, vRw, 
							vDn, vUn, vDs, vUs,
							zbLe, zbRe, zbLw, zbRw,
							zbDn, zbUn, zbDs, zbUs,
							i,j);
					}
					else // wet cell with wet neighbors -> 2nd order
					{
						// inner faces
						recon2inn(eta2, etaLe, etaRw, etaDn, etaUs,
							inVect, limslope,
							dx, dy, i, j);
						recon2inn(h2, hLe, hRw, hDn, hUs,
							inVect, limslope,
							dx, dy, i, j);
						recon2inn(qx2, qxLe, qxRw, qxDn, qxUs,
							inVect, limslope,
							dx, dy, i, j);
						recon2inn(qy2, qyLe, qyRw, qyDn, qyUs,
							inVect, limslope,
							dx, dy, i, j);

						zbLe  = etaLe - hLe;
						zbRw  = etaRw - hRw;
						zbDn  = etaDn - hDn;
						zbUs  = etaUs - hUs;

						// outer faces
						recon2out(eta2, etaRe, etaLw, etaUn, etaDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);
						recon2out(h2, hRe, hLw, hUn, hDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);
						recon2out(qx2, qxRe, qxLw, qxUn, qxDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);
						recon2out(qy2, qyRe, qyLw, qyUn, qyDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);

						zbRe = etaRe - hRe;
						zbLw = etaLw - hLw;
						zbUn = etaUn - hUn;
						zbDs = etaDs - hDs;

						// velocities
						uLe = qxLe/hLe; uRw = qxRw/hRw; uDn = qxDn/hDn; uUs = qxUs/hUs;
						uRe = qxRe/hRe; uLw = qxLw/hLw; uUn = qxUn/hUn; uDs = qxDs/hDs;

						vLe = qyLe/hLe; vRw = qyRw/hRw; vDn = qyDn/hDn; vUs = qyUs/hUs;
						vRe = qyRe/hRe; vLw = qyLw/hLw; vUn = qyUn/hUn; vDs = qyDs/hDs;
					}
				}


				// wetting-drying (Liang 2010) & HLLC Riemann fluxes (Toro 2001)
				// -------------------------------------------------------------------------

				// east face							
				wetDry(etaLe, etaRe, etaLeWD, etaReWD,
					hLe, hRe, hLeWD, hReWD,
					qxLe, qxRe, qxLeWD, qxReWD, 
					qyLe, qyRe, qyLeWD, qyReWD, 
					uLe, uRe, vLe, vRe, 
					zbLe, zbRe, 
					zbe, dz
					);
				xFlux(Fe,
					uLe, uRe, 
					vLe, vRe,
					qxLeWD, qxReWD,
					qyLeWD, qyReWD,
					hLeWD, hReWD, 
					etaLeWD, etaReWD,
					zbe,
					UL, UR, 
					FL, FR, 
					ustar, hstar, 
					SL, SR, Sstar, 
					FstarL, FstarR
					);

				// west face
				wetDry(etaRw, etaLw, etaRwWD, etaLwWD,
					hRw, hLw, hRwWD, hLwWD,
					qxRw, qxLw, qxRwWD, qxLwWD, 
					qyRw, qyLw, qyRwWD, qyLwWD, 
					uRw, uLw, vRw, vLw, 
					zbRw, zbLw, 
					zbw, dz
					);
				xFlux(Fw,
					uLw, uRw, 
					vLw, vRw,
					qxLwWD, qxRwWD,
					qyLwWD, qyRwWD,
					hLwWD, hRwWD, 
					etaLwWD, etaRwWD,
					zbw,
					UL, UR, 
					FL, FR, 
					ustar, hstar, 
					SL, SR, Sstar, 
					FstarL, FstarR
					);

				// north face
				wetDry(etaDn, etaUn, etaDnWD, etaUnWD,
					hDn, hUn, hDnWD, hUnWD,
					qxDn, qxUn, qxDnWD, qxUnWD, 
					qyDn, qyUn, qyDnWD, qyUnWD, 
					uDn, uUn, vDn, vUn, 
					zbDn, zbUn, 
					zbn, dz
					);
				yFlux(Gn, 
					uDn, uUn, 
					vDn, vUn,
					qxDnWD, qxUnWD,
					qyDnWD, qyUnWD,
					hDnWD, hUnWD, 
					etaDnWD, etaUnWD,
					zbn,
					UD, UU, 
					GD, GU, 
					vstar, hstar, 
					SD, SU, Sstar, 
					GstarD, GstarU
					);

				// south face
				wetDry(etaUs, etaDs, etaUsWD, etaDsWD,
					hUs, hDs, hUsWD, hDsWD,
					qxUs, qxDs, qxUsWD, qxDsWD, 
					qyUs, qyDs, qyUsWD, qyDsWD, 
					uUs, uDs, vUs, vDs, 
					zbUs, zbDs, 
					zbs, dz
					);
				yFlux(Gs, 
					uDs, uUs,
					vDs, vUs,
					qxDsWD, qxUsWD,
					qyDsWD, qyUsWD,
					hDsWD, hUsWD, 
					etaDsWD, etaUsWD,
					zbs,
					UD, UU, 
					GD, GU, 
					vstar, hstar, 
					SD, SU, Sstar, 
					GstarD, GstarU
					);


				// Update cell
				// --------------------------------------------

				// get state variable values at t
				U(1) = eta(j,i);
				U(2) = qx(j,i);
				U(3) = qy(j,i);

				// compute source terms
				S(1) = 0.0;
				S(2) = -9.81*0.5*(etaRwWD + etaLeWD)*(zbe-zbw)/dx - Sfxp;
				S(3) = -9.81*0.5*(etaDnWD + etaUsWD)*(zbn-zbs)/dy - Sfyp;

				// Runge-Kutta coefficient
				K31(j,i) = -(Fe(1) - Fw(1))/dx - (Gn(1) - Gs(1))/dy + S(1);
				K32(j,i) = -(Fe(2) - Fw(2))/dx - (Gn(2) - Gs(2))/dy + S(2);
				K33(j,i) = -(Fe(3) - Fw(3))/dx - (Gn(3) - Gs(3))/dy + S(3);

				// compute intermediate quantities
				Up(1) = U(1) + dt*K31(j,i);
				Up(2) = U(2) + dt*K32(j,i);
				Up(3) = U(3) + dt*K33(j,i);

				eta3(j,i) = Up(1);
				qx3(j,i)  = Up(2);
				qy3(j,i)  = Up(3);
				h3(j,i)   = eta3(j,i) - zb(j,i);
				if(h3(j,i) <= 1e-10) {u3(j,i) = 0.0; v3(j,i) = 0.0; qx3(j,i) = 0.0; qy3(j,i) = 0.0; h3(j,i)=0.0; eta3(j,i)=zb(j,i);}
				else {u3(j,i) = qx3(j,i)/h3(j,i); v3(j,i) = qy3(j,i)/h3(j,i);}
			}
		}


		// RK step 4
		// ---------------------------------------------------------------------------------------
		// ---------------------------------------------------------------------------------------

		// ghost cells
		// --------------------------------------------

		// E edge	
		i=nx+2;
		for(j=2; j<=(ny+1); j++)
		{
			indCurrHydr = int(bcnE(j) + fope(j));
			v3(j,i)   = v3(j,i-1); // slip boundary condition
			u3(j,i)   = fube(j)*(fope(j)*u3(j,i-1) + fetae(j)*u3(j,i-1) + fqe(j)*currHydr(indCurrHydr)*flowDistrE(j)) + (1.0-fube(j))*(-u3(j,i-1));
			eta3(j,i) = fube(j)*(fope(j)*eta3(j,i-1) + fetae(j)*currHydr(indCurrHydr) + fqe(j)*eta3(j,i-1)) + (1.0-fube(j))*(eta3(j,i-1));
			h3(j,i)   = eta3(j,i) - zb(j,i);
			qx3(j,i)  = u3(j,i)*h3(j,i);
			qy3(j,i)  = v3(j,i)*h3(j,i);
		}

		// W edge
		i=1;
		for(j=2; j<=(ny+1); j++)
		{
			indCurrHydr = int(bcnW(j) + fopw(j));
			v3(j,i)   = v3(j,i+1); // slip boundary condition
			u3(j,i)   = fubw(j)*(fopw(j)*u3(j,i+1) + fetaw(j)*u3(j,i+1) + fqw(j)*currHydr(indCurrHydr)*flowDistrW(j)) + (1.0-fubw(j))*(-u3(j,i+1));
			eta3(j,i) = fubw(j)*(fopw(j)*eta3(j,i+1) + fetaw(j)*currHydr(indCurrHydr) + fqw(j)*eta3(j,i+1)) + (1.0-fubw(j))*(eta3(j,i+1));
			h3(j,i)   = eta3(j,i) - zb(j,i);
			qx3(j,i)  = u3(j,i)*h3(j,i);
			qy3(j,i)  = v3(j,i)*h3(j,i);
		}

		// S edge
		// --------------------------------------------
		j=ny+2;
		for(i=2; i<=(nx+1); i++)
		{
			indCurrHydr = int(bcnS(i) + fops(i));
			v3(j,i)   = fubs(i)*(fops(i)*v3(j-1,i) + fetas(i)*v3(j-1,i) + fqs(i)*currHydr(indCurrHydr)*flowDistrS(i)) + (1.0-fubs(i))*(-v3(j-1,i));
			u3(j,i)   = u3(j-1,i); // slip boundary condition
			eta3(j,i) = fubs(i)*(fops(i)*eta3(j-1,i) + fetas(i)*currHydr(indCurrHydr) + fqs(i)*eta3(j-1,i)) + (1.0-fubs(i))*(eta3(j-1,i));
			h3(j,i)   = eta3(j,i) - zb(j,i);
			qx3(j,i)  = u3(j,i)*h3(j,i);
			qy3(j,i)  = v3(j,i)*h3(j,i);
		}

		// N edge
		// --------------------------------------------
		j=1;
		for(i=2; i<=(nx+1); i++)
		{
			indCurrHydr = int(bcnN(i) + fopn(i));
			v3(j,i)   = fubn(i)*(fopn(i)*v3(j+1,i) + fetan(i)*v3(j+1,i) + fqn(i)*currHydr(indCurrHydr)*flowDistrN(i)) + (1.0-fubn(i))*(-v3(j+1,i));
			u3(j,i)   = u3(j+1,i); // slip boundary condition
			eta3(j,i) = fubn(i)*(fopn(i)*eta3(j+1,i) + fetan(i)*currHydr(indCurrHydr) + fqn(i)*eta3(j+1,i)) + (1.0-fubn(i))*(eta3(j+1,i));
			h3(j,i)   = eta3(j,i) - zb(j,i);
			qx3(j,i)  = u3(j,i)*h3(j,i);
			qy3(j,i)  = v3(j,i)*h3(j,i);
		}


		// loop through internal cells
		// ---------------------------------------------------------------------------------------
		// ---------------------------------------------------------------------------------------

		for (j=2; j<=(ny+1); j++)
		{
			for (i=2; i<=(nx+1); i++)
			{
				// implicit friction term
				// ----------------------------------------------------------------------------------
				U(2) = qx3(j,i);
				U(3) = qy3(j,i);

				if(h3(j,i) > 1e-10)
				{
					Cf  = (9.81*sqr(man(j,i)))/(pow(h3(j,i),0.33));
					Sfx = Cf*u3(j,i)*sqrt(sqr(u3(j,i)) + sqr(v3(j,i)));
					Sfy = Cf*v3(j,i)*sqrt(sqr(u3(j,i)) + sqr(v3(j,i)));
					Dx  = 1.0 + 2.0*dt*Cf*fabs(U(2))/(sqr(h3(j,i)));
					Dy  = 1.0 + 2.0*dt*Cf*fabs(U(3))/(sqr(h3(j,i)));
					Fx  = Sfx/Dx;
					Fy  = Sfy/Dy;
				}
				else
				{
					Cf  = 0.0;
					Sfx = 0.0;
					Sfy = 0.0;
					Dx  = 1.0;
					Dy  = 1.0;
					Fx  = 0.0;
					Fy  = 0.0;
				}

				// update fluxes
				Up(2) = U(2) + dt*Fx;
				Up(3) = U(3) + dt*Fy;

				// limit friction force
				if(Up(2)*U(2) < 0.0) Up(2) = U(2);
				if(Up(3)*U(3) < 0.0) Up(3) = U(3);

				// compute Sf(t+1)
				Sfxp = (Up(2) - U(2))/dt;
				Sfyp = (Up(3) - U(3))/dt;


				// reconstruct face values
				// 2nd order for wet cells with all wet neighbors, 1st order otherwise
				// ----------------------------------------------------------------------------------

				if(h3(j,i) <= 1e-10) // dry cell? -> 1st order reconstruction, zero inner depth, velocities and fluxes
				{					
					recon1dry(eta3, h3, qx3, qy3, u3, v3, zb,
						etaLe, etaRe, etaLw, etaRw, 
						etaDn, etaUn, etaDs, etaUs,
						hLe, hRe, hLw, hRw, 
						hDn, hUn, hDs, hUs,
						qxLe, qxRe, qxLw, qxRw,
						qxDn, qxUn, qxDs, qxUs, 
						qyLe, qyRe, qyLw, qyRw, 
						qyDn, qyUn, qyDs, qyUs, 
						uLe, uRe, uLw, uRw,
						uDn, uUn, uDs, uUs,
						vLe, vRe, vLw, vRw, 
						vDn, vUn, vDs, vUs,
						zbLe, zbRe, zbLw, zbRw,
						zbDn, zbUn, zbDs, zbUs,
						i,j
						);
				}
				else // wet cell
				{
					// check h in neighbors
					neighh(1)=h3(j,i+1); neighh(2)=h3(j,i-1); neighh(3)=h3(j-1,i); neighh(4)=h3(j+1,i);

					if(doubVectMin(neighh) <= 1e-10) // wet cell with at least one dry neighbor? -> 1st order reconstruction
					{
						recon1dryWet(eta3, h3, qx3, qy3, u3, v3, zb,
							etaLe, etaRe, etaLw, etaRw, 
							etaDn, etaUn, etaDs, etaUs,
							hLe, hRe, hLw, hRw, 
							hDn, hUn, hDs, hUs,
							qxLe, qxRe, qxLw, qxRw,
							qxDn, qxUn, qxDs, qxUs, 
							qyLe, qyRe, qyLw, qyRw, 
							qyDn, qyUn, qyDs, qyUs, 
							uLe, uRe, uLw, uRw,
							uDn, uUn, uDs, uUs,
							vLe, vRe, vLw, vRw, 
							vDn, vUn, vDs, vUs,
							zbLe, zbRe, zbLw, zbRw,
							zbDn, zbUn, zbDs, zbUs,
							i,j);
					}
					else // wet cell with wet neighbors -> 2nd order
					{
						// inner faces
						recon2inn(eta3, etaLe, etaRw, etaDn, etaUs,
							inVect, limslope,
							dx, dy, i, j);
						recon2inn(h3, hLe, hRw, hDn, hUs,
							inVect, limslope,
							dx, dy, i, j);
						recon2inn(qx3, qxLe, qxRw, qxDn, qxUs,
							inVect, limslope,
							dx, dy, i, j);
						recon2inn(qy3, qyLe, qyRw, qyDn, qyUs,
							inVect, limslope,
							dx, dy, i, j);

						zbLe  = etaLe - hLe;
						zbRw  = etaRw - hRw;
						zbDn  = etaDn - hDn;
						zbUs  = etaUs - hUs;

						// outer faces
						recon2out(eta3, etaRe, etaLw, etaUn, etaDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);
						recon2out(h3, hRe, hLw, hUn, hDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);
						recon2out(qx3, qxRe, qxLw, qxUn, qxDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);
						recon2out(qy3, qyRe, qyLw, qyUn, qyDs,
							inVect, limslope, 
							dx, dy, i, j, nx, ny);

						zbRe = etaRe - hRe;
						zbLw = etaLw - hLw;
						zbUn = etaUn - hUn;
						zbDs = etaDs - hDs;

						// velocities
						uLe = qxLe/hLe; uRw = qxRw/hRw; uDn = qxDn/hDn; uUs = qxUs/hUs;
						uRe = qxRe/hRe; uLw = qxLw/hLw; uUn = qxUn/hUn; uDs = qxDs/hDs;

						vLe = qyLe/hLe; vRw = qyRw/hRw; vDn = qyDn/hDn; vUs = qyUs/hUs;
						vRe = qyRe/hRe; vLw = qyLw/hLw; vUn = qyUn/hUn; vDs = qyDs/hDs;
					}
				}


				// wetting-drying (Liang 2010) & HLLC Riemann fluxes (Toro 2001)
				// -------------------------------------------------------------------------

				// east face							
				wetDry(etaLe, etaRe, etaLeWD, etaReWD,
					hLe, hRe, hLeWD, hReWD,
					qxLe, qxRe, qxLeWD, qxReWD, 
					qyLe, qyRe, qyLeWD, qyReWD, 
					uLe, uRe, vLe, vRe, 
					zbLe, zbRe, 
					zbe, dz
					);
				xFlux(Fe,
					uLe, uRe, 
					vLe, vRe,
					qxLeWD, qxReWD,
					qyLeWD, qyReWD,
					hLeWD, hReWD, 
					etaLeWD, etaReWD,
					zbe,
					UL, UR, 
					FL, FR, 
					ustar, hstar, 
					SL, SR, Sstar, 
					FstarL, FstarR
					);

				// west face
				wetDry(etaRw, etaLw, etaRwWD, etaLwWD,
					hRw, hLw, hRwWD, hLwWD,
					qxRw, qxLw, qxRwWD, qxLwWD, 
					qyRw, qyLw, qyRwWD, qyLwWD, 
					uRw, uLw, vRw, vLw, 
					zbRw, zbLw, 
					zbw, dz
					);
				xFlux(Fw,
					uLw, uRw, 
					vLw, vRw,
					qxLwWD, qxRwWD,
					qyLwWD, qyRwWD,
					hLwWD, hRwWD, 
					etaLwWD, etaRwWD,
					zbw,
					UL, UR, 
					FL, FR, 
					ustar, hstar, 
					SL, SR, Sstar, 
					FstarL, FstarR
					);

				// north face
				wetDry(etaDn, etaUn, etaDnWD, etaUnWD,
					hDn, hUn, hDnWD, hUnWD,
					qxDn, qxUn, qxDnWD, qxUnWD, 
					qyDn, qyUn, qyDnWD, qyUnWD, 
					uDn, uUn, vDn, vUn, 
					zbDn, zbUn, 
					zbn, dz
					);
				yFlux(Gn, 
					uDn, uUn, 
					vDn, vUn,
					qxDnWD, qxUnWD,
					qyDnWD, qyUnWD,
					hDnWD, hUnWD, 
					etaDnWD, etaUnWD,
					zbn,
					UD, UU, 
					GD, GU, 
					vstar, hstar, 
					SD, SU, Sstar, 
					GstarD, GstarU
					);

				// south face
				wetDry(etaUs, etaDs, etaUsWD, etaDsWD,
					hUs, hDs, hUsWD, hDsWD,
					qxUs, qxDs, qxUsWD, qxDsWD, 
					qyUs, qyDs, qyUsWD, qyDsWD, 
					uUs, uDs, vUs, vDs, 
					zbUs, zbDs, 
					zbs, dz
					);
				yFlux(Gs, 
					uDs, uUs,
					vDs, vUs,
					qxDsWD, qxUsWD,
					qyDsWD, qyUsWD,
					hDsWD, hUsWD, 
					etaDsWD, etaUsWD,
					zbs,
					UD, UU, 
					GD, GU, 
					vstar, hstar, 
					SD, SU, Sstar, 
					GstarD, GstarU
					);


				// Update cell
				// --------------------------------------------

				// get state variable values at t
				U(1) = eta(j,i);
				U(2) = qx(j,i);
				U(3) = qy(j,i);

				// compute source terms
				S(1) = 0.0;
				S(2) = -9.81*0.5*(etaRwWD + etaLeWD)*(zbe-zbw)/dx - Sfxp;
				S(3) = -9.81*0.5*(etaDnWD + etaUsWD)*(zbn-zbs)/dy - Sfyp;

				// Runge-Kutta coefficient
				K41(j,i) = -(Fe(1) - Fw(1))/dx - (Gn(1) - Gs(1))/dy + S(1);
				K42(j,i) = -(Fe(2) - Fw(2))/dx - (Gn(2) - Gs(2))/dy + S(2);
				K43(j,i) = -(Fe(3) - Fw(3))/dx - (Gn(3) - Gs(3))/dy + S(3);

				// compute final quantities
				Up(1) = U(1) + (1.0/6.0)*dt*(K11(j,i) + 2.0*K21(j,i) + 2.0*K31(j,i) + K41(j,i));
				Up(2) = U(2) + (1.0/6.0)*dt*(K12(j,i) + 2.0*K22(j,i) + 2.0*K32(j,i) + K42(j,i));
				Up(3) = U(3) + (1.0/6.0)*dt*(K13(j,i) + 2.0*K23(j,i) + 2.0*K33(j,i) + K43(j,i));


				// update state variables taking dry and negative cells into account
				updPrim(u, v, h, eta, qx, qy, zb, hcarry, i, j, Up, neighh, nhmax);


				// compute Courant number and update it if new cell > previous cell
				Cx(j-1,i-1) = dt/(dx/(fabs(u(j,i)) + sqrt(9.81*fabs(h(j,i))))); // courant number in x
				Cy(j-1,i-1)  = dt/(dy/(fabs(v(j,i)) + sqrt(9.81*fabs(h(j,i))))); // courant number in y

				/*/ stop execution if simulation becomes unstable
				if(Cx(j,i) != Cx(j,i))
				{
					s_o << "\n\n...simulation became unstable at row " << j-1 << ", column " << i-1 << "\n...aborting execution\n\n";
					exit(0);
				}
				if(Cy(j,i) != Cy(j,i))
				{
					s_o << "\n\n...simulation became unstable at row " << j-1 << ", column " << i-1 << "\n...aborting execution\n\n";
					exit(0);
				}*/
				
			}
		}

		// update time
		// ---------------------------------------------------------------------------------------
		// ---------------------------------------------------------------------------------------

		// get maximum courant number
		Cxmax = doubArrMax(Cx,nx, ny); Cymax = doubArrMax(Cy,nx, ny);				
		Cmax = max(Cxmax, Cymax);

		if(Cmax >= Ctolup) 
		{
			dt=dt*0.5;
		}
		else
		{
			if(Cmax	<= Ctollow)
			{
				dt=2.0*dt;
			}
		}

		t += dt;
	}
}

