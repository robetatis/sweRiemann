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

// write output files
inline void writeScal(std::string& filename, std::string& outdir, std::string& tString, std::string& rfpre, std::string& fileext,
					  int& i, int& j, dpreal& dx, int& nx, int& ny, dpreal& xll, dpreal& yll, dpreal& noData,
					  ArrayGen(dpreal)& h, ArrayGen(dpreal)& arr
					  )
{
	filename = outdir + "/" + rfpre + "-" + tString.c_str() + "." + fileext;
	Os fstream (filename.c_str(), NEWFILE);
	fstream << "NCOLS\t"        << nx << "\n";
	fstream << "NROWS\t"        << ny << "\n";
	fstream << "XLLCORNER\t"    << xll << "\n";
	fstream << "YLLCORNER\t"    << yll << "\n";
	fstream << "CELLSIZE\t"     << dx << "\n";
	fstream << "NODATA_value\t" << noData << "\n";
	for (j=2; j<=(ny+1); ++j) 
	{
		for (i=2; i<=(nx+1); i++) 
		{
			if(h(j,i) <= 1e-10) 
				fstream << noData << "\t";
			else 
				fstream << arr(j,i) << "\t";
		}
		fstream << "\n";
	}
}

inline void writeVmag(std::string& filename, std::string& outdir, std::string& tString, std::string& rfpre, std::string& fileext,
					  int& i, int& j, dpreal& dx, int& nx, int& ny, dpreal& xll, dpreal& yll, dpreal& noData,
					  ArrayGen(dpreal)& h, ArrayGen(dpreal)& u, ArrayGen(dpreal)& v
					  )
{
	filename = outdir + "/" + rfpre + "-" + tString.c_str() + "." + fileext;
	Os fstream (filename.c_str(), NEWFILE);
	fstream << "NCOLS\t"        << nx << "\n";
	fstream << "NROWS\t"        << ny << "\n";
	fstream << "XLLCORNER\t"    << xll << "\n";
	fstream << "YLLCORNER\t"    << yll << "\n";
	fstream << "CELLSIZE\t"     << dx << "\n";
	fstream << "NODATA_value\t" << noData << "\n";
	for (j=2; j<=(ny+1); ++j) 
	{
		for (i=2; i<=(nx+1); i++) 
		{
			if(h(j,i) <= 1e-10) 
				fstream << noData << "\t";
			else 
				fstream << sqrt(sqr(u(j,i)) +  sqr(v(j,i))) << "\t";
		}
		fstream << "\n";
	}
}

inline void writeFlowDir(std::string& filename, std::string& outdir, std::string& tString, std::string& rfpre, std::string& fileext,
						 int& i, int& j, dpreal& dx, int& nx, int& ny, dpreal& xll, dpreal& yll, dpreal& noData,
						 ArrayGen(dpreal)& h, ArrayGen(dpreal)& v, ArrayGen(dpreal)& u
						 )
{
	filename = outdir + "/" + rfpre + "-" + tString.c_str() + "." + fileext;
	Os fstream (filename.c_str(), NEWFILE);
	fstream << "NCOLS\t"        << nx << "\n";
	fstream << "NROWS\t"        << ny << "\n";
	fstream << "XLLCORNER\t"    << xll << "\n";
	fstream << "YLLCORNER\t"    << yll << "\n";
	fstream << "CELLSIZE\t"     << dx << "\n";
	fstream << "NODATA_value\t" << noData << "\n";
	for (j=2; j<=(ny+1); ++j)
	{
		for (i=2; i<=(nx+1); i++)
		{
			if(h(j,i) <= 1e-10)
				fstream << noData << "\t";
			else 
				fstream << atan2(v(j,i), u(j,i)) << "\t";
		}
		fstream << "\n";
	}
}


// get runtime parameters and write them to console
void runTime(time_t& sttime, String& sttimestr)
{
	time_t endtime;
	time(&endtime);
	String endtimestr = ctime(&endtime);
	double runtime    = difftime(endtime, sttime);
	double divisor    = 0.0;
	String tunit;

	if (runtime <= 60)
	{
		divisor = 1.0; tunit = " seconds";
	}
	else
	{
		if (runtime <= 3600)
		{
			divisor = 60; tunit = " minutes";
		}
		else
		{
			if (runtime <= 86400)
			{
				divisor = 3600; tunit = " hours";
			}
			else
			{
				divisor = 86400; tunit = " days";
			}
		}
	}

	runtime = runtime/divisor;
	String runtimestr = oform("%5.3f", runtime);

	s_o << "\n...done\n";
	s_o << "\n---------------------------------------------------\nsimulation started:\t" << sttimestr ;
	s_o << "simulation ended:\t" << endtimestr ;
	s_o << "total runtime:\t\t" << runtimestr << tunit << "\n---------------------------------------------------\n";
}


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


// read control file
void readConFile(String& dummy, String& demname, String& bcname, String& etaicname, String& uicname,String& vicname, String& manname,
				 std::string& rfpre, std::string& outdir, dpreal& tstop, dpreal& dt, dpreal& outprf, dpreal& outpc, 
				 std::string& cwd, 
				 String& hout,String& etaout,String& vout,String& uout, String& vmagout, String& flowdout,
				 dpreal& Ctolup, dpreal& Ctollow)
{
	ifstream conFile("control.txt");

	s_o << "\n" << "reading control file" << "\n" << "-------------------------------------\n";

	conFile >> dummy >> tstop;
	conFile >> dummy >> dt;
	conFile >> dummy >> Ctolup;
	conFile >> dummy >> Ctollow;
	conFile >> dummy >> outprf;
	conFile >> dummy >> outpc;
	conFile >> dummy >> demname;
	conFile >> dummy >> bcname;
	conFile >> dummy >> etaicname;
	conFile >> dummy >> uicname;
	conFile >> dummy >> vicname;
	conFile >> dummy >> manname;
	conFile >> dummy >> rfpre;
	conFile >> dummy >> outdir;
	conFile >> dummy >> hout;
	conFile >> dummy >> etaout;
	conFile >> dummy >> vout;
	conFile >> dummy >> uout;
	conFile >> dummy >> vmagout;
	conFile >> dummy >> flowdout;

	s_o << "simulation time [s]                             = " << tstop << "\n";
	s_o << "time step [s]                                   = " << dt << "\n"; 
	s_o << "Max. allowed C                                  = " << Ctolup << "\n"; 
	s_o << "Min. allowed C                                  = " << Ctollow << "\n"; 
	s_o << "output frequency, results file [sec.]           = " << outprf << "\n"; 
	s_o << "output frequency, console [sec.]                = " << outpc << "\n"; 
	s_o << "DEM file                                        = " << demname << "\n"; 
	s_o << "BC file                                         = " << bcname << "\n"; 
	s_o << "water level IC file                             = " << etaicname << "\n"; 
	s_o << "u-velocity IC file                              = " << uicname << "\n"; 
	s_o << "v-velocity IC file                              = " << vicname << "\n"; 
	s_o << "manning file                                    = " << manname << "\n"; 

	// make directory for output files in current directory
	outdir = cwd + "/" + outdir;
	mkdir(outdir.c_str());

	conFile.close();
}

// read header of DEM file
void readDEMhead(String& demname, int& nx, int& ny, dpreal& xll, dpreal& yll, dpreal& dx, dpreal& dy, dpreal& noData, String& dummy,
				 int& i, int& j)
{
	ifstream demFile(demname.c_str());

	demFile >> dummy >> nx; 
	demFile >> dummy >> ny; 
	demFile >> dummy >> xll;
	demFile >> dummy >> yll;
	demFile >> dummy >> dx;
	dy = dx;
	demFile >> dummy >> noData;

	// write this data out to console
	s_o << "\n" << "reading DEM" << "\n" << "-------------------------------------\n";
	s_o << "ncol (x-coord.)                  = " << nx << "\n";
	s_o << "nrow (y-coord.)                  = " << ny << "\n";
	s_o << "xll (x-coord. lower left corner) = " << xll << "\n";
	s_o << "yll (y-coord. lower left corner) = " << yll << "\n";
	s_o << "cell size (dx=dy)                = " << dx << "\n";
	s_o << "NoData value                     = " << noData << "\n";

	demFile.close();
}

// read DEM elevations
void readDEMz(String& demname, int& nx, int& ny, int& i, int& j, ArrayGen(dpreal)& zb)
{
	ifstream demFile(demname.c_str());

	// loop past header
	String dummy;
	for(int k=1; k<=12; k++) demFile >> dummy;

	for(j=2; j<=(ny+1); j++)
		for(i=2; i<=(nx+1); i++)
			demFile >> zb(j,i);

	// populate ghost cells with zb of adjacent internal cells
	i=1;
	for(j=2; j<=(ny+1); j++) zb(j,i) = zb(j,i+1);
	i=nx+2;
	for(j=2; j<=(ny+1); j++) zb(j,i) = zb(j,i-1);
	j=1;
	for(i=2; i<=(nx+1); i++) zb(j,i) =  zb(j+1,i);
	j=ny+2;
	for(i=2; i<=(nx+1); i++) zb(j,i) =  zb(j-1,i);
}

// read Manning raster
void readMan(String& manname, int& nx, int& ny,  String& dummy, int& i, int& j, ArrayGen(dpreal)& man)
{
	s_o << "\n" << "reading Manning file" << "\n" << "-------------------------------------\n";

	ifstream manFile(manname.c_str());

	// read file header
	for (i=1; i<=6; i++)
		manFile >> dummy >> dummy;

	// read Manning values and populate array man
	for(j=2; j<=(ny+1); j++)
		for(i=2; i<=(nx+1); i++)
			manFile >> man(j,i);
}


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


// set boundary conditions according to BC file; default mode is closed boundary
void setBC(dpreal& dx, dpreal& dy, dpreal& xll, dpreal& yll, int& nx, int& ny, dpreal& tstop,
		   ArrayGen(dpreal)& u, ArrayGen(dpreal)& v, ArrayGen(dpreal)& h, ArrayGen(dpreal)& eta, ArrayGen(dpreal)& zb, ArrayGen(dpreal)& man,
		   int& i, int& j,
		   String& bcname, String& dummy, 
		   Vec(dpreal)& fube, Vec(dpreal)& fubw, Vec(dpreal)& fubs, Vec(dpreal)& fubn,
		   Vec(dpreal)& fope, Vec(dpreal)& fopw, Vec(dpreal)& fops, Vec(dpreal)& fopn,
		   Vec(dpreal)& fetae, Vec(dpreal)& fetaw, Vec(dpreal)& fetas, Vec(dpreal)& fetan,
		   Vec(dpreal)& fqe, Vec(dpreal)& fqw, Vec(dpreal)& fqs, Vec(dpreal)& fqn,
		   Vec(dpreal)& flowDistrE, Vec(dpreal)& flowDistrW, Vec(dpreal)& flowDistrS, Vec(dpreal)& flowDistrN,
		   ArrayGen(dpreal)& hydrgArr, int& nrowHydrgArr, int& ncolHydrgArr,
		   VecSimple(int)& hydrgNoSteps,
		   Vec(dpreal)& bcnE, Vec(dpreal)& bcnW, Vec(dpreal)& bcnS, Vec(dpreal)& bcnN
		   )
{
	s_o << "\nsetting boundary conditions\n-------------------------------------\n";


	// declare variables to hold values in BC file
	// --------------------------------------------------------------------------------------------------------------
	String btype, bvaluefile;
	char   bside;
	dpreal bstart=0.0, bend=0.0;
	int	   bstartc=0, bendc=0, bcelln=0, k=0, lc=0, hydrgNo=0;

	// hydrograph counter
	int	   hydrgCount=0; 


	// this section reads all hydrograph files and stores them in array hydrgArr
	// --------------------------------------------------------------------------------------------------------------
	ifstream bcFile(bcname.c_str());

	for(k=1; k<=5; k++) bcFile >> dummy;
	while(bcFile >> bside >> bstart >> bend >> btype >> bvaluefile) 
		if(btype == "eta" || btype == "Q")
			hydrgNo++;
	bcFile.close();

	// if hydrographs exist in BC file (either eta or Q hydrographs), loop through it again, go to each hydrograph file, get its no. of lines
	// and re-size the hydrographs array
	if(hydrgNo > 0)
	{
		hydrgNoSteps.redim(hydrgNo);
		bcFile.open(bcname.c_str());
		for(k=1; k<=5; k++) bcFile >> dummy;
		k=0;
		while(bcFile >> bside >> bstart >> bend >> btype >> bvaluefile)
		{	
			if(btype == "eta" || btype == "Q")
			{
				k++;
				ifstream hydrgFile(bvaluefile.c_str());
				int currHydrgSteps=0;
				while(hydrgFile >> dummy >> dummy) currHydrgSteps++; // count no. of lines of current hydrograph file (incl. column titles)
				hydrgNoSteps(k)=currHydrgSteps-1;
				hydrgFile.close();
			}
		}
		// re-size hydrograph array
		nrowHydrgArr=intVectMax(hydrgNoSteps);
		ncolHydrgArr=2*hydrgNo;
		hydrgArr.redim(nrowHydrgArr, ncolHydrgArr);
		hydrgArr.fill(0.0);
	}
	bcFile.close();

	// fill the hydrographs array
	bcFile.open(bcname.c_str());
	if(hydrgNo > 0)
	{
		for(k=1; k<=5; k++) bcFile >> dummy;
		k=0;
		while(bcFile >> bside >> bstart >> bend >> btype >> bvaluefile)
		{	
			if(btype == "eta" || btype == "Q")
			{
				k++;
				ifstream hydrgFile(bvaluefile.c_str());
				hydrgFile >> dummy >> dummy;
				for(int ll=1; ll<=hydrgNoSteps(k); ll++)
				{
					hydrgFile >> hydrgArr(ll, 2*k-1) >> hydrgArr(ll, 2*k);

					// stop execution if one of the hydrographs is shorter than simulation time
					if(ll == hydrgNoSteps(k)) 
					{
						if(hydrgArr(ll, 2*k-1) < tstop)
						{
							s_o << "\n\nERROR: hydrograph no. " << 2*k-1 << " is shorter than simulation time, aborting execution\n";
							exit(1);
						}
					}

					// stop execution if time axis in hydrograph is not increasing
					if(ll>1)
					{
						if(hydrgArr(ll, 2*k-1) <= hydrgArr(ll-1, 2*k-1))
						{
							s_o << "\n\nERROR: time axis in hydrograph no. " << 2*k-1 << " is not increasing, aborting execution\n";
							exit(1);
						}
					}
				}
				hydrgFile.close();
			}	
		}
		bcFile.close();
	}


	// allocate boundary IDs, boundary type flags and flow distributor (only for Q BCs) in ghost cells
	// --------------------------------------------------------------------------------------------------------------
	try
	{
		ifstream bcFile(bcname.c_str());

		for(k=1; k<=5; k++) bcFile >> dummy;
		while(bcFile >> bside >> bstart >> bend >> btype >> bvaluefile)
		{
			lc++;

			// southern boundary
			// ----------------------------------
			if(bside == 'S')
			{
				j=ny+2;

				if(bstart >= bend)
				{
					s_o << "\n\nERROR: BCs must be entered in increasing x direction, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
					throw 1;
				}

				bstartc = int((bstart - xll)/dx + 2);
				bendc   = int((bend - xll)/dx +1);

				// make sure user-defined boundaries never extend beyond the DEM's edges, if they do, they are automatically shifted 
				if(bstartc <= 2) bstartc = 2;
				if(bendc >= nx+1) bendc = nx+1;

				bcelln = bendc - bstartc + 1;

				// water level boundary condition
				if(btype == "eta")
				{
					hydrgCount++;
					s_o << "water level [m] BC on S edge, from x = " << bstart << " m to x = " << bend << " m\n";

					for(i=bstartc; i<=bendc; i++)
					{
						if(fubs(i) != 0.0)
						{
							s_o << "\n\nERROR: overlapping boundary conditions, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
							throw 1;
						}
						else
						{							
							bcnS(i)  = hydrgCount;
							fubs(i)  = 1.0;
							fetas(i) = 1.0;
						}
					}
				}

				// flow boundary condition
				else
				{
					if(btype == "Q")
					{
						hydrgCount++;
						s_o << "flow [m3/s] BC on S edge, from x = " << bstart << " m to x = " << bend << " m\n";

						// sum water depths, used for distributing imposed flow along boundary 
						// proportionally to each cell's depth
						dpreal sumh = 0.0;								 
						for(i=bstartc; i<=bendc; i++) sumh += h(j-1,i);  
						for(i=bstartc; i<=bendc; i++)
						{
							if(h(j-1,i) <= 1e-10)
							{
								s_o << "\n\nERROR: trying to impose flow BC on a dry boundary cell, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
								throw 1;
							}
							if(fubs(i) != 0.0)
							{
								s_o << "\n\nERROR: overlapping boundary conditions, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
								throw 1;
							}
							bcnS(i)  = hydrgCount;
							fubs(i) = 1.0; 
							fqs(i) = 1.0;
							flowDistrS(i) = h(j-1,i)/sumh;
						}
					}

					// open boundary condition
					else
					{
						if(btype == "open")
						{
							s_o << "open BC on S edge, from x = " << bstart << " m to x = " << bend << " m\n";
							for(i=bstartc; i<=bendc; i++)
							{
								if(fubs(i) != 0.0)
								{
									s_o << "\n\nERROR: overlapping boundary conditions, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
									throw 1;
								}
								else
								{ 
									fubs(i) = 1.0; 
									fops(i) = 1.0;
								}
							}
						}
						else
						{
							s_o << "\n\nERROR: wrong boundary identifier, see line no. " << lc+1 << " in BC file, must be 'eta', 'Q' or 'open', aborting execution" << "\n"; 
							throw 1;
						}
					}
				}
			}
			else
			{

				// eastern boundary
				// ----------------------------------
				if(bside == 'E')
				{
					i=nx+2;

					if(bstart >= bend)
					{				
						s_o << "\n\nERROR: BCs must be entered in increasing y direction, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
						throw 1;
					}

					bstartc = int((bstart - yll)/dy + 2);
					bendc   = int((bend - yll)/dy + 1);

					// make sure user-defined boundaries never extend beyond the DEM's edges, if they do, they are automatically shifted 
					if(bstartc <= 2){ bstartc = ny+1;} else {bstartc = ny - bstartc + 3;}
					if(bendc >= ny+1){ bendc = 2;} else {bendc = ny - bendc + 3;}

					bcelln = bstartc - bendc + 1;

					// water level boundary condition
					if(btype == "eta")
					{
						hydrgCount++;
						s_o << "water level [m] BC on E edge, from y = " << bstart << " m to y = " << bend << " m\n";
						for(j=bstartc; j>=bendc; --j)
						{
							if(fube(j) != 0.0)
							{
								s_o << "ERROR: overlapping boundary conditions, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
								throw 1;
							}
							else
							{
								bcnE(j)  = hydrgCount;
								fube(j)  = 1.0; 
								fetae(j) = 1.0;
							}
						}
					}

					// flow boundary condition
					else
					{
						if(btype == "Q")
						{
							hydrgCount++;
							s_o << "flow [m3/s] BC on E edge, from y = " << bstart << " m to y = " << bend << " m\n";

							// sum water depths, used for distributing imposed flow along boundary 
							// proportionally to each cell's depth
							dpreal sumh = 0.0;								 
							for(j=bstartc; j>=bendc; --j) sumh += h(j,i-1);  
							for(j=bstartc; j>=bendc; --j)
							{
								if(h(j,i-1)<=1e-10)
								{
									s_o << "\n\nERROR: trying to impose flow BC on a dry boundary cell, see line no. " << lc+1 << " in BC file, aborting execution\n";									       
									throw 1;
								}
								if(fube(j) != 0.0)
								{
									s_o << "\n\nERROR: overlapping boundary conditions, see line no. " << lc+1 << " in BC file, aborting execution\n";
									throw 1;
								}
								bcnE(j)  = hydrgCount;
								fube(j) = 1.0; 
								fqe(j) = 1.0;
								flowDistrE(j) = h(j,i-1)/sumh;
							}
						}

						// open boundary condition
						else
						{
							if(btype == "open")
							{
								s_o << "open BC on E edge, from y = " << bstart << " m to y = " << bend << " m\n";
								for(j=bstartc; j>=bendc; --j)
								{
									if(fube(j) != 0.0){s_o << "\n\nERROR: overlapping boundary conditions, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
									throw 1;
									}
									else
									{ 
										fube(j) = 1.0; 
										fope(j) = 1.0;
									}
								}
							}
							else
							{
								s_o << "\n\nERROR: wrong boundary identifier, see line no. " << lc+1 << " in BC file, must be 'eta', 'Q' or 'open', aborting execution" << "\n"; 
								throw 1;
							}
						}
					}
				}
				else
				{

					// western boundary
					// --------------------------
					if(bside == 'W')
					{
						i=1;

						if(bstart >= bend)
						{
							s_o << "\n\nERROR: BCs must be entered in increasing y direction, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
							throw 1;
						}

						bstartc = int((bstart - yll)/dy + 2);
						bendc   = int((bend - yll)/dy + 1);

						// make sure user-defined boundaries never extend beyond the DEM's edges, if they do, they are automatically shifted 
						if(bstartc <= 2){ bstartc = ny+1;} else {bstartc = ny - bstartc + 3;}
						if(bendc >= ny+1){ bendc = 2;} else {bendc = ny - bendc + 3;}

						bcelln = bstartc - bendc + 1;

						// water level boundary condition
						if(btype == "eta")
						{
							hydrgCount++;
							s_o << "water level [m] BC on W edge, from y = " << bstart << " m to y = " << bend << " m\n";
							for(j=bstartc; j>=bendc; --j)
							{
								if(fubw(j) != 0.0)
								{
									s_o << "\n\nERROR: overlapping boundary conditions, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
									throw 1;
								}
								else
								{
									bcnW(j)  = hydrgCount;
									fubw(j)  = 1.0; 
									fetaw(j) = 1.0;
								}
							}
						}

						// flow boundary condition
						else
						{
							if(btype == "Q")
							{
								hydrgCount++;
								s_o << "flow [m3/s] BC on W edge, from y = " << bstart << " m to y = " << bend << " m\n";

								// sum water depths, used for distributing imposed flow along boundary 
								// proportionally to each cell's depth
								dpreal sumh = 0.0;								 
								for(j=bstartc; j>=bendc; --j) sumh += h(j,i+1);  
								for(j=bstartc; j>=bendc; --j)
								{
									if(h(j,i+1) <= 1e-10)
									{
										s_o << "\n\nERROR: trying to impose flow BC on a dry boundary cell, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
										throw 1;
									}
									if(fubw(j) != 0.0) 
									{
										s_o << "\n\nERROR: overlapping boundary conditions, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
										throw 1;
									}
									bcnW(j)  = hydrgCount;
									fubw(j)  = 1.0; 
									fqw(j)   = 1.0;
									flowDistrW(j) = h(j,i+1)/sumh;
								}
							}

							// open boundary condition
							else
							{
								if(btype == "open")
								{
									s_o << "open BC on W edge, from y = " << bstart << " m to y = " << bend << " m\n";
									for(j=bstartc; j>=bendc; --j)
									{
										if(fubw(j) != 0.0)
										{
											s_o << "\n\nERROR: overlapping boundary conditions, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
											throw 1;
										}
										else
										{ 
											fubw(j) = 1.0; fopw(j) = 1.0;
										}
									}
								}
								else{s_o << "\n\nERROR: wrong boundary identifier on boundary no. " << lc << " in BC file, must be 'eta', 'Q' or 'open', aborting execution" << "\n";
								throw 1;
								}
							}
						}
					}
					else
					{

						// northern boundary
						// ----------------------------------
						if(bside == 'N')
						{
							j=1;

							if(bstart >= bend)
							{
								s_o << "\n\nERROR: BCs must be entered in increasing x direction, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
								throw 1;
							}

							bstartc = int((bstart - xll)/dx + 2);
							bendc   = int((bend - xll)/dx +1);

							// make sure user-defined boundaries never extend beyond the DEM's edges, if they do, they are automatically shifted 
							if(bstartc <= 2) bstartc = 2;
							if(bendc >= nx+1) bendc = nx+1;

							bcelln = bendc - bstartc + 1;

							// water level boundary condition
							if(btype == "eta")
							{
								hydrgCount++;
								s_o << "water level [m] BC on N edge, from x = " << bstart << " m to x = " << bend << " m\n";
								for(i=bstartc; i<=bendc; i++)
								{
									if(fubn(i) != 0.0)
									{
										s_o << "\n\nERROR: overlapping boundary conditions, see line no. " << lc+1 << " in BC file, aborting execution\n";
										throw 1;
									}
									else
									{
										bcnN(i)  = hydrgCount;
										fubn(i)  = 1.0; 
										fetan(i) = 1.0;
									}
								}
							}

							// flow boundary condition
							else
							{
								if(btype == "Q")
								{
									hydrgCount++;
									s_o << "flow [m3/s] BC on N edge, from x = " << bstart << " m to x = " << bend << " m\n";

									// sum water depths, used for distributing imposed flow along boundary 
									// proportionally to each cell's depth
									dpreal sumh = 0.0; 
									for(i=bstartc; i<=bendc; i++) sumh += h(j+1,i); 
									for(i=bstartc; i<=bendc; i++)
									{
										if(h(j+1,i)<=1e-10)
										{
											s_o << "\n\nERROR: trying to impose flow BC on a dry boundary cell, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
											throw 1;
										}
										if(fubn(i) != 0.0)
										{
											s_o << "\n\nERROR: overlapping boundary conditions, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
											throw 1;
										}
										bcnN(i)  = hydrgCount;
										fubn(i) = 1.0; 
										fqn(i) = 1.0;
										flowDistrN(i) = h(j+1,i)/sumh;
									}
								}

								// open boundary condition
								else
								{
									if(btype == "open")
									{
										s_o << "open BC on N edge, from x = " << bstart << " m to x = " << bend << " m\n";
										for(i=bstartc; i<=bendc; i++)
										{
											if(fubn(i) != 0.0)
											{
												s_o << "\n\nERROR: overlapping boundary conditions, see line no. " << lc+1 << " in BC file, aborting execution\n"; 
												throw 1;
											}
											else
											{ 
												fubn(i) = 1.0; fopn(i) = 1.0;
											}
										}
									}
									else
									{
										s_o << "\n\nERROR: wrong boundary identifier, see line no. " << lc+1 << " in BC file, must be 'eta', 'Q' or 'open', aborting execution" << "\n"; 
										throw 1;
									}
								}
							}
						}
						else 
						{
							s_o << "\n\nERROR: wrong edge identifier on line " << lc+1 << " of BC file, aborting execution\n"; 
							throw 1;
						}
					}
				}
			}
		}
	} 
	catch(...) {exit(1);}
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


// build time labels for output files
inline void timeLabel(dpreal& tstop, dpreal& dt, 
					  dpreal& outprf, dpreal& outpc, 
					  int& wwfile, int& wwcons, 
					  int& precfile, int& preccons)
{
	if(tstop > 1)
	{
		if(outprf < 1)
		{
			precfile = int(abs(log10(outprf)));
			wwfile   = int(abs(log10(tstop)))+2 + precfile;
		}
		else
		{
			wwfile   = int(abs(log10(tstop)))+1;
			precfile = 0;
		}
	}
	else
	{
		wwfile   = int(abs(log10(outprf)))+2;
		precfile = int(abs(log10(outprf)));
	}

	if(tstop > 1)
	{
		if(outpc < 1)
		{
			preccons = int(abs(log10(outpc)));
			wwcons   = int(abs(log10(tstop)))+2 + preccons;
		}
		else
		{
			wwcons   = int(abs(log10(tstop)))+1;
			preccons = 0;
		}
	}
	else
	{
		wwcons   = int(abs(log10(outpc)))+2;
		preccons = int(abs(log10(outpc)));
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
