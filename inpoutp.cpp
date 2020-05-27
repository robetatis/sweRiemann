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

