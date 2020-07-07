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

# temp change

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

