#include "EvolutionData.hxx"

#include "LogFile.hxx"
#include "StringLine.hxx"

#include <TGraph.h>
#include <TString.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <algorithm>

	//________________________________________________________________________
	//
	//		EvolutionData
	//
	//
	//
	//
	//________________________________________________________________________



struct my_tolower
{
	char operator()(char c) const
	{
		return std::tolower(static_cast<unsigned char>(c));
	}
};

	//To Lower Case, convert any string in lower case
string tlc(string data)
{
	transform(data.begin(), data.end(), data.begin(), my_tolower());
	return data;
}

	//________________________________________________________________________
	//________________________________________________________________________
	//________________________________________________________________________
EvolutionData operator*(EvolutionData const& evol, double F)
{
	
	EvolutionData evoltmp;
	
	map<ZAI ,TGraph* > EvolutionData = evol.GetEvolutionData();
	map<ZAI ,TGraph* >::iterator it;
	for(it = EvolutionData.begin(); it != EvolutionData.end(); it++)
	{
		double X[(*it).second->GetN()];
		double Y[(*it).second->GetN()];
		
		for(int i = 0; i < (*it).second->GetN(); i++)
		{
			double y;
			(*it).second->GetPoint( i, X[i], y );
			Y[i] = y*F;
		}
		evoltmp.NucleiInsert( pair<ZAI, TGraph*> ( (*it).first,new TGraph((*it).second->GetN(), X, Y) ) );
		
	}
	evoltmp.SetPower(evol.GetPower()*F);
	return evoltmp;
	
}

	//________________________________________________________________________
EvolutionData operator*(double F, EvolutionData const& evol)
{
	
	return evol*F;
	
}
	//________________________________________________________________________
EvolutionData operator/(EvolutionData const& evol, double F)
{
	
	return evol*(1./F);
	
}

	//________________________________________________________________________
	//________________________________________________________________________


ClassImp(EvolutionData)


EvolutionData::EvolutionData()
{
}

	//________________________________________________________________________
EvolutionData::EvolutionData(LogFile* Log)
{
	
	SetLog(Log);
	fIsCrossSection = false;
	
}

	//________________________________________________________________________
EvolutionData::EvolutionData(LogFile* Log, string DB_file, bool oldread, ZAI zai)
{
	
	SetLog(Log);
	fIsCrossSection = false;
	fDB_file = DB_file;
	
	if(zai != ZAI(0,0,0))
		AddAsStable(zai);
	else
		ReadDB( fDB_file, oldread);		// Read Evolution Produc DB file name
	
	
	
	
}

	//________________________________________________________________________
EvolutionData::~EvolutionData()
{
	
	
}


bool EvolutionData::NucleiInsert(pair<ZAI, TGraph*> zaitoinsert)
{
	
	pair<map<ZAI, TGraph*>::iterator, bool> IResult;
	IResult = fEvolutionData.insert( zaitoinsert);
	return IResult.second;
	
}

bool EvolutionData::FissionXSInsert(pair<ZAI, TGraph*> zaitoinsert)
{
	
	pair<map<ZAI, TGraph*>::iterator, bool> IResult;
	IResult = fFissionXS.insert( zaitoinsert);
	return IResult.second;
	
}

bool EvolutionData::CaptureXSInsert(pair<ZAI, TGraph*> zaitoinsert)
{
	
	pair<map<ZAI, TGraph*>::iterator, bool> IResult;
	IResult = fCaptureXS.insert( zaitoinsert);
	return IResult.second;
	
}

bool EvolutionData::n2nXSInsert(pair<ZAI, TGraph*> zaitoinsert)
{
	
	pair<map<ZAI, TGraph*>::iterator, bool> IResult;
	IResult = fn2nXS.insert( zaitoinsert);
	return IResult.second;
	
}

	//________________________________________________________________________
void EvolutionData::AddAsStable(ZAI zai)
{
	
	double time[2] = {0, (500*365.25*3600*24)};
	double quantity[2] = {1., 1.};
	
	fEvolutionData.insert(pair<ZAI ,TGraph* >(zai, new TGraph(2, time, quantity) ) );
	
}

	//________________________________________________________________________
Double_t EvolutionData::Interpolate(double t, TGraph& EvolutionGraph)
{
	
	TString fOption;
	return (double)EvolutionGraph.Eval(t,0x0,fOption);
	
}

	//________________________________________________________________________
TGraph*	EvolutionData::GetEvolutionTGraph(const ZAI& zai)
{
	
	map<ZAI ,TGraph *>::iterator it = GetEvolutionData().find(zai) ;
	
	if ( it != GetEvolutionData().end() )
		return it->second;
	else
		return new TGraph();
	
	
}

	//________________________________________________________________________
IsotopicVector	EvolutionData::GetIsotopicVectorAt(double t)
{
	
	IsotopicVector IsotopicVectorTmp;
	map<ZAI ,TGraph* >::iterator it;
	for( it = fEvolutionData.begin(); it != fEvolutionData.end(); it++ )
	{
		IsotopicVectorTmp.Add( (*it).first, Interpolate(t, *((*it).second)) );
	}
	
	
	return IsotopicVectorTmp;
}

	//________________________________________________________________________
double	EvolutionData::GetGetXSForAt(double t, ZAI zai, int ReactionId)
{
	
	map<ZAI ,TGraph* > XSEvol;
	switch(ReactionId)
	{
		case 1: XSEvol = GetFissionXS();
			break;
		case 2: XSEvol = GetCaptureXS();
			break;
		case 3: XSEvol = Getn2nXS();
			break;
		default:cout << "!!Error!! !!!EvolutionData!!! \n Wrong ReactionId !!" << endl;
			GetLog()->fLog << "!!Error!! !!!EvolutionData!!! \n Wrong ReactionId !!" << endl;
			exit(1);
	}
	
	map<ZAI ,TGraph* >::iterator it = XSEvol.find(zai);
	
	
	if (it == XSEvol.end())
		return 0.;
	else
		return Interpolate(t, *((*it).second) );
	
}











	//________________________________________________________________________

	//________________________________________________________________________




void EvolutionData::ReadDB(string DBfile, bool oldread)
{
	
	if(oldread)
	{
		OldReadDB(DBfile);
		return;
	}

	ifstream DecayDB(DBfile.c_str());	// Open the File
	if(!DecayDB)				//check if file is correctly open
	{
		cout << "!!Warning!! !!!EvolutionData!!! \n Can't open \"" << DBfile << "\"\n" << endl;
		GetLog()->fLog << "!!Warning!! !!!EvolutionData!!! \n Can't open \"" << DBfile << "\"\n" << endl;
	}
	vector<double> vTime;
	
	string line;
	int start = 0;
	
	getline(DecayDB, line);
	if( tlc(StringLine::NextWord(line, start, ' ')) != "time")
	{
		cout << "!!Bad Trouble!! !!!EvolutionData!!! Bad Database file : " <<  DBfile << endl;
		cout << "!!Bad Trouble!! !!!EvolutionData!!! The first Line MUST be the time line !!!" << endl;
		GetLog()->fLog << "!!Bad Trouble!! !!!EvolutionData!!! Bad Database file : " <<  DBfile << endl;
		GetLog()->fLog << "!!Bad Trouble!! !!!EvolutionData!!! The first Line MUST be the time line !!!" << endl;
		exit (1);
	}
	
	while(start < (int)line.size())
		vTime.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
	
	fFinalTime = vTime.back();
	double Time[vTime.size()];
	for(int i=0; i < (int)vTime.size();i++)
		Time[i] = vTime[i];
	const int NTimeStep = sizeof(Time)/sizeof(double);

	
	enum { Keff=1, Flux, Inv, XSFis, XSCap, XSn2n };
	
	map<string, int> keyword_map;
	keyword_map["keff"] = Keff;
	keyword_map["flux"] = Flux;
	keyword_map["xsfis"] = XSFis;
	keyword_map["xscap"] = XSCap;
	keyword_map["xsn2n"] = XSn2n;
	keyword_map["inv"] = Inv;

	getline(DecayDB, line);

	do
	{
		start = 0;
		switch (keyword_map[tlc(StringLine::NextWord(line, start, ' '))])
		{
			case Keff:
				ReadKeff(line, Time, NTimeStep);
				break;
				
			case Flux:
				ReadFlux(line, Time, NTimeStep);
				break;

			case Inv:
				ReadInv(line, Time, NTimeStep);
				break;
				
			case XSFis:
				ReadXSFis(line, Time, NTimeStep);
				break;
				
			case XSCap:
				ReadXSCap(line, Time, NTimeStep);
				break;
				
			case XSn2n:
				ReadXSn2n(line, Time, NTimeStep);
				break;
				
				
			default:
				break;
		}
		
		getline(DecayDB, line);

		
	}while ( !DecayDB.eof() );
		
		
		
}

void EvolutionData::ReadKeff(string line, double* time, int NTimeStep)
{
	

	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "keff" )	// Check the keyword
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"keff\" not found !" << endl;
		GetLog()->fLog << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"keff\" not found !" << endl;
		exit(1);
	}
	
	
	double Keff[NTimeStep];

	int i = 0;
	while(start < (int)line.size() && i < NTimeStep )		// Read the Data
	{
		Keff[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
		i++;
	}
	
	fFlux = new TGraph(NTimeStep, time, Keff);			// Add the TGraph
	
	
	
	
}

void EvolutionData::ReadFlux(string line, double* time, int NTimeStep)
{
	
	
	int start = 0;
	
	if( tlc(StringLine::NextWord(line, start, ' ')) != "flux" )	// Check the keyword
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"flux\" not found !" << endl;
		GetLog()->fLog << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"flux\" not found !" << endl;
		exit(1);
	}
	
	double Flux[NTimeStep];
	
	int i = 0;
	while(start < (int)line.size() && i < NTimeStep )		// Read the Data
	{
		Flux[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
		i++;
	}
	
	
	fFlux = new TGraph(NTimeStep, time, Flux);			// Add the TGraph
	
	ReadInfo();							// Read the .info associeted

	
}


void	EvolutionData::ReadInv(string line, double* time, int NTimeStep)
{
	
	
	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "inv" )	// Check the keyword
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"inv\" not found !" << endl;
		GetLog()->fLog << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"inv\" not found !" << endl;
		exit(1);
	}
		// Read the Z A I
	int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
	
	if(A!=0 && Z!=0)
	{
		double Inv[NTimeStep];
		
		int i = 0;
		while(start < (int)line.size() && i < NTimeStep )	// Read the Data
		{
			Inv[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
			i++;
		}
			// Add the TGraph
		fEvolutionData.insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph(NTimeStep, time, Inv) ) );
	}
	
	
}


void	EvolutionData::ReadXSFis(string line, double* time, int NTimeStep)
{
	

	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "xsfis" )	// Check the keyword
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"xsfis\" not found !" << endl;
		GetLog()->fLog << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"xsfis\" not found !" << endl;
		exit(1);
	}
		// Read the Z A I
	int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
	
	if(A!=0 && Z!=0)
	{
		double XSFis[NTimeStep];
		
		int i = 0;
		while(start < (int)line.size() && i < NTimeStep )	// Read the Data
		{
			XSFis[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
			i++;
		}

			// Add the TGraph
		fFissionXS.insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph(NTimeStep, time, XSFis) ) );
	}
	
	
}

void	EvolutionData::ReadXSCap(string line, double* time, int NTimeStep)
{
	

	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "xscap" )	// Check the keyword
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"xscap\" not found !" << endl;
		GetLog()->fLog << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"xscap\" not found !" << endl;
		exit(1);
	}
		// Read the Z A I
	int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
	
	if(A!=0 && Z!=0)
	{
		double XSCap[NTimeStep];
		
		int i = 0;
		while(start < (int)line.size() && i < NTimeStep )	// Read the Data
		{
			XSCap[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
			i++;
		}

			// Add the TGraph
		fCaptureXS.insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph(NTimeStep, time, XSCap) ) );
	}
	
	
}

void	EvolutionData::ReadXSn2n(string line, double* time, int NTimeStep)
{
	

	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "xsn2n" )	// Check the keyword
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"xsn2n\" not found !" << endl;
		GetLog()->fLog << "!!ERROR!! !!!EvolutionData!!! \n Bad keyword : \"xsn2n\" not found !" << endl;
		exit(1);
	}
		// Read the Z A I
	int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
	
	if(A!=0 && Z!=0)
	{
		double XSn2n[NTimeStep];
		
		int i = 0;
		while(start < (int)line.size() && i < NTimeStep )	// Read the Data
		{
			XSn2n[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
			i++;
		}
			// Add the TGraph
		fn2nXS.insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph(NTimeStep, time, XSn2n) ) );
	}
	
	
}



void EvolutionData::ReadInfo()
{
	string InfoDBFile  = fDB_file.erase(fDB_file.size()-3,fDB_file.size());
	InfoDBFile += "Info";
	ifstream InfoDB(InfoDBFile.c_str());				// Open the File
	if(!InfoDB)
	{
		cout << "!!ERROR!! !!!EvolutionData!!! \n Can't open \"" << InfoDBFile << "\"\n" << endl;
		GetLog()->fLog << "!!ERROR!! !!!EvolutionData!!! \n Can't open \"" << InfoDBFile << "\"\n" << endl;
		exit(1);
	}
	
	int start = 0;
	string line;
	getline(InfoDB, line);
	if ( tlc(StringLine::NextWord(line, start, ' ')) == "reactor")
		fReactorType =  StringLine::NextWord(line, start, ' ');
	start = 0;
	getline(InfoDB, line);
	if (tlc(StringLine::NextWord(line, start, ' ')) == "fueltype")
		fFuelType =  StringLine::NextWord(line, start, ' ');
	start = 0;
	getline(InfoDB, line);
	if ( tlc(StringLine::NextWord(line, start, ' ')) == "cycletime")
		fCycleTime =  atof(StringLine::NextWord(line, start, ' ').c_str());;
	getline(InfoDB, line); // Assembly HM Mass
	start = 0;
	getline(InfoDB, line);
	if ( tlc(StringLine::NextWord(line, start, ' ')) == "constantpower")
		fPower =  atof(StringLine::NextWord(line, start, ' ').c_str());

	
	getline(InfoDB, line); // cutoff
	getline(InfoDB, line); // NUmber of Nuclei
	start = 0;
	getline(InfoDB, line);
	if ( tlc(StringLine::NextWord(line, start, ' ')) == "normalizationfactor")
	{
		double NormFactor = atof(StringLine::NextWord(line, start, ' ').c_str());
		fPower = fPower * NormFactor;
	}
	start = 0;
	getline(InfoDB, line);
	if ( tlc(StringLine::NextWord(line, start, ' ')) == "finalheavymetalmass")
		fHMMass =  atof(StringLine::NextWord(line, start, ' ').c_str());
}

void EvolutionData::OldReadDB(string DBfile)
{
	
	{
		ifstream DecayDB(DBfile.c_str());							// Open the File
		if(!DecayDB)
		{
			cout << "!!Warning!! !!!EvolutionData!!! \n Can't open \"" << DBfile << "\"\n" << endl;
			GetLog()->fLog << "!!Warning!! !!!EvolutionData!!! \n Can't open \"" << DBfile << "\"\n" << endl;
		}
		vector<double> vTime;
		vector<double> vTimeErr;
		
		string line;
		int start = 0;
		
		getline(DecayDB, line);
		if( StringLine::NextWord(line, start, ' ') != "time")
		{
			cout << "!!Bad Trouble!! !!!EvolutionData!!! Bad Database file : " <<  DBfile << endl;
			GetLog()->fLog << "!!Bad Trouble!! !!!EvolutionData!!! Bad Database file : " <<  DBfile << endl;
			exit (1);
		}
		
		while(start < (int)line.size())
			vTime.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
		
		fFinalTime = vTime.back();
		double Time[vTime.size()];
		for(int i=0; i < (int)vTime.size();i++)
			Time[i] = vTime[i];
		vector<double> vFlux;
		start = 0;
		getline(DecayDB, line);
		string tmp = StringLine::NextWord(line, start, ' ');
		if ( tmp == "keff"  || tmp == "Keff" )
		{
			vector<double> vKeff;
			while(start < (int)line.size())
				vKeff.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
			
			double Keff[vKeff.size()];
			for(int i=0; i < (int)vKeff.size();i++)
				Keff[i] = vKeff[i];
			
			fKeff = new TGraph(vTime.size(), Time, Keff);
			
			start = 0;
			getline(DecayDB, line);
			if (StringLine::NextWord(line, start, ' ') == "flux")
			{
				
				
				while(start < (int)line.size())
					vFlux.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
				
				double Flux[vFlux.size()];
				for(int i=0; i < (int)vFlux.size();i++)
					Flux[i] = vFlux[i];
				
				fFlux = new TGraph(vTime.size(), Time, Flux);
				
			}
		}
		
		
		do
		{
			
			
			start = 0;
			int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
			int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
			
			if(A!=0 && Z!=0)
			{
				double DPQuantity[vTime.size()];
				for(int k = 0; k < (int)vTime.size(); k++ )
					DPQuantity[k] = 0;
				
				
				ZAI zaitmp(Z, A, I);
				int i=0;
				while(start < (int)line.size())
				{
					double DPQuantityTmp = atof(StringLine::NextWord(line, start, ' ').c_str());
					DPQuantity[i] = (double)DPQuantityTmp;
					i++;
					
				}
				fEvolutionData.insert(pair<ZAI ,TGraph* >(zaitmp, new TGraph((int)vTime.size()-1, Time, DPQuantity) ) );
			}
			
			getline(DecayDB, line);
			if(line == "" || line == "CrossSection" ) break;
		}while (!DecayDB.eof() );
		
		if(line == "CrossSection")
		{
			fIsCrossSection = true;
			getline(DecayDB, line);
			
			if (line == "Fission")
			{
				getline(DecayDB, line);
				
				do
				{
					double DPQuantity[vTime.size()];
					for(int k = 0; k < (int)vTime.size(); k++ )
						DPQuantity[k] = 0;
					
					start = 0;
					int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
					int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
					int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
					if(A!=0 && Z!=0)
					{
						
						
						ZAI zaitmp(Z, A, I);
						int i=0;
						while(start < (int)line.size())
						{
							long double DPQuantityTmp = atof(StringLine::NextWord(line, start, ' ').c_str());
							DPQuantity[i] = (double)DPQuantityTmp;
							i++;
							
						}
						fFissionXS.insert(pair<ZAI ,TGraph* >(zaitmp, new TGraph(vTime.size()-1, Time, DPQuantity) ) );
					}
					getline(DecayDB, line);
					if(line == "" || line == "Capture" ) break;
				}while (  !DecayDB.eof() );
			}
			
			if (line == "Capture")
			{
				getline(DecayDB, line); // Nuclei is given with "A Z"
				
				do
				{
					double DPQuantity[vTime.size()];
					for(int k = 0; k < (int)vTime.size(); k++ )
						DPQuantity[k] = 0;
					
					
					start = 0;
					int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
					int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
					int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
					
					if(A!=0 && Z!=0)
					{
						
						
						ZAI zaitmp(Z, A, I);
						int i=0;
						while(start < (int)line.size())
						{
							long double DPQuantityTmp = atof(StringLine::NextWord(line, start, ' ').c_str());
							DPQuantity[i] = (double)DPQuantityTmp;
							i++;
							
						}
						fCaptureXS.insert(pair<ZAI ,TGraph* >(zaitmp, new TGraph(vTime.size()-1, Time, DPQuantity) ) );
					}
					getline(DecayDB, line); // Nuclei is given with "A Z"
					if(line == "" || line == "n2n" ) break;
				}while ( !DecayDB.eof() );
				
			}
			
			if (line == "n2n")
			{
				
				getline(DecayDB, line); // Nuclei is given with "A Z"
				
				do
				{
					double DPQuantity[vTime.size()];
					for(int k = 0; k < (int)vTime.size(); k++ )
						DPQuantity[k] = 0;
					
					start = 0;
					int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
					int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
					int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
					
					if(A!=0 && Z!=0)
					{
						
						
						ZAI zaitmp(Z, A, I);
						int i=0;
						while(start < (int)line.size())
						{
							long double DPQuantityTmp = atof(StringLine::NextWord(line, start, ' ').c_str());
							DPQuantity[i] = (double)DPQuantityTmp;
							i++;
							
						}
						fn2nXS.insert(pair<ZAI ,TGraph* >(zaitmp, new TGraph(vTime.size()-1, Time, DPQuantity) ) );
					}
					getline(DecayDB, line); // Nuclei is given with "A Z"
					if(line == "" ) break;
					
				}while ( !DecayDB.eof() );
			}
			
		}
		DecayDB.close();
		start = 0;
		
		string InfoDBFile  = DBfile.erase(DBfile.size()-3,DBfile.size());
		InfoDBFile += "info";
		ifstream InfoDB(InfoDBFile.c_str());							// Open the File
		if(!InfoDB)
		{
			GetLog()->fLog << "!!Warning!! !!!EvolutionData!!! \n Can't open \"" << InfoDBFile << "\"\n" << endl;
			return;
		}
		
		start = 0;
		getline(InfoDB, line);
		if (StringLine::NextWord(line, start, ' ') == "Reactor")
			fReactorType =  StringLine::NextWord(line, start, ' ');
		start = 0;
		getline(InfoDB, line);
		if (StringLine::NextWord(line, start, ' ') == "Fueltype")
			fFuelType =  StringLine::NextWord(line, start, ' ');
		start = 0;
		getline(InfoDB, line);
		if (StringLine::NextWord(line, start, ' ') == "CycleTime")
			fCycleTime =  atof(StringLine::NextWord(line, start, ' ').c_str());;
		getline(InfoDB, line); // Assembly HM Mass
		start = 0;
		getline(InfoDB, line);
		if (StringLine::NextWord(line, start, ' ') == "ConstantPower")
			fPower =  atof(StringLine::NextWord(line, start, ' ').c_str());
		getline(InfoDB, line); // cutoff
		getline(InfoDB, line); // NUmber of Nuclei
		start = 0;
		getline(InfoDB, line);
		if (StringLine::NextWord(line, start, ' ') == "NormalizationFactor")
		{
			double NormFactor = atof(StringLine::NextWord(line, start, ' ').c_str());
			fPower = fPower * NormFactor;
			double Flux[vFlux.size()];
			for(int i=0; i < (int)vFlux.size();i++)
				Flux[i] = vFlux[i];
			
			fFlux = new TGraph(vTime.size()-1, Time, Flux);
		}
		start = 0;
		getline(InfoDB, line);
		if (StringLine::NextWord(line, start, ' ') == "FinalHeavyMetalMass")
			fHMMass =  atof(StringLine::NextWord(line, start, ' ').c_str());
		
	}
	
}





