#include "EvolutionData.hxx"
#include "CLASSMethod.hxx"

#include "CLASSLogger.hxx"
#include "CLASSConstante.hxx"
#include "external/StringLine.hxx"

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


double 	Distance(IsotopicVector IV1, EvolutionData Evd1 )
{
	
	double d2 = 0;
	IsotopicVector IV2 = Evd1.GetIsotopicVectorAt(0.);
	
	IsotopicVector IVtmp = IV1;
	map<ZAI ,double> IVtmpIsotopicQuantity = IVtmp.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	
	double SumOfXs = 0;
	for( it = IVtmpIsotopicQuantity.begin(); it != IVtmpIsotopicQuantity.end(); it++)
	{
		
		SumOfXs += Evd1.GetXSForAt(0., (*it).first, 1);
		SumOfXs += Evd1.GetXSForAt(0., (*it).first, 2);
		SumOfXs += Evd1.GetXSForAt(0., (*it).first, 3);
		
	}
	
	
	for( it = IVtmpIsotopicQuantity.begin(); it != IVtmpIsotopicQuantity.end(); it++)
	{
		double Z1 = 0.0;
		double Z2 = 0.0;
		double XS = 0.0;
		
		Z1 = IV1.GetZAIIsotopicQuantity( (*it).first );
		Z2 = IV2.GetZAIIsotopicQuantity( (*it).first );
		XS = Evd1.GetXSForAt(0., (*it).first, 1)
		+ Evd1.GetXSForAt(0., (*it).first, 2)
		+ Evd1.GetXSForAt(0., (*it).first, 3);
		
		d2 += pow(Z1-Z2 , 2 ) * pow(XS,2);
		
	}
	
	return sqrt(d2)/SumOfXs;
	
}

double 	Distance(EvolutionData Evd1, IsotopicVector IV1 )
{
	return Distance(IV1,Evd1);
}
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
EvolutionData operator*(EvolutionData const& evol, double F)
{
	
	EvolutionData evoltmp;
	
	map<ZAI ,TGraph* > EvolutionData = evol.GetInventoryEvolution();
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
	evoltmp.SetFissionXS(evol.GetFissionXS());
	evoltmp.SetCaptureXS(evol.GetCaptureXS());
	evoltmp.Setn2nXS(evol.Getn2nXS());
	
	
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

EvolutionData Multiply(EvolutionData const& evol, double F)
{
	
	EvolutionData evoltmp;
	map<ZAI ,TGraph* > EvolutionData = evol.GetInventoryEvolution();
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
	
	
	EvolutionData = evol.GetFissionXS();
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
		evoltmp.FissionXSInsert( pair<ZAI, TGraph*> ( (*it).first,new TGraph((*it).second->GetN(), X, Y) ) );
		
	}
	
	EvolutionData = evol.GetCaptureXS();
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
		evoltmp.CaptureXSInsert( pair<ZAI, TGraph*> ( (*it).first,new TGraph((*it).second->GetN(), X, Y) ) );
		
	}
	EvolutionData = evol.Getn2nXS();
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
		evoltmp.n2nXSInsert( pair<ZAI, TGraph*> ( (*it).first,new TGraph((*it).second->GetN(), X, Y) ) );
		
	}
	
	evoltmp.SetPower(evol.GetPower()*F);
	
	
	
	return evoltmp;
	
}

//________________________________________________________________________
EvolutionData Multiply(double F, EvolutionData const& evol)
{
	
	return Multiply(evol,F);
	
}

EvolutionData Sum(EvolutionData const& evol1, EvolutionData const& evol2)
{
	EvolutionData EvolSum = evol1;
	map<ZAI ,TGraph* > EvolutionData1 = EvolSum.GetInventoryEvolution();
	map<ZAI ,TGraph* > EvolutionData2 = evol2.GetInventoryEvolution();
	map<ZAI ,TGraph* >::iterator it;
	
	for(it = EvolutionData2.begin(); it != EvolutionData2.end(); it++)
	{
		pair<map<ZAI, TGraph*>::iterator, bool> IResult;
		
		IResult  = EvolutionData1.insert( pair<ZAI, TGraph*> ( *it ) );
		if(!(IResult.second) )
		{
			double X[(*it).second->GetN()];
			double Y[(*it).second->GetN()];
			map<ZAI ,TGraph* >::iterator it2 = EvolutionData1.find( (*it).first );
			
			
			for(int i = 0; i < (*it).second->GetN(); i++)
			{
				double y;
				(*it).second->GetPoint( i, X[i], y );
				Y[i] = y + (*it2).second->Eval(X[i]);
			}
			IResult.first->second = new TGraph((*it).second->GetN(), X, Y);
			
		}
		
	}
	EvolSum.SetInventoryEvolution(EvolutionData1);
	
	
	EvolutionData1 = evol1.GetFissionXS();
	EvolutionData2 = evol2.GetFissionXS();
	
	for(it = EvolutionData2.begin(); it != EvolutionData2.end(); it++)
	{
		pair<map<ZAI, TGraph*>::iterator, bool> IResult;
		
		IResult  = EvolutionData1.insert( pair<ZAI, TGraph*> ( *it ) );
		
		if(!(IResult.second) )
		{
			double X[(*it).second->GetN()];
			double Y[(*it).second->GetN()];
			map<ZAI ,TGraph* >::iterator it2 = EvolutionData1.find( (*it).first );
			
			
			for(int i = 0; i < (*it).second->GetN(); i++)
			{
				double y;
				(*it).second->GetPoint( i, X[i], y );
				Y[i] = y + (*it2).second->Eval(X[i]);
			}
			IResult.first->second = new TGraph((*it).second->GetN(), X, Y);
		}
	}
	EvolSum.SetFissionXS(EvolutionData1);
	
	
	EvolutionData1 = EvolSum.GetCaptureXS();
	EvolutionData2 = evol2.GetCaptureXS();
	
	for(it = EvolutionData2.begin(); it != EvolutionData2.end(); it++)
	{
		pair<map<ZAI, TGraph*>::iterator, bool> IResult;
		
		IResult  = EvolutionData1.insert( pair<ZAI, TGraph*> ( *it ) );
		
		if(!(IResult.second) )
		{
			double X[(*it).second->GetN()];
			double Y[(*it).second->GetN()];
			map<ZAI ,TGraph* >::iterator it2 = EvolutionData1.find( (*it).first );
			
			
			for(int i = 0; i < (*it).second->GetN(); i++)
			{
				double y;
				(*it).second->GetPoint( i, X[i], y );
				Y[i] = y + (*it2).second->Eval(X[i]);
			}
			IResult.first->second = new TGraph((*it).second->GetN(), X, Y);
		}
	}
	EvolSum.SetCaptureXS(EvolutionData1);
	
	
	EvolutionData1 = EvolSum.Getn2nXS();
	EvolutionData2 = evol2.Getn2nXS();
	
	for(it = EvolutionData2.begin(); it != EvolutionData2.end(); it++)
	{
		pair<map<ZAI, TGraph*>::iterator, bool> IResult;
		
		IResult  = EvolutionData1.insert( pair<ZAI, TGraph*> ( *it ) );
		
		if(!(IResult.second) )
		{
			double X[(*it).second->GetN()];
			double Y[(*it).second->GetN()];
			map<ZAI ,TGraph* >::iterator it2 = EvolutionData1.find( (*it).first );
			
			
			for(int i = 0; i < (*it).second->GetN(); i++)
			{
				double y;
				(*it).second->GetPoint( i, X[i], y );
				Y[i] = y + (*it2).second->Eval(X[i]);
			}
			IResult.first->second = new TGraph((*it).second->GetN(), X, Y);
		}
	}
	EvolSum.Setn2nXS(EvolutionData1);
	
	return EvolSum;
}




//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________


ClassImp(EvolutionData)



EvolutionData::EvolutionData():CLASSObject()
{
	fIsCrossSection = false;
	fPower = 0;
	fCycleTime = 0;
	fKeff = 0;
	fFlux = 0;
	fDB_file = "";
}

//________________________________________________________________________
EvolutionData::EvolutionData(CLASSLogger* log):CLASSObject(log)
{
	fIsCrossSection = false;
	fPower = 0;
	fCycleTime = 0;
	fKeff = 0;
	fFlux = 0;
	fDB_file = "";
}

//________________________________________________________________________
EvolutionData::EvolutionData(CLASSLogger* log, string DB_file, bool isDecay, ZAI zai):CLASSObject(log)
{
	
	fIsCrossSection = false;
	fDB_file = DB_file;
	fPower = 0;
	fCycleTime = 0;
	fKeff = 0;
	fFlux = 0;
	fisDecay = isDecay;
	
	if(zai != ZAI(0,0,0))
		AddAsStable(zai);
	else
		ReadDB( fDB_file, false);		// Read Evolution Produc DB file name
	
	
	
	
}


//________________________________________________________________________
EvolutionData::EvolutionData(bool oldread, CLASSLogger* log, string DB_file, bool isDecay, ZAI zai):CLASSObject(log)
{
	
	fIsCrossSection = false;
	fDB_file = DB_file;
	fPower = 0;
	fCycleTime = 0;
	fKeff = 0;
	fFlux = 0;
	fisDecay = isDecay;
	
	if(zai != ZAI(0,0,0))
		AddAsStable(zai);
	else
		ReadDB( fDB_file, oldread);		// Read Evolution Produc DB file name
	
	
	
	
}

//________________________________________________________________________
EvolutionData::~EvolutionData()
{

}
//________________________________________________________________________
void EvolutionData::DeleteEvolutionData()
{
	
	map<ZAI ,TGraph* >::iterator it_del;
	
	for( it_del = fInventoryEvolution.begin(); it_del != fInventoryEvolution.end(); it_del++)
	{
		delete (*it_del).second;
		(*it_del).second = 0;
	}
	for( it_del = fFissionXS.begin(); it_del != fFissionXS.end(); it_del++)
	{
		delete (*it_del).second;
		(*it_del).second = 0;
	}
	for( it_del = fCaptureXS.begin(); it_del != fCaptureXS.end(); it_del++)
	{
		delete (*it_del).second;
		(*it_del).second = 0;
	}
	for( it_del = fn2nXS.begin(); it_del != fn2nXS.end(); it_del++)
	{
		delete (*it_del).second;
		(*it_del).second = 0;
	}
	
	
	delete	fKeff;
	delete	fFlux;
	
	fInventoryEvolution.clear();
	fFissionXS.clear();
	fCaptureXS.clear();
	fn2nXS.clear();
	
	fFlux = 0;
	fKeff = 0;

	
}


//________________________________________________________________________
void EvolutionData::DeleteEvolutionDataCopy()
{
	
	if(fDB_file == "")
	{
		map<ZAI ,TGraph* >::iterator it_del;
		
		for( it_del = fInventoryEvolution.begin(); it_del != fInventoryEvolution.end(); it_del++)
		{
			delete (*it_del).second;
			(*it_del).second = 0;
		}
		for( it_del = fFissionXS.begin(); it_del != fFissionXS.end(); it_del++)
		{
			delete (*it_del).second;
			(*it_del).second = 0;
		}
		for( it_del = fCaptureXS.begin(); it_del != fCaptureXS.end(); it_del++)
		{
			delete (*it_del).second;
			(*it_del).second = 0;
		}
		for( it_del = fn2nXS.begin(); it_del != fn2nXS.end(); it_del++)
		{
			delete (*it_del).second;
			(*it_del).second = 0;
		}
		
		
		delete	fKeff;
		delete	fFlux;
		
		fInventoryEvolution.clear();
		fFissionXS.clear();
		fCaptureXS.clear();
		fn2nXS.clear();
		
		fFlux = 0;
		fKeff = 0;
	}
	
	
}

bool EvolutionData::NucleiInsert(pair<ZAI, TGraph*> zaitoinsert)
{
	
	pair<map<ZAI, TGraph*>::iterator, bool> IResult;
	IResult = fInventoryEvolution.insert( zaitoinsert);
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
	
	double time[2] = {0, (500*cYear)};
	double quantity[2] = {1., 1.};
	
	fInventoryEvolution.insert(pair<ZAI ,TGraph* >(zai, new TGraph(2, time, quantity) ) );
	
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
	
	map<ZAI ,TGraph *>::iterator it = GetInventoryEvolution().find(zai) ;
	
	if ( it != GetInventoryEvolution().end() )
		return it->second;
	else
		return new TGraph();
	
	
}

//________________________________________________________________________
IsotopicVector	EvolutionData::GetIsotopicVectorAt(double t)
{
	
	IsotopicVector IsotopicVectorTmp;
	map<ZAI ,TGraph* >::iterator it;
	for( it = fInventoryEvolution.begin(); it != fInventoryEvolution.end(); it++ )
	{
		IsotopicVectorTmp.Add( (*it).first, Interpolate(t, *((*it).second)) );
	}
	
	
	return IsotopicVectorTmp;
}

//________________________________________________________________________
double	EvolutionData::GetXSForAt(double t, ZAI zai, int ReactionId)
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
		default:ERROR << " Wrong ReactionId !!" << endl;
			exit(1);
	}
	
	map<ZAI ,TGraph* >::iterator it = XSEvol.find(zai);
	
	
	if (it ==  XSEvol.end())
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
	
	if(!fisDecay)
		ReadInfo();							// Read the .info associeted
	
	ifstream DecayDB(DBfile.c_str());	// Open the File
	if(!DecayDB)				//check if file is correctly open
	{
		ERROR << " \n Can't open \"" << DBfile << "\"\n" << endl;
		ERROR << "\t-> Hint : If loading .dat files using a .idx file (like for a decay data base)\nmake sure the paths in it are correct" << endl;;
		exit(1);
	}
	vector<double> vTime;
	
	string line;
	int start = 0;
	
	getline(DecayDB, line);
	if( tlc(StringLine::NextWord(line, start, ' ')) != "time")
	{
		ERROR << " Bad Database file : " <<  DBfile << endl;
		ERROR << " The first Line MUST be the time line !!!" << endl;
		exit (1);
	}
	
	while(start < (int)line.size())
		vTime.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
	
	fFinalTime = vTime.back();
	double Time[vTime.size()];
	for(int i = 0; i < (int)vTime.size();i++)
		Time[i] = vTime[i];
	const int NTimeStep = sizeof(Time)/sizeof(double);
	
	
	enum { Keff = 1, Flux, Inv, XSFis, XSCap, XSn2n };
	
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
	
	fHeavyMetalMass = cZAIMass.GetMass( GetIsotopicVectorAt(0.).GetActinidesComposition() );
	
	DecayDB.close();
	
	
}

void EvolutionData::ReadKeff(string line, double* time, int NTimeStep)
{
	if(fKeff != 0)
		delete fKeff;
	
	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "keff" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"keff\" not found !" << endl;
		exit(1);
	}
	
	
	double Keff[NTimeStep];
	
	int i = 0;
	while(start < (int)line.size() && i < NTimeStep )		// Read the Data
	{
		Keff[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
		i++;
	}
	
	
	fKeff = new TGraph(NTimeStep, time, Keff);			// Add the TGraph
	
	
	
	
}

void EvolutionData::ReadFlux(string line, double* time, int NTimeStep)
{
	
	if(fFlux != 0)
		delete fFlux;
	
	int start = 0;
	
	if( tlc(StringLine::NextWord(line, start, ' ')) != "flux" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"flux\" not found !" << endl;
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
	
	
	
}


void	EvolutionData::ReadInv(string line, double* time, int NTimeStep)
{
	
	
	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "inv" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"inv\" not found !" << endl;
		exit(1);
	}
	// Read the Z A I
	int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
	
	if(A != 0 && Z != 0)
	{
		double Inv[NTimeStep];
		
		int i = 0;
		while(start < (int)line.size() && i < NTimeStep )	// Read the Data
		{
			Inv[i] = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
			i++;
		}
		// Add the TGraph
		fInventoryEvolution.insert(pair<ZAI ,TGraph* >(ZAI(Z,A,I), new TGraph(NTimeStep, time, Inv) ) );
	}
	
	
}


void	EvolutionData::ReadXSFis(string line, double* time, int NTimeStep)
{
	
	
	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "xsfis" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"xsfis\" not found !" << endl;
		exit(1);
	}
	// Read the Z A I
	int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
	
	if(A != 0 && Z != 0)
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
		ERROR << " Bad keyword : \"xscap\" not found !" << endl;
		exit(1);
	}
	// Read the Z A I
	int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
	
	if(A != 0 && Z != 0)
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
		ERROR << " Bad keyword : \"xsn2n\" not found !" << endl;
		exit(1);
	}
	// Read the Z A I
	int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
	
	if(A != 0 && Z != 0)
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
	string InfoDBFile = fDB_file.erase(fDB_file.size()-3,fDB_file.size());
	InfoDBFile += "Info";
	ifstream InfoDB_tmp(InfoDBFile.c_str());				// Open the File
	
	if(!InfoDB_tmp)
	{
		InfoDBFile  = InfoDBFile.erase(InfoDBFile.size()-4,InfoDBFile.size());
		InfoDBFile += "info";
		
	}
	InfoDB_tmp.close();
	
	ifstream InfoDB(InfoDBFile.c_str());				// Open the File
	if(!InfoDB)
	{
		WARNING << "!!ERROR!! !!!EvolutionData!!! \n Can't open \"" << InfoDBFile << "\"\n" << endl;
	}
	
	int start = 0;
	string line;
	getline(InfoDB, line);
	if ( tlc(StringLine::NextWord(line, start, ' ')) ==  "reactor")
		fReactorType = StringLine::NextWord(line, start, ' ');
	
	start = 0;
	getline(InfoDB, line);
	if (tlc(StringLine::NextWord(line, start, ' ')) ==  "fueltype")
		fFuelType = StringLine::NextWord(line, start, ' ');
	start = 0;
	getline(InfoDB, line);
	if ( tlc(StringLine::NextWord(line, start, ' ')) ==  "cycletime")
		fCycleTime = atof(StringLine::NextWord(line, start, ' ').c_str());;
	getline(InfoDB, line); // Assembly HM Mass
	start = 0;
	getline(InfoDB, line);
	if ( tlc(StringLine::NextWord(line, start, ' ')) ==  "constantpower")
		fPower = atof(StringLine::NextWord(line, start, ' ').c_str());
	InfoDB.close();
}

void EvolutionData::OldReadDB(string DBfile)
{
	
	
	ifstream DecayDB(DBfile.c_str());							// Open the File
	if(!DecayDB)
	{
		ERROR << " Can't open \"" << DBfile << "\"\n" << endl;
		ERROR  << "\t-> Hint : If loading .dat files using a .idx file (like for a decay data base)\nmake sure the paths in it are correct";
		
		exit(1);
	}
	vector<double> vTime;
	vector<double> vTimeErr;
	
	string line;
	int start = 0;
	
	getline(DecayDB, line);
	if( StringLine::NextWord(line, start, ' ') != "time")
	{
		ERROR << " Bad Database file : " <<  DBfile << endl;
		exit (1);
	}
	
	while(start < (int)line.size())
		vTime.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
	
	fFinalTime = vTime.back();
	double Time[vTime.size()];
	for(int i = 0; i < (int)vTime.size();i++)
		Time[i] = vTime[i];
	vector<double> vFlux;
	start = 0;
	getline(DecayDB, line);
	string tmp = StringLine::NextWord(line, start, ' ');
	if ( tmp ==  "keff"  || tmp ==  "Keff" )
	{
		vector<double> vKeff;
		while(start < (int)line.size())
			vKeff.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
		
		double Keff[vKeff.size()];
		for(int i = 0; i < (int)vKeff.size();i++)
			Keff[i] = vKeff[i];
		
		fKeff = new TGraph(vTime.size(), Time, Keff);
		
		start = 0;
		getline(DecayDB, line);
		if (StringLine::NextWord(line, start, ' ') ==  "flux")
		{
			
			
			while(start < (int)line.size())
				vFlux.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
			
			double Flux[vFlux.size()];
			for(int i = 0; i < (int)vFlux.size();i++)
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
		
		if(A != 0 && Z != 0)
		{
			double DPQuantity[vTime.size()];
			for(int k = 0; k < (int)vTime.size(); k++ )
				DPQuantity[k] = 0;
			
			
			ZAI zaitmp(Z, A, I);
			int i = 0;
			while(start < (int)line.size())
			{
				double DPQuantityTmp = atof(StringLine::NextWord(line, start, ' ').c_str());
				DPQuantity[i] = (double)DPQuantityTmp;
				i++;
				
			}
			TGraph* tgraphtmp = new TGraph((int)vTime.size()-1, Time, DPQuantity);
			fInventoryEvolution.insert(pair<ZAI ,TGraph* >(zaitmp, tgraphtmp) );
		}
		
		getline(DecayDB, line);
		if(line ==  "" || line ==  "CrossSection" ) break;
	}while (!DecayDB.eof() );
	
	if(line ==  "CrossSection")
	{
		fIsCrossSection = true;
		getline(DecayDB, line);
		
		if (line ==  "Fission")
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
				if(A != 0 && Z != 0)
				{
					
					
					ZAI zaitmp(Z, A, I);
					int i = 0;
					while(start < (int)line.size())
					{
						long double DPQuantityTmp = atof(StringLine::NextWord(line, start, ' ').c_str());
						DPQuantity[i] = (double)DPQuantityTmp;
						i++;
						
					}
					fFissionXS.insert(pair<ZAI ,TGraph* >(zaitmp, new TGraph(vTime.size()-1, Time, DPQuantity) ) );
				}
				getline(DecayDB, line);
				if(line ==  "" || line ==  "Capture" ) break;
			}while (  !DecayDB.eof() );
		}
		
		if (line ==  "Capture")
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
				
				if(A != 0 && Z != 0)
				{
					
					
					ZAI zaitmp(Z, A, I);
					int i = 0;
					while(start < (int)line.size())
					{
						long double DPQuantityTmp = atof(StringLine::NextWord(line, start, ' ').c_str());
						DPQuantity[i] = (double)DPQuantityTmp;
						i++;
						
					}
					fCaptureXS.insert(pair<ZAI ,TGraph* >(zaitmp, new TGraph(vTime.size()-1, Time, DPQuantity) ) );
				}
				getline(DecayDB, line); // Nuclei is given with "A Z"
				if(line ==  "" || line ==  "n2n" ) break;
			}while ( !DecayDB.eof() );
			
		}
		
		if (line ==  "n2n")
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
				
				if(A != 0 && Z != 0)
				{
					
					
					ZAI zaitmp(Z, A, I);
					int i = 0;
					while(start < (int)line.size())
					{
						long double DPQuantityTmp = atof(StringLine::NextWord(line, start, ' ').c_str());
						DPQuantity[i] = (double)DPQuantityTmp;
						i++;
						
					}
					fn2nXS.insert(pair<ZAI ,TGraph* >(zaitmp, new TGraph(vTime.size()-1, Time, DPQuantity) ) );
				}
				getline(DecayDB, line); // Nuclei is given with "A Z"
				if(line ==  "" ) break;
				
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
		INFO << " Can't open \"" << InfoDBFile << "\"\n" << endl;
		return;
	}
	
	start = 0;
	getline(InfoDB, line);
	if (StringLine::NextWord(line, start, ' ') ==  "Reactor")
		fReactorType = StringLine::NextWord(line, start, ' ');
	start = 0;
	getline(InfoDB, line);
	if (StringLine::NextWord(line, start, ' ') ==  "Fueltype")
		fFuelType = StringLine::NextWord(line, start, ' ');
	start = 0;
	getline(InfoDB, line);
	if (StringLine::NextWord(line, start, ' ') ==  "CycleTime")
		fCycleTime = atof(StringLine::NextWord(line, start, ' ').c_str());;
	getline(InfoDB, line); // Assembly HM Mass
	start = 0;
	getline(InfoDB, line);
	if (StringLine::NextWord(line, start, ' ') ==  "ConstantPower")
		fPower = atof(StringLine::NextWord(line, start, ' ').c_str());
	getline(InfoDB, line); // cutoff
	getline(InfoDB, line); // NUmber of Nuclei
	start = 0;
	getline(InfoDB, line);
	if (StringLine::NextWord(line, start, ' ') ==  "NormalizationFactor")
	{
		double NormFactor = atof(StringLine::NextWord(line, start, ' ').c_str());
		fPower = fPower * NormFactor;
		double Flux[vFlux.size()];
		for(int i = 0; i < (int)vFlux.size();i++)
			Flux[i] = vFlux[i];
		
		fFlux = new TGraph(vTime.size()-1, Time, Flux);
	}
	InfoDB.close();
}


//________________________________________________________________________
void EvolutionData::Print(string filename)
{	

	map<ZAI ,TGraph* >::iterator iterator;
	ofstream out(filename.c_str());


	out<<"time ";
	for(int t = 0 ; t < fInventoryEvolution.begin()->second->GetN() ; t++ )
			out<<fInventoryEvolution.begin()->second->GetX()[t]<< " ";
	out<<endl;
	
	if(fFlux)	
	{	
		out<<"flux ";
		for(int t = 0 ; t < fFlux->GetN() ; t++ )
			out<<fFlux->GetY()[t]<< " ";
		out<<endl;
	}

	if(fKeff)	
	{	out<<"keff ";
		for(int t = 0 ; t < fKeff->GetN() ; t++ )
			out<<fKeff->GetY()[t]<< " ";
		out<<endl;
	}
	for( iterator = fInventoryEvolution.begin(); iterator != fInventoryEvolution.end(); iterator++)
	{
		int N = (*iterator).second->GetN();

		if((*iterator).second->GetY()[N-1] != 0)
		{	out<<"Inv "<<(*iterator).first.Z()<<" "<<(*iterator).first.A()<<" "<<(*iterator).first.I()<<" ";

		for(int t = 0 ; t < N ; t++ )
			out<<(*iterator).second->GetY()[t]<< " ";
		out<<endl;
		}
	}

	for( iterator = fFissionXS.begin(); iterator != fFissionXS.end(); iterator++)
	{
		int N = (*iterator).second->GetN();
		if((*iterator).second->GetY()[N-1] != 0)
		{	out<<"XSFis "<<(*iterator).first.Z()<<" "<<(*iterator).first.A()<<" "<<(*iterator).first.I()<<" ";

		for(int t = 0 ; t < (*iterator).second->GetN() ; t++ )
			out<<-(*iterator).second->GetY()[t]*1e24<< " ";

		out<<endl;
		}
	}

	for( iterator = fCaptureXS.begin(); iterator != fCaptureXS.end(); iterator++)
	{
		int N = (*iterator).second->GetN();
		if((*iterator).second->GetY()[N-1] != 0)
		{			out<<"XSCap "<<(*iterator).first.Z()<<" "<<(*iterator).first.A()<<" "<<(*iterator).first.I()<<" ";

		for(int t = 0 ; t < (*iterator).second->GetN() ; t++ )
			out<<-(*iterator).second->GetY()[t]*1e24<< " ";

		out<<endl;
		}
	}

	for( iterator = fn2nXS.begin(); iterator != fn2nXS.end(); iterator++)
	{
		int N = (*iterator).second->GetN();
		if((*iterator).second->GetY()[N-1] != 0)
		{			out<<"XSn2n "<<(*iterator).first.Z()<<" "<<(*iterator).first.A()<<" "<<(*iterator).first.I()<<" ";

		for(int t = 0 ; t < (*iterator).second->GetN() ; t++ )
			out<<-(*iterator).second->GetY()[t]*1e24<< " ";

		out<<endl;
		}
	}


	
}

