#include "EvolutiveProduct.hxx"
#include "IsotopicVector.hxx"
#include "ZAI.hxx"
#include "LogFile.hxx"
#include "Defines.hxx"
#include "StringLine.hxx"

#include <TGraphErrors.h>
#include <TString.h>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <TGraphErrors.h>

using namespace std;
//________________________________________________________________________
//
//		EvolutiveProduct
//
//
//
//
//________________________________________________________________________

EvolutiveProduct::EvolutiveProduct(LogFile* Log, int Z, int A, int I, string DBindexfile)
{
DBGL;
	fLog = Log;
	ifstream DB_index(DBindexfile.c_str());
	if( !DB_index)
	{
		cout << "!!!EVOLUTIVE DB!!!! Can't open \"" << DBindexfile << "\"" << endl;
		fLog->fLog << "!!!EVOLUTIVE DB!!!! Can't open \"" << DBindexfile << "\"" << endl;
		exit (1);
	}
	bool zaifind = false;		
	string tmp;		
	getline(DB_index,tmp);							// Read first line
	while (!DB_index.eof()) 
	{
		string line;			
		int start=0;				
		getline(DB_index,line);							// Read first line
		string first=StringLine::NextWord(line,start);				// Read first word
		
		if(first.size()==0) break;						// If First word is null.... quit
		
		int rZ=atoi(first.c_str());						// Get Z
		int rA=atoi(StringLine::NextWord(line,start).c_str());			// Get A
		int rI=atoi(StringLine::NextWord(line,start).c_str());			// Get Isomeric State

		if(rZ == Z && rA == A && rI == I)
		{
			string file_name = StringLine::NextWord(line,start);
			ReadDB(file_name);		// Read Decay produc DB file name
			zaifind = true;
		}
	}
	if(zaifind == false) 
	{
		#pragma omp critical(LOGupdate)
		{
		fLog->fLog << "!!Warning!! !!!EVOLUTIVE DB!!! Oups... Can't Find the ZAI : " 
			   << Z << " " << A << " "	<< I << "!!! It will be considered as stable !!" << endl;
		AddAsStable(Z, A, I);
		}
	}
DBGL;

}

//________________________________________________________________________

EvolutiveProduct::EvolutiveProduct(string DB_reactor_file)
{
	DBGL;
	string file_name = DB_reactor_file;
	ReadDB(file_name);		// Read Evolution Produc DB file name
	DBGL;

}

//________________________________________________________________________
EvolutiveProduct::~EvolutiveProduct()
{
	DBGL;
	DBGL;
}

//________________________________________________________________________
//________________________________________________________________________

void EvolutiveProduct::AddAsStable(int Z,int  A,int I)
{
	DBGL;
	double time[2] = {0, (int)(500*365.4*3600*24)};
	double Err[2]	= {0, 0};
	double quantity[2] = {1., 1.};
	ZAI zaitmp(Z, A, I);
	
	fEvolutiveProduct.insert(pair<ZAI ,TGraphErrors* >(zaitmp, new TGraphErrors(2, time, quantity, Err, Err) ) );
	DBGL;
}
//________________________________________________________________________

//________________________________________________________________________
void EvolutiveProduct::ReadDB(string DBfile)
{
	DBGL;
	ifstream DecayDB(DBfile.c_str());							// Open the File
	if(!DecayDB)
	{
		cout << "!!Warning!! !!!EvolutiveProduct!!! \n Can't open \"" << DBfile << "\"\n" << endl;
		fLog->fLog << "!!Warning!! !!!EvolutiveProduct!!! \n Can't open \"" << DBfile << "\"\n" << endl;
	}
	vector<double> vTime;
	vector<double> vTimeErr;

	string line;
	int start = 0;
	
	getline(DecayDB, line); // Nuclei is given with "A Z"
	if( StringLine::NextWord(line, start, ' ') != "time") 
	{
		cout << "!!Bad Trouble!! !!!EvolutiveProduct!!! Bad Database file : " <<  DBfile << endl;
		fLog->fLog << "!!Bad Trouble!! !!!EvolutiveProduct!!! Bad Database file : " <<  DBfile << endl;
		exit (1);
	}
	
	while(start < (int)line.size())
	{
		vTime.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
		vTimeErr.push_back(0); // Set to Zero because not yet error in DB
		
	}

	double Time[vTime.size()];
	double TimeErr[vTime.size()];
	for(int i=0; i < (int)vTime.size();i++)
		{Time[i] = vTime[i]; TimeErr[i] = vTimeErr[i];} 


	while (!DecayDB.eof())
	{
		double DPQuantity[vTime.size()];
		double DPQuantityErr[vTime.size()];

		getline(DecayDB, line); // Nuclei is given with "A Z"
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
				DPQuantityErr[i]=0; // Set to Zero because not yet error in DB
				i++;

			}
			fEvolutiveProduct.insert(pair<ZAI ,TGraphErrors* >(zaitmp, new TGraphErrors(vTime.size(), Time, DPQuantity, TimeErr, DPQuantityErr) ) );
		}

	}
	DBGL;
}

//________________________________________________________________________
Double_t EvolutiveProduct::Interpolate(double t, TGraphErrors& EvolutionGraph)
{
	DBGL;
	TString fOption; 
	return (double)EvolutionGraph.Eval(t,0x0,fOption);
	DBGL;
}

//________________________________________________________________________
TGraphErrors*	EvolutiveProduct::GetEvolutionTGraphErrors(const ZAI& zai)
{
	DBGL;
	map<ZAI ,TGraphErrors *>::iterator it = GetEvolutiveProduct().find(zai) ;
	
	if ( it != GetEvolutiveProduct().end() ) 
		return it->second;	
	else
		return new TGraphErrors();
	
	DBGL;
}

//________________________________________________________________________
IsotopicVector	EvolutiveProduct::GetIsotopicVectorAt(double t)
{
	DBGL;
	IsotopicVector IsotopicVectorTmp;
	map<ZAI ,TGraphErrors* >::iterator it;
	for( it = fEvolutiveProduct.begin(); it != fEvolutiveProduct.end(); it++ )
	{
		IsotopicVectorTmp.Add( (*it).first, Interpolate(t, *((*it).second)) );
	}

	DBGL;
	return IsotopicVectorTmp;
}

//________________________________________________________________________






