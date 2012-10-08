#include "EvolutionDataBase.hxx"

#include "IsotopicVector.hxx"
#include "EvolutiveProduct.hxx"
#include "LogFile.hxx"
#include "Defines.hxx"
#include "StringLine.hxx"

#include <sstream>
#include <TGraph.h>
#include <TString.h>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>


using namespace std;
//________________________________________________________________________
//
//		EvolutionDataBase
//
//
//
//
//________________________________________________________________________

template<>
void EvolutionDataBase<IsotopicVector>::ReadDataBase();

template<>
EvolutionDataBase<ZAI>::EvolutionDataBase(LogFile* Log, string DB_index_file)
{
	DBGL;
	fLog = Log;
	fDataBaseIndex = DB_index_file;
	DBGL;
}
//________________________________________________________________________
template<>
EvolutionDataBase<ZAI>::~EvolutionDataBase()
{
	DBGL;
	DBGL;
}

template<>
IsotopicVector	EvolutionDataBase<ZAI>::Evolution(const ZAI& zai, double dt)
{
	DBGL;
	IsotopicVector	returnIV;

	map<ZAI ,EvolutiveProduct >::iterator it = fEvolutionDataBase.find(zai);

	if (it == fEvolutionDataBase.end() )
	{
		ifstream DB_index(fDataBaseIndex.c_str());
		if( !DB_index)
		{
			cout << "!!!EVOLUTIVE DB!!!! Can't open \"" << fDataBaseIndex << "\"" << endl;
			fLog->fLog << "!!!EVOLUTIVE DB!!!! Can't open \"" << fDataBaseIndex << "\"" << endl;
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

			if(rZ == zai.Z() && rA == zai.A() && rI == zai.I() )
			{
				string file_name = StringLine::NextWord(line,start);
				EvolutiveProduct evolutionproduct = EvolutiveProduct(fLog, file_name);
				#pragma omp critical(DBupdate)
				{fEvolutionDataBase.insert( pair<ZAI ,EvolutiveProduct >(zai, evolutionproduct) );}
				returnIV = evolutionproduct.GetIsotopicVectorAt(dt);
				zaifind = true;
			}
		}
	
		if(zaifind == false) 
		{
			{
			fLog->fLog << "!!Warning!! !!!EVOLUTIVE DB!!! Oups... Can't Find the ZAI : " 
				   << zai.Z() << " " << zai.A() << " "	<< zai.I() << "!!! It will be considered as stable !!" << endl;
				EvolutiveProduct evolutionproduct = EvolutiveProduct(fLog," " , true, zai);
				{fEvolutionDataBase.insert( pair<ZAI, EvolutiveProduct >(zai, evolutionproduct) );}
				returnIV = evolutionproduct.GetIsotopicVectorAt(dt);

			}
		}
	

	}
	else	returnIV = (*it).second.GetIsotopicVectorAt(dt);
	
	return returnIV;
}

template<>
bool EvolutionDataBase<ZAI>::IsDefine(const ZAI& zai) const
{
	DBGL;
	map<ZAI ,EvolutiveProduct > evolutiondb = (*this).GetEvolutionDataBase();
	if (evolutiondb.find(zai) != evolutiondb.end()) 
		return true;
	else	
		return false;

}
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________

template<>
EvolutionDataBase<IsotopicVector>::~EvolutionDataBase()
{
	DBGL;
	DBGL;
}

template<>
EvolutionDataBase<IsotopicVector>::EvolutionDataBase(LogFile* Log, string DB_index_file)
{
	DBGL;
	fLog = Log;
	fDataBaseIndex = DB_index_file;
	ReadDataBase();
	DBGL;
}

template<>
void EvolutionDataBase<IsotopicVector>::ReadDataBase()
{
	DBGL;
	ifstream DataDB(fDataBaseIndex.c_str());							// Open the File
	if(!DataDB)
	{
		cout << "!!Warning!! !!!EvolutionDataBase!!! \n Can't open \"" << fDataBaseIndex << "\"\n" << endl;
		fLog->fLog << "!!Warning!! !!!EvolutionDataBase!!! \n Can't open \"" << fDataBaseIndex << "\"\n" << endl;
	}
	vector<double> vTime;
	vector<double> vTimeErr;

	string line;
	int start = 0;
	
	
	
	// First Get Fuel Type
	getline(DataDB, line); 
	if( StringLine::NextWord(line, start, ' ') != "TYPE") 
	{
		cout << "!!Bad Trouble!! !!!EvolutionDataBase!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the type of the DataBase"<< endl;
		fLog->fLog << "!!Bad Trouble!! !!!EvolutionDataBase!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the type of the DataBase"<< endl;
		exit (1);
	}
	fFuelType = StringLine::NextWord(line, start, ' ');
	// First Get Fuel Parameter
	getline(DataDB, line); 
	start = 0;
	if( StringLine::NextWord(line, start, ' ') != "PARAM") 
	{
		cout << "!!Bad Trouble!! !!!EvolutionDataBase!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the Parameter of the DataBase"<< endl;
		fLog->fLog << "!!Bad Trouble!! !!!EvolutionDataBase!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the Parameter of the DataBase"<< endl;
		exit (1);
	}
	while(start < (int)line.size())
		fFuelParameter.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));

	
	//Then Get All the Database

	while (!DataDB.eof())
	{
		getline(DataDB, line); 
		if(line != "")
		{
			EvolutiveProduct evolutionproduct = EvolutiveProduct(fLog, line);
			IsotopicVector ivtmp  = evolutionproduct.GetIsotopicVectorAt(0.);
			if(ivtmp.GetZAIIsotopicQuantity(94, 242, 0) == 0 ) 
				ivtmp.Print();
			fEvolutionDataBase.insert( pair<IsotopicVector, EvolutiveProduct >(ivtmp , evolutionproduct) );
		}
	}
	DBGL;
}
//________________________________________________________________________



template<>
map<double, EvolutiveProduct> EvolutionDataBase<IsotopicVector>::GetDistances(IsotopicVector isotopicvector) const
{
	DBGL;
	map<double, EvolutiveProduct> distances;
	map<IsotopicVector, EvolutiveProduct >::iterator it;
	map<IsotopicVector, EvolutiveProduct > evolutiondb = fEvolutionDataBase;
	
	
	for( it = evolutiondb.begin(); it != evolutiondb.end(); it++ )
	{
		pair<map<double, EvolutiveProduct>::iterator, bool> IResult;
		IResult = distances.insert( pair<double, EvolutiveProduct>( RelativDistance(isotopicvector, (*it).first / Norme( (*it).first )*Norme(isotopicvector) ), (*it).second ) );
	}

	return distances;
	DBGL;
}

