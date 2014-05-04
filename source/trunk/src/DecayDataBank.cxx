#include "DecayDataBank.hxx"

#include "IsotopicVector.hxx"
#include "CLASSHeaders.hxx"
#include "LogFile.hxx"
#include "StringLine.hxx"

#include <TGraph.h>
#include <TString.h>


#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>


//________________________________________________________________________
//
//		DecayDataBank
//
//
//
//
//________________________________________________________________________

DecayDataBank::DecayDataBank()
{
}

//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________

DecayDataBank::DecayDataBank(LogFile* Log, string DB_index_file, bool setlog, bool olfreadmethod)
{
	
	SetLog(Log);
	fDataBaseIndex = DB_index_file;
	
	fOldReadMethod = olfreadmethod;
	
	// Warning
	if(IsLog())
	{
		cout	<< "!!INFO!! !!!DecayDataBank!!! A EvolutionData has been define :" << endl;
		cout	<< "\t His index is : \"" << DB_index_file << "\"" << endl << endl;
		
		GetLog()->fLog 	<< "!!INFO!! !!!DecayDataBank!!! A EvolutionData has been define :" << endl;
		GetLog()->fLog	<< "\t His index is : \"" << DB_index_file << "\"" << endl << endl;
	}
	
}

//________________________________________________________________________
DecayDataBank::~DecayDataBank()
{
	
}

//________________________________________________________________________
IsotopicVector	DecayDataBank::Evolution(const ZAI& zai, double dt)
{
	
	IsotopicVector	returnIV;
	
	map<ZAI ,EvolutionData >::iterator it = fDecayDataBank.find(zai);
	
	if (it == fDecayDataBank.end() )
	{
		ifstream DB_index(fDataBaseIndex.c_str());
		if( !DB_index)
		{
			cout << "!!!EVOLUTIVE DB!!!! Can't open \"" << fDataBaseIndex << "\"" << endl;
			GetLog()->fLog << "!!!EVOLUTIVE DB!!!! Can't open \"" << fDataBaseIndex << "\"" << endl;
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
				EvolutionData evolutionproduct = EvolutionData(GetLog(), file_name);
#pragma omp critical(DBupdate)
				{fDecayDataBank.insert( pair<ZAI ,EvolutionData >(zai, evolutionproduct) );}
				returnIV = evolutionproduct.GetIsotopicVectorAt(dt);
				zaifind = true;
			}
		}
		
		if(zaifind == false)
		{
			GetLog()->fLog << "!!Warning!! !!!EVOLUTIVE DB!!! Oups... Can't Find the ZAI : " ;
			GetLog()->fLog << zai.Z() << " " << zai.A() << " "	<< zai.I() << "!!! It will be considered as stable !!" << endl;
			
			EvolutionData evolutionproduct = EvolutionData(GetLog()," " , false, zai);
			{fDecayDataBank.insert( pair<ZAI, EvolutionData >(zai, evolutionproduct) );}
			returnIV = evolutionproduct.GetIsotopicVectorAt(dt);
			
			
		}
		
		
	}
	else	returnIV = (*it).second.GetIsotopicVectorAt(dt);
	
	return returnIV;
}

bool DecayDataBank::IsDefine(const ZAI& zai) const
{
	
	map<ZAI ,EvolutionData > evolutiondb = (*this).GetDecayDataBank();
	if (evolutiondb.find(zai) != evolutiondb.end())
		return true;
	else
		return false;
	
}

//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
