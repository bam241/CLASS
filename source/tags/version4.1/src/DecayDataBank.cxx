#include "DecayDataBank.hxx"

#include "CLASSLogger.hxx"
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

DecayDataBank::DecayDataBank():CLASSObject(new CLASSLogger("DecayDataBank.log"))
{
	
	string  CLASSPATH = getenv("CLASS_PATH");
	string	DB_index_file = CLASSPATH + "/DATA_BASES/DECAY/ALL/Decay.idx";
	fDataBaseIndex = DB_index_file;
	fFastCalculation = true;
	
	// Warning
	INFO 	 <<  " A EvolutionData has been define :" << endl;
	INFO	 <<  "\t His index is : \"" << DB_index_file << "\"" << endl << endl;

}

//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________

//________________________________________________________________________

DecayDataBank::DecayDataBank(string DB_index_file):CLASSObject(new CLASSLogger("DecayDataBank.log"))
{

	fDataBaseIndex = DB_index_file;
	fFastCalculation = true;

	// Warning
	INFO 	 <<  " A EvolutionData has been define :" << endl;
	INFO	 <<  "\t His index is : \"" << DB_index_file << "\"" << endl << endl;

}
//________________________________________________________________________

DecayDataBank::DecayDataBank(CLASSLogger* log, string DB_index_file):CLASSObject(log)
{
	
	fDataBaseIndex = DB_index_file;
	fFastCalculation = true;
	
	// Warning
	INFO 	 <<  " A EvolutionData has been define :" << endl;
	INFO	 <<  "\t His index is : \"" << DB_index_file << "\"" << endl << endl;
	
}

//________________________________________________________________________
DecayDataBank::~DecayDataBank()
{
	
}

//________________________________________________________________________
IsotopicVector	DecayDataBank::Evolution(const ZAI& zai, double dt)
{
DBGL
	IsotopicVector	returnIV;
	
	map<ZAI ,EvolutionData >::iterator it = fDecayDataBank.find(zai);
	
	if (it ==  fDecayDataBank.end() )
	{
		ifstream DB_index(fDataBaseIndex.c_str());
		if( !DB_index)
		{
			ERROR << " Can't open \"" << fDataBaseIndex << "\"" << endl;
			exit (1);
		}
		bool zaifind = false;
		string tmp;
		getline(DB_index,tmp);							// Read first line
		while (!DB_index.eof())
		{
			string line;
			int start = 0;
			getline(DB_index,line);							// Read first line
			string first = StringLine::NextWord(line,start);				// Read first word
			
			if(first.size() == 0) break;						// If First word is null.... quit
			
			int rZ = atoi(first.c_str());						// Get Z
			int rA = atoi(StringLine::NextWord(line,start).c_str());			// Get A
			int rI = atoi(StringLine::NextWord(line,start).c_str());			// Get Isomeric State
			
			if(rZ ==  zai.Z() && rA ==  zai.A() && rI ==  zai.I() )
			{
				string file_name = StringLine::NextWord(line,start);
				EvolutionData evolutionproduct = EvolutionData(GetLog(), file_name, true);
#pragma omp critical(DBupdate)
				{fDecayDataBank.insert( pair<ZAI ,EvolutionData >(zai, evolutionproduct) );}
				returnIV = evolutionproduct.GetIsotopicVectorAt(dt);
				zaifind = true;
			}
		}
		
		if(!zaifind)
		{
			WARNING << " Oups... Can't Find the ZAI : "
			 <<  zai.Z() << " " << zai.A() << " "	 <<  zai.I() << "!!! It will be considered as stable !!" << endl;
			
			EvolutionData evolutionproduct = EvolutionData(GetLog()," " , true, zai);
			{fDecayDataBank.insert( pair<ZAI, EvolutionData >(zai, evolutionproduct) );}
			returnIV = evolutionproduct.GetIsotopicVectorAt(dt);
			
			
		}
		
		
	}
	else	returnIV = (*it).second.GetIsotopicVectorAt(dt);

DBGL
	return returnIV;
}

//________________________________________________________________________
bool DecayDataBank::IsDefine(const ZAI& zai) const
{
	
	map<ZAI ,EvolutionData > evolutiondb = (*this).GetDecayDataBank();
	if (evolutiondb.find(zai) != evolutiondb.end())
		return true;
	else
		return false;
	
}



//________________________________________________________________________
//	Get Decay
//________________________________________________________________________

IsotopicVector DecayDataBank::GetDecay(IsotopicVector isotopicvector, cSecond t)
{
	DBGL
	IsotopicVector IV;
	
	// Time slicing
	
	if( t > 1e16)
	{
		ERROR << " To long decay ... cut it !!! (max = 1e16 s ~ 300 000 000 years )" << endl;
		exit(1);
	}
	
	if(fFastCalculation)		// If using fast calculation simply interpolate between 2 timesteps
	{
		map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
		map<ZAI ,double >::iterator it;
		for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		{
			if((*it).second > 0)
			{
 			IsotopicVector ivtmp = Evolution(it->first, t) * (*it).second ;
				IV += ivtmp;
			}
		}
		
	}
	else			// If using precise calculation calculate the evolution at each digit
	{
		int evolutionDecade[17];
		cSecond remainingTime = t;
		for(int i = 16; i >= 0; i--)
		{
			evolutionDecade[i] = (cSecond)remainingTime/pow(10,i);
			remainingTime -=  evolutionDecade[i]*pow(10,i);
		}
		
		
		IV = isotopicvector;
		
		for (int i = 16; i >= 0; i--)
		{
			if(evolutionDecade[i] != 0)
			{
				map<ZAI ,double> isotopicquantity = IV.GetIsotopicQuantity();
				map<ZAI ,double >::iterator it;
				
				IV  = IsotopicVector();
				for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
					IV += Evolution(it->first, evolutionDecade[i]*pow(10,i) ) * (*it).second ;
			}
		}
	}
	DBGL
	return IV;
}
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
