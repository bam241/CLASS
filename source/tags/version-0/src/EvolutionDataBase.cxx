#include "EvolutionDataBase.hxx"
#include "IsotopicVector.hxx"
#include "EvolutiveProduct.hxx"
#include "LogFile.hxx"
#include "Defines.hxx"

#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;
//________________________________________________________________________
//
//		EvolutionDataBase
//
//
//
//
//________________________________________________________________________

EvolutionDataBase::EvolutionDataBase(LogFile* Log, string DB_index_file)
{
	DBGL;
	fLog = Log;
	fDataBaseIndex = DB_index_file;
	DBGL;
}

//________________________________________________________________________
EvolutionDataBase::~EvolutionDataBase()
{	
	DBGL;
}


//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
IsotopicVector	EvolutionDataBase::DecayProduction(const ZAI& zai, double dt)
{
	DBGL;
	IsotopicVector	returnIV;
	map<ZAI ,EvolutiveProduct* >::iterator it;
		it = fEvolutionDataBase.find(zai);

		if (it == fEvolutionDataBase.end() )
		{
		
			EvolutiveProduct* evolutionproduct = new EvolutiveProduct(fLog, zai.Z(), zai.A(), zai.I(), fDataBaseIndex);
			#pragma omp critical(DBupdate)
			{fEvolutionDataBase.insert( pair<ZAI ,EvolutiveProduct *>(zai, evolutionproduct) );}
			returnIV = evolutionproduct->GetIsotopicVectorAt(dt);
		}
		else	returnIV = (*it).second->GetIsotopicVectorAt(dt);
	
	return returnIV;
}

//________________________________________________________________________
bool EvolutionDataBase::IsDefine(const ZAI& zai) const
{
	DBGL;
	map<ZAI ,EvolutiveProduct* > evolutiondb = (*this).GetEvolutionDataBase();
	if (evolutiondb.find(zai) != evolutiondb.end()) 
		return true;
	else	
		return false;

}

//________________________________________________________________________
bool EvolutionDataBase::AddEvolutiveProduct(const ZAI& zai)
{
	DBGL;
	pair<map<ZAI, EvolutiveProduct *>::iterator, bool> IResult;
	
	IResult = fEvolutionDataBase.insert(
		pair<ZAI ,EvolutiveProduct *>(zai, new EvolutiveProduct(fLog, zai.Z(), zai.A(), zai.I(), GetDataBaseIndex())) );
	
	DBGL;
	return IResult.second;
	
}

//________________________________________________________________________
void EvolutionDataBase::Print() const 
{


}


