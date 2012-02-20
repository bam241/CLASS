#include "CLASSHeaders.hxx"

using namespace std;
//________________________________________________________________________
//
//		EvolutionDataBase
//
//
//
//
//________________________________________________________________________
EvolutionDataBase::EvolutionDataBase(string DB_index_file)
{
	DBGL;
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
//	if{dt > fDatabaseEndTime}
	map<ZAI ,EvolutiveProduct* >::iterator it;
	it = fEvolutionDataBase.find(zai);
	if (it == fEvolutionDataBase.end() )
	{
		EvolutiveProduct* evolutionproduct = new EvolutiveProduct(zai.Z(), zai.A(), zai.I(), fDataBaseIndex );
		fEvolutionDataBase.insert( pair<ZAI ,EvolutiveProduct *>(zai, evolutionproduct) );
		return evolutionproduct->GetIsotopicVectorAt(dt);
	}
	else	return (*it).second->GetIsotopicVectorAt(dt);

}

//________________________________________________________________________
bool EvolutionDataBase::IsDefine(const ZAI& zai) const
{
	DBGL;
	map<ZAI ,EvolutiveProduct* > evolutiondb = (*this).GetEvolutionDataBase();
	if (  evolutiondb.find(zai) != evolutiondb.end() ) 
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
		pair<ZAI ,EvolutiveProduct *>(zai, new EvolutiveProduct(zai.Z(), zai.A(), zai.I(), GetDataBaseIndex())) );
	
	DBGL;
	return IResult.second;
	
}

//________________________________________________________________________
void EvolutionDataBase::Print() const 
{


}


