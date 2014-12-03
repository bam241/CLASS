#include "CLASSNucleiFiliation.hxx"
#include "ZAI.hxx"
#include "IsotopicVector.hxx"

#include <map>
#include <vector>
#include "stdlib.h"

using namespace std;

//const string DEFAULTDATABASE = "DecayBase.dat";
//________________________________________________________________________
//
//		CLASSNucleiFiliation
//
//
//
//
//________________________________________________________________________
//____________________________InClass Operator____________________________
//________________________________________________________________________
ClassImp(CLASSNucleiFiliation)


CLASSNucleiFiliation::CLASSNucleiFiliation():CLASSObject()
{
}


CLASSNucleiFiliation::CLASSNucleiFiliation(CLASSNucleiFiliation CNF):CLASSObject()
{
	fNucleiFIliation = CNF.GetNucleiFIliation();
}



//________________________________________________________________________
CLASSNucleiFiliation::~CLASSNucleiFiliation()
{
	
}


//________________________________________________________________________
void CLASSNucleiFiliation::AddDaughterToZAI(ZAI Mother, IsotopicVector Daughter )
{
	DBGL
	
	pair< map<ZAI, IsotopicVector>::iterator, bool> IResult;
	
	
	IResult = fNucleiFiliation.insert( pair<ZAI, IsotopicVector> ( Mother, Daughter ) );
	
	if( !IResult.second)
		(*IResult.first).second += Daughter ;
	
	
	DBGL
}


//________________________________________________________________________
IsotopicVector CLASSNucleiFiliation::GetFiliation(ZAI Mother)
{
	DBGL
	map<ZAI, IsotopicVector>::iterator it_Filiation;
	
	it_Filiation = fNucleiFiliation.find(Mother);	// search for the ZAI in the map

	DBGL
	if(it_Filiation != fNucleiFiliation.end())	// test if it is present in the map
		return it_Filiation->second;
	else
		return ZAI(-1,-1,-1)*1;		// return -1 -1 _1 ZAI if unknown nuclei....
	
}

//________________________________________________________________________
void CLASSNucleiFiliation::FiliationCleanUp(map<ZAI, int> GoodNuclei, CLASSNucleiFiliation CuttedNuclei)
{
	DBGL
	map<ZAI, IsotopicVector>::iterator it_Filiation;
	for(it_Filiation = fNucleiFiliation.begin(); it_Filiation != fNucleiFiliation.end(); it_Filiation++) // loop on the filiation map (on the Mother ZAI)
	{
		vector<ZAI> DautherList = it_Filiation->second.GetZAIList();		// recover for each mother ZAI, the list of its daughter ZAI
		
		for (int i = 0; i < (int)DautherList.size(); i++)			// Loop on the Daughter ZAI
		{
			if(GoodNuclei.find(DautherList[i]) == GoodNuclei.end() ) // if the ZAI is not in a dealed nuclei (cutted or unknown)
			{
				double Daughter_BR = it_Filiation->second.GetQuantity(DautherList[i]);	// Get the quantity of the ZAI
				it_Filiation->second -= Daughter_BR * DautherList[i];			// Remove it from the daughter list
				
				
				IsotopicVector FastDecayChain = CuttedNuclei.GetFiliation(DautherList[i]); // Get the fast decay chain of it
				
				if(FastDecayChain.GetQuantity(-1, -1, -1) != 0) // Check if the FastDecayChain is known
					it_Filiation->second += Daughter_BR * FastDecayChain; // Add the FastDecayCHain & apply the BR for the cutted Daughter
				else
					it_Filiation->second += Daughter_BR * ZAI(-3,-3,-3); // Add a TMP nuclei the daughter nuclei is not known at all...
				
			}
		}
		
	}
	DBGL
}

//________________________________________________________________________
void CLASSNucleiFiliation::SelfFiliationCleanUp(map<ZAI, int> GoodNuclei)
{
	DBGL

	bool Cleaned = false;
	
	while (!Cleaned)	// loop until all the map is cleaned (all cut have been done)
	{
		
		Cleaned = true;
	
		map<ZAI, IsotopicVector>::iterator it_Filiation;
		for(it_Filiation = fNucleiFiliation.begin(); it_Filiation != fNucleiFiliation.end(); it_Filiation++) // Loop on the mother ZAI
		{
			vector<ZAI> DautherList = it_Filiation->second.GetZAIList(); // Get the list of daughter ZAI
			
			for (int i = 0; i < (int)DautherList.size(); i++)		//Loop on daughter
			{
				if(GoodNuclei.find(DautherList[i]) == GoodNuclei.end() ) // if the ZAI is not in a dealed nuclei (cutted or unknown)
				{
					Cleaned = false;
					
					double Daughter_BR = it_Filiation->second.GetQuantity(DautherList[i]);	// Get the quantity of the ZAI
					it_Filiation->second -= Daughter_BR * DautherList[i];			// Remove it from the daughter list
					
					
					IsotopicVector FastDecayChain = (*this).GetFiliation(DautherList[i]); // Get the fast decay chain of it
					
					if(FastDecayChain.GetQuantity(-1, -1, -1) != 0) // Check if the FastDecayChain is known
						it_Filiation->second += Daughter_BR * FastDecayChain; // Add the FastDecayCHain & apply the BR for the cutted Daughter
					else
						it_Filiation->second += Daughter_BR * ZAI(-3,-3,-3); // Add a TMP nuclei the daughter nuclei is not known at all...
					
				}
			}
			
		}
	
	}
	DBGL
}



//________________________________________________________________________
void CLASSNucleiFiliation::NormalizeBranchingRatio(double Value)
{
	DBGL
	map<ZAI, IsotopicVector>::iterator it_Filiation;
	for(it_Filiation = fNucleiFiliation.begin(); it_Filiation != fNucleiFiliation.end(); it_Filiation++)
	{
		it_Filiation->second *= Value/it_Filiation->second.GetSumOfAll();
	}
	
	DBGL
}


//________________________________________________________________________
void CLASSNucleiFiliation::NormalizeBranchingRatio(ZAI Mother, double Value)
{
	DBGL
	map<ZAI, IsotopicVector>::iterator it_Filiation = fNucleiFiliation.find(Mother);
	
	if( it_Filiation != fNucleiFiliation.end())
		it_Filiation->second *= Value/it_Filiation->second.GetSumOfAll();
	else
		WARNING << "Trying to normaliza a Branching Ratio for a Mother wich are not present in the Filiatiuon List...." << endl;
	
	DBGL
}

























































