#include "CLASSNucleiFiliation.hxx"
#include "ZAI.hxx"
#include "IsotopicVector.hxx"

#include <map>
#include <vector>
#include <cmath>

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

CLASSNucleiFiliation::CLASSNucleiFiliation():CLASSObject()
{
}
CLASSNucleiFiliation::CLASSNucleiFiliation(CLASSLogger* log):CLASSObject(log)
{
}

CLASSNucleiFiliation::CLASSNucleiFiliation(const CLASSNucleiFiliation& CNF):CLASSObject()
{
	fNucleiFiliation = CNF.GetNucleiFIliation();
}
//________________________________________________________________________
CLASSNucleiFiliation::~CLASSNucleiFiliation()
{
	fNucleiFiliation.clear() ;
}
//________________________________________________________________________
void CLASSNucleiFiliation::Add( ZAI Mother, IsotopicVector Daughter )
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
	DBGV(Mother.Z() << " " << Mother.A() << " " << Mother.I());
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
				
				if(FastDecayChain.GetQuantity(-1, -1, -1) == 0) // Check if the FastDecayChain is known
					it_Filiation->second += Daughter_BR * FastDecayChain; // Add the FastDecayCHain & apply the BR for the cutted Daughter
				else
				{
					
					ZAI Mother = DautherList[i];
					while (FastDecayChain.GetQuantity(-1, -1, -1) != 0 && GoodNuclei.find(Mother) == GoodNuclei.end())
					{
						Mother = GetArtificialDecay(Mother);			// Do an Artifial decay on the nuclei
						FastDecayChain = CuttedNuclei.GetFiliation(Mother); // Get the fast decay chain of it
					}
					
					if(GoodNuclei.find(Mother) != GoodNuclei.end())
						it_Filiation->second += Mother * Daughter_BR;
					
					else if ( FastDecayChain.GetQuantity(-1, -1, -1) == 0)
						it_Filiation->second += FastDecayChain * Daughter_BR;
					
					else
					{
						ERROR << "Problem in Articial Decay!! check it!!" << endl;
						exit(1);
					}

				}
			}
		}
		
	}
	DBGL
}
//________________________________________________________________________

ZAI CLASSNucleiFiliation::GetArtificialDecay(ZAI Mother)
{
	DBGL

	int A = Mother.A();
	int Z = Mother.Z();
	int I = Mother.I();
	
	if(I!=0)
		return ZAI(Z,A,I-1);
	else
	{
		//Coef Ac & As of Bette & Weisacker are approximativ but precise enough for this application....
		double Ac = 0.695;
		double As = 23.2;
		
		double ZTh = A/2 * ( 1 )/ ( 1 + Ac / (4*As) * pow(A,2/3) );  // Stable Z from isobarn calculation using Bette & Weisacker formula.
	
		
		if( Z > ZTh )		// Then Beta+
			return ZAI(Z-1,A,I);
		else			// Then Beta-
			return ZAI(Z+1,A,I);
		
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
		int count = 0;
		map<ZAI, IsotopicVector> CopyfNucleiFiliation = fNucleiFiliation;
		map<ZAI, IsotopicVector>::iterator it_Filiation;
		for(it_Filiation = fNucleiFiliation.begin(); it_Filiation != fNucleiFiliation.end(); it_Filiation++) // Loop on the mother ZAI
		{
			vector<ZAI> DautherList = it_Filiation->second.GetZAIList(); // Get the list of daughter ZAI
			
			for (int i = 0; i < (int)DautherList.size(); i++)		//Loop on daughter
			{
				if(GoodNuclei.find(DautherList[i]) == GoodNuclei.end() ) // if the ZAI is not in a dealed nuclei (cutted or unknown)
				{	count++;
					map<ZAI, IsotopicVector>::iterator it_FiliationCopy = CopyfNucleiFiliation.find(it_Filiation->first)  ;

					double Daughter_BR = it_FiliationCopy->second.GetQuantity(DautherList[i]);	// Get the quantity of the ZAI
					it_FiliationCopy->second -= Daughter_BR * DautherList[i];			// Remove it from the daughter list					
					IsotopicVector FastDecayChain = (*this).GetFiliation(DautherList[i]); // Get the fast decay chain of it

					if(FastDecayChain.GetQuantity(-1, -1, -1) == 0) // Check if the FastDecayChain is known
						it_FiliationCopy->second += Daughter_BR * FastDecayChain; // Add the FastDecayCHain & apply the BR for the cutted Daughter
					else
					{					
						ZAI Mother = DautherList[i];
						while (FastDecayChain.GetQuantity(-1, -1, -1) != 0 && GoodNuclei.find(Mother) == GoodNuclei.end())
						{
							Mother = GetArtificialDecay(Mother);			// Do an Artifial decay on the nuclei
							FastDecayChain = (*this).GetFiliation(Mother); // Get the fast decay chain of it
						}

						if(GoodNuclei.find(Mother) != GoodNuclei.end())
							it_FiliationCopy->second += Mother * Daughter_BR;
						
						else if ( FastDecayChain.GetQuantity(-1, -1, -1) == 0)
							it_FiliationCopy->second += FastDecayChain * Daughter_BR;
						
						else
						{
							ERROR << "Problem in Articial Decay!! check it!!" << endl;
							exit(1);
						}

					}
				}
			}
			
		}
		fNucleiFiliation = CopyfNucleiFiliation;
		Cleaned = true;
		
		for(it_Filiation = fNucleiFiliation.begin(); it_Filiation != fNucleiFiliation.end(); it_Filiation++) // Loop on the mother ZAI
		{
			vector<ZAI> DautherList = it_Filiation->second.GetZAIList(); // Get the list of daughter ZAI
			
			for (int i = 0; i < (int)DautherList.size(); i++)		//Loop on daughter
				if(GoodNuclei.find(DautherList[i]) == GoodNuclei.end() ) // if the ZAI is not in a dealed nuclei (cutted or unknown)
					Cleaned = false;
		}

	}
	
	NormalizeBranchingRatio();
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


//________________________________________________________________________

vector<ZAI> CLASSNucleiFiliation::GetZAIList() const
{
	
	map<ZAI ,IsotopicVector > IsotopicQuantity = GetNucleiFIliation();
	map<ZAI ,IsotopicVector >::iterator it;
	vector<ZAI> zailist;
	for( it = IsotopicQuantity.begin(); it != IsotopicQuantity.end(); it++)
		zailist.push_back( (*it).first );
	
	return zailist;
	
}



















































