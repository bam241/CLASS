#include "TreatmentFactory.hxx"
#include "EvolutionDataBase.hxx"
#include "LogFile.hxx"
#include "Defines.hxx"


#include <sstream>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

//________________________________________________________________________
//
//		TreatmentFactory
//
//
//
//
//________________________________________________________________________
ClassImp(TreatmentFactory)

TreatmentFactory::TreatmentFactory()
{
DBGL;
DBGL;
}

//________________________________________________________________________
TreatmentFactory::TreatmentFactory(EvolutionDataBase* evolutivedb,
				   long int creation,
				   long int coolingtime ,
				   long int separationtime)
{
DBGL;
	fStockManagement = true;
	fCoolingTime = coolingtime;
	fInternalTime = 0;
	fDecayDataBase = evolutivedb;
	fSeparationTime = separationtime;
	fCreationTime = creation;
	IsStarted = false;
	fCoolingLastIndex = 0;
DBGL;
}

//________________________________________________________________________
TreatmentFactory::~TreatmentFactory()
{
DBGL;
DBGL;
}



//________________________________________________________________________
void	TreatmentFactory::AddValorisableIV(ZAI zai, double factor)
{
DBGL;
	pair<map<ZAI, double>::iterator, bool> IResult;
	if(factor > 1) factor = 1;
	
	if(factor > 0)
	{
		IResult = fValorisableIV.insert( pair<ZAI ,double>(zai, factor));
		if(IResult.second == false)
			IResult.first->second = factor;
	}
DBGL;
}


//________________________________________________________________________
//	Get Decay
//________________________________________________________________________
IsotopicVector TreatmentFactory::GetDecayProduct(IsotopicVector isotopicvector, long int t)
{
DBGL;
	IsotopicVector IV;
	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	for(map<ZAI ,double >::iterator it = isotopicquantity.begin(); 
			it != isotopicquantity.end(); it++)
	{
		if((*it).second > 0)
		{
 			IsotopicVector ivtmp = fDecayDataBase->DecayProduction(it->first, t) * (*it).second ;
			IV += ivtmp;
		}
	}
DBGL;
	return IV;
DBGL;
}


//________________________________________________________________________
//	Add Temporary IV : 
//		Cooling
//		Separation
//________________________________________________________________________
void TreatmentFactory::AddIVCooling(IsotopicVector IV)
{ 
DBGL;
	fIVCooling.push_back(IV);
	fCoolingStartingTime.push_back(fInternalTime);
	fCoolingLastIndex++;
	fCoolingIndex.push_back(fCoolingLastIndex);
DBGL;
}


//________________________________________________________________________
void TreatmentFactory::RemoveIVCooling(int i)		//!< Remove a Cooling IsotopicVector
{
DBGL;
	fIVCooling.erase(fIVCooling.begin()+i);
	fCoolingStartingTime.erase(fCoolingStartingTime.begin()+i);
	fCoolingIndex.erase(fCoolingIndex.begin()+i); 
DBGL;
}

void TreatmentFactory::ClearIVStock()
{
DBGL;
	IsotopicVector EmptyIV;
	fIVFullStock = EmptyIV;
	fIVStock.clear();
DBGL;
}


void TreatmentFactory::AddIVStock(IsotopicVector isotopicvector)
{
DBGL;

	if(fParc->GetStockManagement() == true) fIVStock.push_back(isotopicvector);
	AddIVFullStock(isotopicvector);

DBGL;
}


//________________________________________________________________________
//	Time Action with the reste of the Universe : 
//		In/Out stock
//		Evolution : 
//			Waste, Stock, Separating, Cooling
//		Separation
//________________________________________________________________________
void TreatmentFactory::TakeFromStock(IsotopicVector isotopicvector, int index)
{
DBGL;

	if(fStockManagement == true )
		fIVStock[index] -= isotopicvector;

	fIVFullStock -= isotopicvector;


DBGL;
}

//________________________________________________________________________
pair<IsotopicVector, IsotopicVector> TreatmentFactory::Separation(IsotopicVector isotopicvector)
{
DBGL;
	//[0] = stock ; [1] = waste
	pair<IsotopicVector, IsotopicVector>	IVTmp;
	
	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	for(map<ZAI ,double >::iterator it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	{
		map<ZAI ,double>::iterator it2;
		it2 = fValorisableIV.find((*it).first);
		if ( it2 != fValorisableIV.end() )
		{
			IVTmp.first.Add(	(*it).first, (*it).second * (*it2).second );		//stock
			IVTmp.second.Add(	(*it).first, (*it).second * (1-(*it2).second) );	//waste
		}
		else IVTmp.second.Add(	(*it).first, (*it).second );	//waste
	}
DBGL;
	return IVTmp;
}


//________________________________________________________________________
void TreatmentFactory::WasteEvolution(long int t)
{
DBGL;
	// Check if the TF has been created ...
	if(t<fCreationTime) return;

	long int EvolutionTime = t - fInternalTime;
	fIVWaste = GetDecayProduct(fIVWaste , EvolutionTime);
DBGL;
}

//________________________________________________________________________
void TreatmentFactory::StockEvolution(long int t)
{
DBGL;
	if(t == fInternalTime && t!=0) return;
	int EvolutionTime = t - fInternalTime;
	fIVFullStock = GetDecayProduct(fIVFullStock , EvolutionTime);

#pragma omp parallel for
	for (int i=0; i <(int) fIVStock.size() ; i++)
	{
		fIVStock[i] = GetDecayProduct(fIVStock[i] , EvolutionTime);
	}
	
	
DBGL;
}

//________________________________________________________________________
void TreatmentFactory::CoolingEvolution(long int t)
{
DBGL;
	if(t == fInternalTime && t!=0) return;
	int i;
	int RemainingCoolingTime;
	long int EvolutionTime = t - fInternalTime;
#pragma omp parallel for
	for ( i = 0 ; i < (int)fIVCooling.size() ; i++)
	{
		if (t - fCoolingStartingTime[i] >= fCoolingTime) // ">" should not append, only "=" is normal...
		{
			if (t - fCoolingStartingTime[i] > fCoolingTime) // Warning & Quit
			{
				cout		<< "!!Warning!! !!!TreamtmentFactory!!! Cooling Step : " << t/3600./24/365.25<< " :"
						<< " An evolution Step is probably missing ! " << " " << endl;
				cout << t << " " << fCoolingStartingTime[i] << " " << fCoolingTime << endl;
				fLog->fLog 	<< "!!Warning!! !!!TreamtmentFactory!!! Cooling Step : "<< t << " :"
						<< " An evolution Step is probably missing ! " << endl;
				exit (1);
			}

			RemainingCoolingTime = fCoolingTime - (fInternalTime - fCoolingStartingTime[i]);
			//Cooling Decay
			fIVCooling[i] = GetDecayProduct( fIVCooling[i], RemainingCoolingTime);

#pragma omp critical(DeleteCoolingIVPB)
			{fCoolingEndOfCycle.push_back(i);}

		}
		else if ( fCoolingStartingTime[i] != t )
		{
			fIVCooling[i] = GetDecayProduct( fIVCooling[i] , EvolutionTime);
		}
	}

	sort (fCoolingEndOfCycle.begin(), fCoolingEndOfCycle.end());
DBGL;
}


//________________________________________________________________________
void TreatmentFactory::Evolution(long int t)
{
DBGL;
	// Check if the TF has been created ...
	if(t<fCreationTime) return;
	if(t == fInternalTime && t!=0) return;
	
	if(fInternalTime == 0 && IsStarted == false)
	{
		fInternalTime = fCreationTime;
		IsStarted = true;
	}
	// Make the evolution for the Waste ...
	WasteEvolution(t);
	
	// ... the Stock ...
	StockEvolution(t);
	
	// ... Then Deal the Cooling IV ...
	CoolingEvolution(t);
	
	// ... And Finaly update the AbsoluteInternalTime
	fInternalTime = t;
	

DBGL;
}


//________________________________________________________________________
void TreatmentFactory::Dump()
{
//------ Cooling ------//
	for(int i = (int)fCoolingEndOfCycle.size()-1; i >=0 ; i--)	// IV End Of Cooling
	{
		int idx = fCoolingEndOfCycle[i];			// Get Index number
		AddIVStock(fIVCooling[idx]);
		
		fCoolingEndOfCycle.erase(fCoolingEndOfCycle.begin()+i);	// Remove index entry
		RemoveIVCooling(idx);					// Remove IVcooling
	}
	
	if((int)fCoolingEndOfCycle.size() != 0 )// Control
	{
		cout		<< "Problem while Dumping Cooling"<< endl;
		fLog->fLog 	<< "Problem while Dumping Cooling"<< endl;
		exit (1);
	}
}


