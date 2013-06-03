#include "Storage.hxx"

#include "EvolutionDataBase.hxx"
#include "CLASS.hxx"
#include "Defines.hxx"
#include "LogFile.hxx"


#include <sstream>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

//________________________________________________________________________
//
//		Storage
//
//
//
//
//________________________________________________________________________
ClassImp(Storage)

Storage::Storage()
{
DBGL;
DBGL;
}

//________________________________________________________________________
Storage::Storage(EvolutionDataBase<ZAI>* evolutivedb)
{
DBGL;
	fInternalTime = 0;
	fDecayDataBase = evolutivedb;
DBGL;
}

//________________________________________________________________________
Storage::~Storage()
{
DBGL;
DBGL;
}


//________________________________________________________________________
//	Get Decay
//________________________________________________________________________
IsotopicVector Storage::GetDecay(IsotopicVector isotopicvector, cSecond t)
{
DBGL;
	IsotopicVector IV;

	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	{
		if((*it).second > 0)
		{
 			IsotopicVector ivtmp = fDecayDataBase->Evolution(it->first, t) * (*it).second ;
			IV += ivtmp;
		}
	}
DBGL;
	return IV;
DBGL;
}

//________________________________________________________________________
void Storage::ClearStock()
{
DBGL;
	IsotopicVector EmptyIV;
	fIVFullStock = EmptyIV;
	fIVStock.clear();
DBGL;
}

//________________________________________________________________________
void Storage::AddToStock(IsotopicVector isotopicvector)
{
DBGL;
	if(fParc->GetStockManagement() == true)
		fIVStock.push_back(isotopicvector);
	AddToFullStock(isotopicvector);
DBGL;
}

//________________________________________________________________________
void Storage::TakeFractionFromStock(int IVId,double fraction)
{
DBGL;
	if(fParc->GetStockManagement() == true)
	{
		if(fraction > 1 || fraction < 0)
		{
			cout << fraction << endl;
			cout << "!!Warning!! !!!Storage!!! You try to remove fraction superior than 1 or a negative one..." << endl;
			fLog->fLog << "!!Warning!! !!!Storage!!! You try to remove fraction superior than 1 or a negative one..." << endl;
		}
		else 
		{
			fIVFullStock -= fIVStock[IVId]*fraction;
			fIVStock[IVId] = fIVStock[IVId]*(1-fraction);
			
		}

	}
	else
	{
		cout << "!!Warning!! !!!Storage!!! TakeFractionFromStock can't be DEFINE without REAL stock management" << endl;
		fLog->fLog << "!!Warning!! !!!Storage!!! TakeFractionFromStock can't be DEFINE without REAL stock management" << endl;
		exit(1);

	}
	
	
DBGL;
}

void Storage::TakeFromStock(IsotopicVector isotopicvector)
{
DBGL;
	if(fParc->GetStockManagement() == false)
		fIVFullStock -= isotopicvector;
	else
	{
		cout << "!!Warning!! !!!Storage!!! TakeFromStock can't be DEFINE WITH REAL stock management" << endl;
		fLog->fLog << "!!Warning!! !!!Storage!!! TakeFromStock can't be DEFINE WITH REAL stock management" << endl;
		exit(1);
	}
DBGL;
}

//________________________________________________________________________
void Storage::StorageEvolution(cSecond t)
{
DBGL;

	if(t == fInternalTime && t !=0 ) return;

	for(int i = (int)fIVStock.size()-1 ; i >=0; i--)
		if(Norme(fIVStock[i]) == 0)
			fIVStock.erase(fIVStock.begin()+i); 
	
	

	int EvolutionTime = t - fInternalTime;

	fIVFullStock = 	GetDecay(fIVFullStock , EvolutionTime);

#pragma omp parallel for
	for (int i=0; i <(int) fIVStock.size() ; i++)
	{
		fIVStock[i] = GetDecay(fIVStock[i] , EvolutionTime);
	}
	

	
DBGL;
}

//________________________________________________________________________
void Storage::Evolution(cSecond t)
{
DBGL;
	// Check if the Storage has been created ...
	if(t == fInternalTime && t!=0) return;
	// Make the evolution for the Storage ...
	StorageEvolution(t);
	// ... And Finaly update the AbsoluteInternalTime
	fInternalTime = t;
	

DBGL;
}

void Storage::Write(string filename, cSecond date)
{
DBGL;
	for(int i=0;i < (int)fIVStock.size(); i++)
	{
		
		fIVStock[i].Write(filename, date);
	}
DBGL;
}

