#include "Storage.hxx"

#include "Scenario.hxx"
#include "CLASSLogger.hxx"


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




Storage::Storage():CLASSBackEnd(-1)
{
	DBGL
	SetIsStorageType();
	SetName("S_Storage.");

	DBGL
}

Storage::Storage(CLASSLogger* log):CLASSBackEnd(log, -1)
{
	DBGL
	SetIsStorageType();

	SetName("S_Storage.");

	INFO << " A new Storage has been define." << endl;
	
}
//________________________________________________________________________
Storage::Storage(CLASSLogger* log, DecayDataBank* evolutivedb):CLASSBackEnd(log, -1)
{
	DBGL
	SetIsStorageType();

	SetDecayDataBank(evolutivedb);
	
	SetName("S_Storage.");

	INFO << " A new Storage has been define." << endl;


	DBGL
}

//________________________________________________________________________
Storage::~Storage()
{


}

//________________________________________________________________________
void Storage::AddIV(IsotopicVector isotopicvector)
{

	AddCumulativeIVIn(isotopicvector);

	if(GetParc())
	{
		if(GetParc()->GetStockManagement() )
		{
			fIVArray.push_back(isotopicvector);
			fIVArrayArrivalTime.push_back(fInternalTime);
		}
	}
	else
	{
		fIVArray.push_back(isotopicvector);
		fIVArrayArrivalTime.push_back(fInternalTime);
	}
	UpdateInsideIV();
}

//________________________________________________________________________
void Storage::TakeFractionFromStock(int IVId,double fraction)
{
DBGL
	if(GetParc())
	{
		if(GetParc()->GetStockManagement() )
		{
			if(fraction > 1 || fraction < 0)
			{
				WARNING << " You try to remove fraction superior than 1 or a negative one..." << endl;
			}
			else
			{
				AddCumulativeIVOut(fIVArray[IVId]*fraction);
				fIVArray[IVId]  -=  fIVArray[IVId] * fraction;
			}

		}
		else
		{
			ERROR << " TakeFractionFromStock can't be DEFINE without REAL stock management" << endl;
			exit(1);
			
		}
	}
	else
	{
		if(fraction > 1 || fraction < 0)
		{
			WARNING << " You try to remove fraction superior than 1 or a negative one..." << endl;
		}
		else
		{
			AddCumulativeIVOut(fIVArray[IVId]*fraction);
			fIVArray[IVId]  -=  fIVArray[IVId] * fraction;
		}

	}
	UpdateInsideIV();	
	DBGL
}

void Storage::TakeFromStock(IsotopicVector isotopicvector)
{

	if(GetParc())
	{
		if(!GetParc()->GetStockManagement())
		{

			AddCumulativeIVOut(isotopicvector);
			fInsideIV -=  isotopicvector;
		}
		else
		{
			ERROR << " TakeFromStock can't be DEFINE WITH REAL stock management" << endl;
			exit(1);
		}
	}
	else
	{
		AddCumulativeIVOut(isotopicvector);
		fInsideIV -=  isotopicvector;

	}

}

//________________________________________________________________________
void Storage::StorageEvolution(cSecond t)
{
DBGL

	if(t ==  fInternalTime && t != 0 ) return;

	RemoveEmptyStocks();

	cSecond EvolutionTime = t - fInternalTime;

#pragma omp parallel for
	for (int i = 0; i <(int) fIVArray.size() ; i++)
	{
		fIVArray[i] = GetDecay(fIVArray[i] , EvolutionTime);
	}
DBGL
}

//________________________________________________________________________
void Storage::Evolution(cSecond t)
{

	// Check if the Storage has been created ...
	if(t ==  fInternalTime && t != 0) return;
	// Make the evolution for the Storage ...

	StorageEvolution(t);
	// Update Inside IV;
	UpdateInsideIV();

	// ... And Finaly update the AbsoluteInternalTime
	fInternalTime = t;
	
}

void Storage::Write(string filename, cSecond date)
{

	for(int i = 0;i < (int)fIVArray.size(); i++)
	{
		
		fIVArray[i].Write(filename, date);
	}

}
//________________________________________________________________________
void Storage::RemoveEmptyStocks()
{
	for(int i = (int)fIVArray.size()-1 ; i >= 0; i--) //Removing empty Stock
		if(fIVArray[i].GetSumOfAll() ==  0)
		{
			fIVArray.erase(fIVArray.begin()+i);
			fIVArrayArrivalTime.erase(fIVArrayArrivalTime.begin()+i);
            
		}
}
