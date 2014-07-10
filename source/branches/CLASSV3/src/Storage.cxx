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
	SetIsStorageType();
	SetName("S_Storage.");

}

Storage::Storage(CLASSLogger* log):CLASSBackEnd(log, -1)
{
	SetIsStorageType();

	SetName("S_Storage.");

	INFO << " A new Storage has been define." << endl;
	
}
//________________________________________________________________________
Storage::Storage(CLASSLogger* log, DecayDataBank* evolutivedb):CLASSBackEnd(log, -1)
{
	SetIsStorageType();

	SetDecayDataBank(evolutivedb);
	
	SetName("S_Storage.");

	INFO << " A new Storage has been define." << endl;


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
	AddToFullStock(isotopicvector);

}

//________________________________________________________________________
void Storage::TakeFractionFromStock(int IVId,double fraction)
{

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

				fInsideIV	-= fIVArray[IVId] * fraction;
				fIVArray[IVId]  -= fIVArray[IVId] * fraction;
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

			fInsideIV	-= fIVArray[IVId] * fraction;
			fIVArray[IVId]  -= fIVArray[IVId] * fraction;
		}

	}
	
	

}

void Storage::TakeFromStock(IsotopicVector isotopicvector)
{

	if(GetParc())
	{
		if(!GetParc()->GetStockManagement())
		{

			AddCumulativeIVOut(isotopicvector);
			fInsideIV -= isotopicvector;
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
		fInsideIV -= isotopicvector;

	}

}

//________________________________________________________________________
void Storage::StorageEvolution(cSecond t)
{
DBGL

	if(t == fInternalTime && t !=0 ) return;

	for(int i = (int)fIVArray.size()-1 ; i >=0; i--) //Removing empty Stock
		if(Norme(fIVArray[i]) == 0)
		{
			fIVArray.erase(fIVArray.begin()+i);
			fIVArrayArrivalTime.erase(fIVArrayArrivalTime.begin()+i);

		}
	

	int EvolutionTime = t - fInternalTime;

	fInsideIV = 	GetDecay(fInsideIV , EvolutionTime);

#pragma omp parallel for
	for (int i=0; i <(int) fIVArray.size() ; i++)
	{
		fIVArray[i] = GetDecay(fIVArray[i] , EvolutionTime);
	}
DBGL
}

//________________________________________________________________________
void Storage::Evolution(cSecond t)
{

	// Check if the Storage has been created ...
	if(t == fInternalTime && t!=0) return;
	// Make the evolution for the Storage ...
	StorageEvolution(t);
	// ... And Finaly update the AbsoluteInternalTime
	fInternalTime = t;
	


}

void Storage::Write(string filename, cSecond date)
{

	for(int i=0;i < (int)fIVArray.size(); i++)
	{
		
		fIVArray[i].Write(filename, date);
	}

}

