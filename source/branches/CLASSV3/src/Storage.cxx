#include "Storage.hxx"

#include "DecayDataBank.hxx"
#include "CLASS.hxx"
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

Storage::Storage():CLASSBackEnd()
{
	SetFacilityType(-1);
	SetIsStorageType();
	SetName("S_Storage.");

}

Storage::Storage(LogFile* log):CLASSBackEnd(log)
{
	SetFacilityType(-1);
	SetIsStorageType();

	SetName("S_Storage.");

	cout	<< "!!INFO!! !!!Storage!!! A new Storage has been define." << endl;
	
	GetLog()->fLog	<< "!!INFO!! !!!Storage!!! A new Storage has been define." << endl;
	
}
//________________________________________________________________________
Storage::Storage(LogFile* log, DecayDataBank* evolutivedb):CLASSBackEnd(log)
{
	SetFacilityType(-1);
	SetIsStorageType();

	SetDecayDataBank(evolutivedb);
	
	SetName("S_Storage.");

	cout	<< "!!INFO!! !!!Storage!!! A new Storage has been define." << endl;
	
	GetLog()->fLog	<< "!!INFO!! !!!Storage!!! A new Storage has been define." << endl;


}

//________________________________________________________________________
Storage::~Storage()
{


}

//________________________________________________________________________
void Storage::AddIV(IsotopicVector isotopicvector)
{

	AddCumulativeIVIn(isotopicvector);

	if(GetParc()->GetStockManagement() )
		fIVArray.push_back(isotopicvector);
	AddToFullStock(isotopicvector);

}

//________________________________________________________________________
void Storage::TakeFractionFromStock(int IVId,double fraction)
{

	if(GetParc()->GetStockManagement() )
	{
		if(fraction > 1 || fraction < 0)
		{
			cout << "!!Warning!! !!!Storage!!! You try to remove fraction superior than 1 or a negative one..." << endl;
			GetLog()->fLog << "!!Warning!! !!!Storage!!! You try to remove fraction superior than 1 or a negative one..." << endl;
		}
		else 
		{
			AddCumulativeIVOut(fIVArray[IVId]*fraction);

			fInsideIV -= fIVArray[IVId]*fraction;
			fIVArray[IVId] = fIVArray[IVId]*(1-fraction);
		}

	}
	else
	{
		cout << "!!Warning!! !!!Storage!!! TakeFractionFromStock can't be DEFINE without REAL stock management" << endl;
		GetLog()->fLog << "!!Warning!! !!!Storage!!! TakeFractionFromStock can't be DEFINE without REAL stock management" << endl;
		exit(1);

	}
	
	

}

void Storage::TakeFromStock(IsotopicVector isotopicvector)
{

	if(GetParc()->GetStockManagement() == false)
	{

		AddCumulativeIVOut(isotopicvector);
		fInsideIV -= isotopicvector;
	}
	else
	{
		cout << "!!Warning!! !!!Storage!!! TakeFromStock can't be DEFINE WITH REAL stock management" << endl;
		GetLog()->fLog << "!!Warning!! !!!Storage!!! TakeFromStock can't be DEFINE WITH REAL stock management" << endl;
		exit(1);
	}

}

//________________________________________________________________________
void Storage::StorageEvolution(cSecond t)
{


	if(t == fInternalTime && t !=0 ) return;

	for(int i = (int)fIVArray.size()-1 ; i >=0; i--) //Removing empty Stock
		if(Norme(fIVArray[i]) == 0)
			fIVArray.erase(fIVArray.begin()+i); 
	
	

	int EvolutionTime = t - fInternalTime;

	fInsideIV = 	GetDecay(fInsideIV , EvolutionTime);

#pragma omp parallel for
	for (int i=0; i <(int) fIVArray.size() ; i++)
	{
		fIVArray[i] = GetDecay(fIVArray[i] , EvolutionTime);
	}
	

	

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

