#include "Storage.hxx"

#include "DataBank.hxx"
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

Storage::Storage()
{
	fDecayDataBase = 0;

}

Storage::Storage(LogFile* log)
{
	
	SetLog(log);
	fDecayDataBase = 0;
	cout	<< "!!INFO!! !!!Storage!!! A new Storage has been define." << endl;
	
	GetLog()->fLog	<< "!!INFO!! !!!Storage!!! A new Storage has been define." << endl;
	
}
//________________________________________________________________________
Storage::Storage(LogFile* log, DataBank<ZAI>* evolutivedb)
{

	SetLog(log);
	fDecayDataBase = evolutivedb;
	
	cout	<< "!!INFO!! !!!Storage!!! A new Storage has been define." << endl;
	
	GetLog()->fLog	<< "!!INFO!! !!!Storage!!! A new Storage has been define." << endl;


}

//________________________________________________________________________
Storage::~Storage()
{


}


//________________________________________________________________________
//	Get Decay
//________________________________________________________________________
IsotopicVector Storage::GetDecay(IsotopicVector isotopicvector, cSecond t)
{

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

	return IV;

}

//________________________________________________________________________
void Storage::ClearStock()
{

	IsotopicVector EmptyIV;
	fInsideIV = EmptyIV;
	fIVStock.clear();

}

//________________________________________________________________________
void Storage::AddToStock(IsotopicVector isotopicvector)
{

	AddCumulativeIVIn(isotopicvector);

	if(GetParc()->GetStockManagement() == true)
		fIVStock.push_back(isotopicvector);
	AddToFullStock(isotopicvector);

}

//________________________________________________________________________
void Storage::TakeFractionFromStock(int IVId,double fraction)
{

	if(GetParc()->GetStockManagement() == true)
	{
		if(fraction > 1 || fraction < 0)
		{
			cout << "!!Warning!! !!!Storage!!! You try to remove fraction superior than 1 or a negative one..." << endl;
			GetLog()->fLog << "!!Warning!! !!!Storage!!! You try to remove fraction superior than 1 or a negative one..." << endl;
		}
		else 
		{
			fInsideIV -= fIVStock[IVId]*fraction;
			fIVStock[IVId] = fIVStock[IVId]*(1-fraction);

			AddCumulativeIVOut(fIVStock[IVId]*fraction);
			
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

	for(int i = (int)fIVStock.size()-1 ; i >=0; i--) //Removing empty Stock
		if(Norme(fIVStock[i]) == 0)
			fIVStock.erase(fIVStock.begin()+i); 
	
	

	int EvolutionTime = t - fInternalTime;

	fInsideIV = 	GetDecay(fInsideIV , EvolutionTime);

#pragma omp parallel for
	for (int i=0; i <(int) fIVStock.size() ; i++)
	{
		fIVStock[i] = GetDecay(fIVStock[i] , EvolutionTime);
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

	for(int i=0;i < (int)fIVStock.size(); i++)
	{
		
		fIVStock[i].Write(filename, date);
	}

}

