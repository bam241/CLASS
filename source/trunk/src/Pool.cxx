#include "Pool.hxx"

#include "DataBank.hxx"
#include "IsotopicVector.hxx"
#include "Storage.hxx"
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
//		Pool
//
//
//
//
//________________________________________________________________________
ClassImp(Pool)

Pool::Pool(LogFile* log)
{
DBGL;
	fLog = log;
	fCoolingTime = 5*3600.*24.*365.25;

	IsStarted = false;
	fPutToWaste = true;
	fInternalTime = 0;
	fCreationTime = 0;
	fCoolingLastIndex = 0;

	cout	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	cout	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl << endl;
	cout	<< "\t The Cooling Time set at\t " << (double)(fCoolingTime/3600/24/365.25) << " year" << endl;
	cout	<< "!!WARNING!! !!!Pool!!! All Cooled Fuel goes directly to WASTE after cooling !! " << endl;

	fLog->fLog	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	fLog->fLog	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "\t Life time (Operating's Duration) set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl << endl;
	fLog->fLog	<< "\t The Cooling Time set at\t " << (double)(fCoolingTime/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "!!WARNING!! !!!Pool!!! All Cooled Fuel goes directly to WASTE after cooling !! " << endl;

DBGL;
}
	//________________________________________________________________________
Pool::Pool(LogFile* log, double creation,
				   double coolingtime)
{
DBGL;
	fLog = log;
	fCoolingTime = (cSecond)coolingtime;
	fInternalTime = 0;
	fCreationTime = (cSecond)creation;
	IsStarted = false;
	fPutToWaste = true;
	fCoolingLastIndex = 0;
	
	cout	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	cout	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl << endl;
	cout	<< "\t The Cooling Time set at\t " << (double)(fCoolingTime/3600/24/365.25) << " year" << endl;
	cout	<< "!!WARNING!! !!!Pool!!! All Cooled Fuel goes directly to WASTE after cooling !! " << endl;
	
	fLog->fLog	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	fLog->fLog	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "\t Life time (Operating's Duration) set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl << endl;
	fLog->fLog	<< "\t The Cooling Time set at\t " << (double)(fCoolingTime/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "!!WARNING!! !!!Pool!!! All Cooled Fuel goes directly to WASTE after cooling !! " << endl;

DBGL;
}

//________________________________________________________________________
Pool::Pool(LogFile* log, Storage* storage,
				   double creation,
				   double coolingtime)
{
DBGL;
	fLog = log;
	fCoolingTime = (cSecond)coolingtime;
	fInternalTime = 0;
	fStorage = storage;
	fCreationTime = (cSecond)creation;
	IsStarted = false;
	fPutToWaste = false;
	fCoolingLastIndex = 0;
	
	
	cout	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	cout	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl << endl;
	cout	<< "\t The Cooling Time set at\t " << (double)(fCoolingTime/3600/24/365.25) << " year" << endl;
	
	fLog->fLog	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	fLog->fLog	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "\t Life time (Operating's Duration) set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl << endl;
	fLog->fLog	<< "\t The Cooling Time set at\t " << (double)(fCoolingTime/3600/24/365.25) << " year" << endl;

DBGL;
}

//________________________________________________________________________
Pool::~Pool()
{
DBGL;
DBGL;
}




//________________________________________________________________________
//	Get Decay
//________________________________________________________________________
IsotopicVector Pool::GetDecay(IsotopicVector isotopicvector, cSecond t)
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
//	Add Temporary IV : 
//		Cooling
//		Separation
//________________________________________________________________________
void Pool::AddIVCooling(IsotopicVector IV)
{ 
DBGL;
	fIVCooling.push_back(IV);
	fCoolingStartingTime.push_back(fInternalTime);
	fCoolingLastIndex++;
	fCoolingIndex.push_back(fCoolingLastIndex);
DBGL;
}


//________________________________________________________________________
void Pool::RemoveIVCooling(int i)		//!< Remove a Cooling IsotopicVector
{
DBGL;
	fIVCooling.erase(fIVCooling.begin()+i);
	fCoolingStartingTime.erase(fCoolingStartingTime.begin()+i);
	fCoolingIndex.erase(fCoolingIndex.begin()+i); 
DBGL;
}

IsotopicVector Pool::GetFullCooling()
{
DBGL;
	IsotopicVector tmp = 0*ZAI(0,0,0);
	
	for(int i =0; i< (int)fIVCooling.size(); i++)
		tmp += fIVCooling[i];
	
	return tmp;
DBGL;
}



//________________________________________________________________________
//	Time Action with the reste of the Universe : 
//		Out Storage
//		Evolution : 
//			Cooling
//________________________________________________________________________




//________________________________________________________________________
void Pool::CoolingEvolution(cSecond t)
{
DBGL;
	if(t == fInternalTime && t!=0) return;
	int RemainingCoolingTime;
	cSecond EvolutionTime = t - fInternalTime;
#pragma omp parallel for
	for ( int i = 0 ; i < (int)fIVCooling.size() ; i++)
	{
		if ( abs(t - fCoolingStartingTime[i] - fCoolingTime) < 3600 ) // ">" should not append, only "=" is normal...
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
			fIVCooling[i] = GetDecay( fIVCooling[i], RemainingCoolingTime);

#pragma omp critical(DeleteCoolingIVPB)
			{fCoolingEndOfCycle.push_back(i);}

		}
		else if ( fCoolingStartingTime[i] != t )
		{
			fIVCooling[i] = GetDecay( fIVCooling[i] , EvolutionTime);
		}
	}
#pragma omp critical(DeleteCoolingIVPB)
	{sort (fCoolingEndOfCycle.begin(), fCoolingEndOfCycle.end());}
DBGL;
}


//________________________________________________________________________
void Pool::Evolution(cSecond t)
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

	// Make the evolution for the Cooling IV ...
	CoolingEvolution(t);
	
	// ... And Finaly update the AbsoluteInternalTime
	fInternalTime = t;
	

DBGL;
}


//________________________________________________________________________
void Pool::Dump()
{
//------ Cooling ------//
	for(int i = (int)fCoolingEndOfCycle.size()-1; i >=0 ; i--)	// IV End Of Cooling
	{
	
		int idx = fCoolingEndOfCycle[i];			// Get Index number
		
		if(fPutToWaste == false)
			fStorage->AddToStock(fIVCooling[idx]);
		else
			fParc->AddWaste(fIVCooling[idx]);
		
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


