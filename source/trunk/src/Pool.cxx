#include "Pool.hxx"

#include "DataBank.hxx"
#include "IsotopicVector.hxx"
#include "Storage.hxx"
#include "CLASS.hxx"
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

Pool::Pool()
{



}

Pool::Pool(LogFile* log)
{
	
	SetLog(log);
	fCycleTime = 5*3600.*24.*365.25;
	
	fIsStarted = false;
	fPutToWaste = true;
	fInternalTime = 0 ;
	SetCreationTime( 0 );
	fCoolingLastIndex = 0;
	
	cout	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	cout	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl << endl;
	cout	<< "\t The Cooling Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	cout	<< "!!WARNING!! !!!Pool!!! All Cooled Fuel goes directly to WASTE after cooling !! " << endl;
	
	GetLog()->fLog	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	GetLog()->fLog	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t Life time (Operating's Duration) set at \t" << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl << endl;
	GetLog()->fLog	<< "\t The Cooling Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "!!WARNING!! !!!Pool!!! All Cooled Fuel goes directly to WASTE after cooling !! " << endl;
	
	
}
	//________________________________________________________________________
Pool::Pool(LogFile* log, double creation,
				   double coolingtime)
{

	SetLog(log);
	fCycleTime = (cSecond)coolingtime;
	fInternalTime = 0;
	SetCreationTime( (cSecond)creation );
	fIsStarted = false;
	fPutToWaste = true;
	fCoolingLastIndex = 0;
	
	cout	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	cout	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl << endl;
	cout	<< "\t The Cooling Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	cout	<< "!!WARNING!! !!!Pool!!! All Cooled Fuel goes directly to WASTE after cooling !! " << endl;
	
	GetLog()->fLog	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	GetLog()->fLog	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl << endl;
	GetLog()->fLog	<< "\t The Cooling Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "!!WARNING!! !!!Pool!!! All Cooled Fuel goes directly to WASTE after cooling !! " << endl;


}

//________________________________________________________________________
Pool::Pool(LogFile* log, Storage* storage,
				   double creation,
				   double coolingtime)
{

	SetLog(log);
	fCycleTime = (cSecond)coolingtime;
	fInternalTime = 0;
	fStorage = storage;
	SetCreationTime( (cSecond)creation );
	fIsStarted = false;
	fPutToWaste = false;
	fCoolingLastIndex = 0;
	
	
	cout	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	cout	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl << endl;
	cout	<< "\t The Cooling Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	
	GetLog()->fLog	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	GetLog()->fLog	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl << endl;
	GetLog()->fLog	<< "\t The Cooling Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;


}

//________________________________________________________________________
Pool::~Pool()
{


}




//________________________________________________________________________
//	Get Decay
//________________________________________________________________________
IsotopicVector Pool::GetDecay(IsotopicVector isotopicvector, cSecond t)
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
//	Add Temporary IV : 
//		Cooling
//		Separation
//________________________________________________________________________
void Pool::AddIVCooling(IsotopicVector IV)
{ 

	fIVCooling.push_back(IV);
	fCoolingStartingTime.push_back(fInternalTime);
	fCoolingLastIndex++;
	fCoolingIndex.push_back(fCoolingLastIndex);
}


//________________________________________________________________________
void Pool::RemoveIVCooling(int i)		//!< Remove a Cooling IsotopicVector
{

	fIVCooling.erase(fIVCooling.begin()+i);
	fCoolingStartingTime.erase(fCoolingStartingTime.begin()+i);
	fCoolingIndex.erase(fCoolingIndex.begin()+i); 

}

IsotopicVector Pool::GetFullCooling()
{

	IsotopicVector tmp = 0*ZAI(0,0,0);
	
	for(int i =0; i< (int)fIVCooling.size(); i++)
		tmp += fIVCooling[i];
	fInsideIV = tmp;
	return fInsideIV;

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

	if(t == fInternalTime && t!=0) return;
	int RemainingCoolingTime;
	cSecond EvolutionTime = t - fInternalTime;
#pragma omp parallel for
	for ( int i = 0 ; i < (int)fIVCooling.size() ; i++)
	{
		if ( abs(t - fCoolingStartingTime[i] - fCycleTime) < 3600 ) // ">" should not append, only "=" is normal...
		{
			if (t - fCoolingStartingTime[i] > fCycleTime) // Warning & Quit
			{
				cout		<< "!!Warning!! !!!TreamtmentFactory!!! Cooling Step : " << t/3600./24/365.25<< " :"
						<< " An evolution Step is probably missing ! " << " " << endl;
				cout << t << " " << fCoolingStartingTime[i] << " " << fCycleTime << endl;
				GetLog()->fLog 	<< "!!Warning!! !!!TreamtmentFactory!!! Cooling Step : "<< t << " :"
						<< " An evolution Step is probably missing ! " << endl;
				exit (1);
			}
   
			RemainingCoolingTime = fCycleTime - (fInternalTime - fCoolingStartingTime[i]);
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

}


//________________________________________________________________________
void Pool::Evolution(cSecond t)
{

	// Check if the Pool has been created ...
	if(t<GetCreationTime()) return;
	if(t == fInternalTime && t!=0) return;
	
	if(fInternalTime == 0 && fIsStarted == false)
	{
		fInternalTime = GetCreationTime();
		fIsStarted = true;
	}

	// Make the evolution for the Cooling IV ...
	CoolingEvolution(t);
	
	// ... And Finaly update the AbsoluteInternalTime
	fInternalTime = t;
	


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
			GetParc()->AddWaste(fIVCooling[idx]);
		
		fCoolingEndOfCycle.erase(fCoolingEndOfCycle.begin()+i);	// Remove index entry
		RemoveIVCooling(idx);					// Remove IVcooling

	}
	
	if((int)fCoolingEndOfCycle.size() != 0 )// Control
	{
		cout		<< "Problem while Dumping Cooling"<< endl;
		GetLog()->fLog 	<< "Problem while Dumping Cooling"<< endl;
		exit (1);
	}
}


