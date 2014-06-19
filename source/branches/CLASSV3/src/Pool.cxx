#include "Pool.hxx"

#include "DecayDataBank.hxx"
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

Pool::Pool():CLASSBackEnd()
{
	fOutBackEndFacility = 0;
	SetFacilityType(8);
	SetName("P_Pool.");
}

	//________________________________________________________________________
Pool::Pool(LogFile* log, cSecond coolingtime):CLASSBackEnd(log, coolingtime)
{

	SetFacilityType(8);

	fCycleTime = (cSecond)coolingtime;
	fIsStarted = false;
	fPutToWaste = true;
	fCoolingLastIndex = 0;

	fOutBackEndFacility = 0;
	SetName("P_Pool.");

	
	cout	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	cout	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Cooling Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	cout	<< "!!WARNING!! !!!Pool!!! All Cooled Fuel goes directly to WASTE after cooling !! " << endl;
	
	GetLog()->fLog	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	GetLog()->fLog	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t The Cooling Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "!!WARNING!! !!!Pool!!! All Cooled Fuel goes directly to WASTE after cooling !! " << endl;


}

//________________________________________________________________________
Pool::Pool(LogFile* log, CLASSBackEnd* storage, cSecond coolingtime):CLASSBackEnd(log, coolingtime)
{

	SetFacilityType(8);

	fOutBackEndFacility = storage;
	SetIsStorageType(false);

	fIsStarted = false;
	fPutToWaste = false;
	fCoolingLastIndex = 0;
	SetName("P_Pool.");

	
	cout	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	cout	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;


	cout	<< "\t The Cooling Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	
	GetLog()->fLog	<< "!!INFO!! !!!Pool!!! A new Pool has been define :" << endl;
	GetLog()->fLog	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t The Cooling Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;


}

//________________________________________________________________________
Pool::~Pool()
{


}

//________________________________________________________________________
//________________________________________________________________________
void Pool::SetIVArray(vector<IsotopicVector> ivarray)
{
	cout << "this method as no effect !!!" << endl;
	cout << "Use SetIVArray(vector<IsotopicVector> ivarray, vector<cSecond>n timearray) unstead!!!!"<<endl;
}


//________________________________________________________________________
void Pool::SetIVArray(vector<IsotopicVector> ivarray, vector<cSecond> timearray)
{
	fIVArray = ivarray;
	fCoolingStartingTime =  timearray;

}
//________________________________________________________________________
//	Add Temporary IV : 
//		Cooling
//		
//________________________________________________________________________
void Pool::AddIV(IsotopicVector IV)
{ 

	fIVArray.push_back(IV);
	fInsideIV += IV;
	fCoolingStartingTime.push_back(fInternalTime);
	fCoolingLastIndex++;
	fCoolingIndex.push_back(fCoolingLastIndex);

	AddCumulativeIVIn(IV);

}


//________________________________________________________________________
void Pool::RemoveIVCooling(int i)		//!< Remove a Cooling IsotopicVector
{
	AddCumulativeIVOut(fIVArray[i]);

	fIVArray.erase(fIVArray.begin()+i);
	fCoolingStartingTime.erase(fCoolingStartingTime.begin()+i);
	fCoolingIndex.erase(fCoolingIndex.begin()+i); 

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

	fInsideIV = IsotopicVector();

#pragma omp parallel for
	for ( int i = 0 ; i < (int)fIVArray.size() ; i++)
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
			fIVArray[i] = GetDecay( fIVArray[i], RemainingCoolingTime);


#pragma omp critical(DeleteCoolingIVPB)
			{fCoolingEndOfCycle.push_back(i);}

		}
		else if ( fCoolingStartingTime[i] != t )
		{
			fIVArray[i] = GetDecay( fIVArray[i] , EvolutionTime);
			fInsideIV += fIVArray[i];
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
	
	if(fInternalTime == 0 && !fIsStarted)
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
		
		if(!fPutToWaste)
			fOutBackEndFacility->AddIV(fIVArray[idx]);
		else
			GetParc()->AddWaste(fIVArray[idx]);

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


