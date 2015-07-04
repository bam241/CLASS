#include "Pool.hxx"

#include "IsotopicVector.hxx"
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
//		Pool
//
//
//
//
//________________________________________________________________________
ClassImp(Pool)


Pool::Pool():CLASSBackEnd(8)
{
	fOutBackEndFacility = 0;
	SetName("P_Pool.");
}

	//________________________________________________________________________
Pool::Pool(CLASSLogger* log, cSecond coolingtime):CLASSBackEnd(log, coolingtime, 8)
{
	DBGL


	fCycleTime = (cSecond)coolingtime;
	fPutToWaste = true;
	fCoolingLastIndex = 0;

	fOutBackEndFacility = 0;
	SetName("P_Pool.");

	
	INFO	 <<  " A new Pool has been define :" << endl;
	INFO	 <<  "\t The Cooling Time set at\t " << (double)(fCycleTime/cYear) << " year" << endl;
	WARNING	 <<  " All Cooled Fuel goes directly to WASTE after cooling !! " << endl;

	DBGL

}

//________________________________________________________________________
Pool::Pool(CLASSLogger* log, CLASSBackEnd* storage, cSecond coolingtime):CLASSBackEnd(log, coolingtime, 8)
{
	DBGL

	fOutBackEndFacility = storage;
	SetIsStorageType(false);

	fPutToWaste = false;
	fCoolingLastIndex = 0;
	SetName("P_Pool.");

	

	INFO	 <<  " A new Pool has been define :" << endl;
	INFO	 <<  "\t The Cooling Time set at\t " << (double)(fCycleTime/cYear) << " year" << endl;

	DBGL

}

//________________________________________________________________________
Pool::~Pool()
{


}

//________________________________________________________________________
//________________________________________________________________________
void Pool::SetIVArray(vector<IsotopicVector> ivarray)
{
	INFO << "This method as no effect !!!" << endl;
	INFO << "Use SetIVArray(vector<IsotopicVector> ivarray, vector<cSecond>n timearray) unstead!!!!" << endl;
}


//________________________________________________________________________
void Pool::SetIVArray(vector<IsotopicVector> ivarray, vector<cSecond> timearray)
{
	fIVArray = ivarray;
	fIVArrayArrivalTime =  timearray;

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
	fIVArrayArrivalTime.push_back(fInternalTime);
	fCoolingLastIndex++;
	fCoolingIndex.push_back(fCoolingLastIndex);

	AddCumulativeIVIn(IV);

}


//________________________________________________________________________
void Pool::RemoveIVCooling(int i)		//!< Remove a Cooling IsotopicVector
{
	AddCumulativeIVOut(fIVArray[i]);

	fIVArray.erase(fIVArray.begin()+i);
	fIVArrayArrivalTime.erase( fIVArrayArrivalTime.begin()+i);
	fCoolingIndex.erase(fCoolingIndex.begin()+i);
	UpdateInsideIV();
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
DBGL

	if(t == fInternalTime && t!=0) return;
	int RemainingCoolingTime;
	cSecond EvolutionTime = t - fInternalTime;

#pragma omp parallel for
	for ( int i = 0 ; i < (int)fIVArray.size() ; i++)
	{
		if ( abs(t -  fIVArrayArrivalTime[i] - fCycleTime) < 3600 ) // ">" should not append, only "=" is normal...
		{
			if (t -  fIVArrayArrivalTime[i] > fCycleTime) // Warning & Quit
			{
				ERROR << " Cooling Step : " << t/cYear <<  " :" << " an evolution Step is probably missing ! " << " " << endl;
				exit (1);
			}
   
			RemainingCoolingTime = fCycleTime - (fInternalTime -  fIVArrayArrivalTime[i]);
			//Cooling Decay
			fIVArray[i] = GetDecay( fIVArray[i], RemainingCoolingTime);


#pragma omp critical(DeleteCoolingIVPB)
			{fCoolingEndOfCycle.push_back(i);}

		}
		else if (  fIVArrayArrivalTime[i] != t )
		{
			fIVArray[i] = GetDecay( fIVArray[i] , EvolutionTime);
		}
	}
#pragma omp critical(DeleteCoolingIVPB)
	{sort (fCoolingEndOfCycle.begin(), fCoolingEndOfCycle.end());}

DBGL
}


//________________________________________________________________________
void Pool::Evolution(cSecond t)
{

	// Check if the Pool has been created ...
	if(t == fInternalTime && t!=0) return;
	// Make the evolution for the Cooling IV ...
	CoolingEvolution(t);
	// Update Inside IV
	UpdateInsideIV();
	// ... And Finaly update the AbsoluteInternalTime
	fInternalTime = t;
	


}


//________________________________________________________________________
void Pool::Dump()
{
DBGL
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
		ERROR << "Problem while Dumping Cooling" <<  endl;
		exit (1);
	}
DBGL
}


