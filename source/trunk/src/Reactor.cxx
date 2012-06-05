#include "Reactor.hxx"
#include "IsotopicVector.hxx"
#include "EvolutiveProduct.hxx"
#include "CLASS.hxx"
#include "LogFile.hxx"
#include "TreatmentFactory.hxx"
#include "Defines.hxx"

#include <iostream>
#include <algorithm>
#include <omp.h>

//________________________________________________________________________
//
//		Reactor
//
//
//
//
//________________________________________________________________________
ClassImp(Reactor)

Reactor::Reactor()
{
	DBGL;
	DBGL;
}

Reactor::Reactor(EvolutiveProduct* evolutivedb, 
 		 TreatmentFactory* TreatmentFactory,
 		 long int creationtime,
 		 long int lifetime )
{
	DBGL;
	fCycleTime = (long int)(3600*24*365.4)*5;
	fIsStarted = false;
	fShutDown = false;
	fEndOfCycle = false;
	fEvolutionDB = evolutivedb; 
	fAssociedTreatmentFactory = TreatmentFactory;

	fInternalTime = 0;
	fInCycleTime = 0;
	fCreationTime = creationtime;
	fLifeTime = lifetime;

	fIVBeginCycle = fEvolutionDB->GetIsotopicVectorAt(0);
	fIVInCycle = fEvolutionDB->GetIsotopicVectorAt(0);
	fIVOutCycle = fEvolutionDB->GetIsotopicVectorAt(fCycleTime);
	DBGL;
}
Reactor::Reactor(EvolutiveProduct* evolutivedb, 
 		 TreatmentFactory* TreatmentFactory,
 		 long int creationtime,
 		 long int lifetime,
 		 long int cycletime)
{
	DBGL;
	fCycleTime = cycletime;
	fIsStarted = false;
	fShutDown = false;
	fEndOfCycle = false;
	fEvolutionDB = evolutivedb; 
	fAssociedTreatmentFactory = TreatmentFactory;

	fInternalTime = 0;
	fInCycleTime = 0;
	fCreationTime = creationtime;
	fLifeTime = lifetime;

	fIVBeginCycle = fEvolutionDB->GetIsotopicVectorAt(0);
	fIVInCycle = fEvolutionDB->GetIsotopicVectorAt(0);
	fIVOutCycle = fEvolutionDB->GetIsotopicVectorAt(fCycleTime);
	DBGL;
}

Reactor::Reactor(EvolutiveProduct* evolutivedb, 
 		 TreatmentFactory* TreatmentFactory,
 		 long int creationtime,
 		 long int lifetime,
 		 long int cycletime,
 		 IsotopicVector IVincyle,
 		 IsotopicVector IVoutcycle )
{
	DBGL;
	fCycleTime = cycletime;
	fIsStarted = false;
	fShutDown = false;
	fEndOfCycle = false;
	fEvolutionDB = evolutivedb; 
	fAssociedTreatmentFactory = TreatmentFactory;

	fInternalTime = 0;
	fInCycleTime = 0;
	fCreationTime = creationtime;
	fLifeTime = lifetime;

	fIVBeginCycle = fEvolutionDB->GetIsotopicVectorAt(0);
	fIVInCycle = IVincyle;
	fIVOutCycle = IVoutcycle;
	DBGL;
}
//________________________________________________________________________
Reactor::~Reactor()
{
	DBGL;
	DBGL;
}
	void Reactor::SetCycleTime(long int cycletime)
{
	fCycleTime = cycletime;
	fIVOutCycle = fEvolutionDB->GetIsotopicVectorAt(fCycleTime);
}

void Reactor::SetEvolutionDB(EvolutiveProduct* evolutionDB)
{
	DBGL;
	fEvolutionDB = evolutionDB;
	fIVOutCycle = fEvolutionDB->GetIsotopicVectorAt(fCycleTime);
	DBGL;
}



//________________________________________________________________________
void Reactor::Evolution(long int t)
{
	DBGL;
	if(fShutDown == true) return; // Reactor stop...
	
	// Check if the Reactor has been created ...
	if(t<fCreationTime) return;
	DBGL;
	
	
	if(fInternalTime == 0 && fIsStarted == false) // Start of the Reactor
	{
		fEndOfCycle = true;
		fIVReactor  = fIVBeginCycle;
		fInternalTime = t;
	}

	// Check if the Reactor if started ...
	if(fIsStarted == false) return; 

	long int EvolutionTime = t - fInternalTime; // Calculation of the evolution time (relativ)

	if(EvolutionTime + fInCycleTime < fCycleTime )			//During Cycle
	{

		fIVReactor = fEvolutionDB->GetIsotopicVectorAt(EvolutionTime + fInCycleTime);	// update the fuel composition
		fInternalTime += EvolutionTime;							// Update Internal Time
		fInCycleTime += EvolutionTime;							// Update InCycleTime
		if(t == fCreationTime + fLifeTime)	fShutDown = true;
		else
		{
		}
	}
	else if(EvolutionTime + fInCycleTime == fCycleTime )		//End of Cycle
	{
		fEndOfCycle = true;
		fInternalTime += EvolutionTime; 						// Update Internal Time
		
		if(t < fCreationTime + fLifeTime)			//if the Next Cycle Exist...
		{
			fIVReactor  = fIVBeginCycle;
			fInCycleTime = 0;
		}
		else 							//if Not....
		{
			fIVReactor.Clear();
			fShutDown = true;
		}

	}
	else
	{
		// This is so bad!! You will probably unsynchronize all the reactor....
		cout << "!!Warning!! !!!Reactor!!!"
		     << " Evolution is too long! This is a Bad way to deal the evolution of the reactor..."
		     << t/365.4/3600/24 << " :" << endl;
				
		fLog->fLog << "!!Warning!! !!!Reactor!!!"
		           << " Evolution is too long! This is a Bad way to deal the evolution of the reactor..."
		           << t/365.4/3600/24 << " :" << endl;
		exit(1);
	}
	DBGL;
}

//________________________________________________________________________
void Reactor::Dump()
{
DBGL;

	if(fShutDown == true && fIsStarted == false) return; // Reactor stopped...
	
	if(fEndOfCycle == true && fShutDown == false )
	{

		
		fEndOfCycle = false;
		if(fParc->GetStockManagement() == true)
		{
			IsotopicVector BuildIVtmp  = fParc->BuildIsotopicVector(fIVInCycle);
			fAssociedTreatmentFactory->AddIVGodIncome(fIVInCycle - BuildIVtmp);

		}
		else
		{

			IsotopicVector BuildIVtmp ;
			IsotopicVector GodPart;
			BuildIVtmp.Add(fAssociedTreatmentFactory->GetIVFullStock().GetIsotopicQuantity());
			BuildIVtmp -= fIVInCycle;
			GodPart.Add(BuildIVtmp.GetIsotopicQuantityNeeded()) ;
			fAssociedTreatmentFactory->TakeFromStock( fIVInCycle - GodPart, -1);
			fAssociedTreatmentFactory->AddIVGodIncome(GodPart);

		}

		if(fIsStarted == true )					// A Cycle has already been done
			fAssociedTreatmentFactory->AddIVCooling(fIVOutCycle);
		else fIsStarted = true;					// Just start the first cycle

	}
	else if (fEndOfCycle == true && fShutDown == true)		//shutdown at end of Cycle
	{

		fAssociedTreatmentFactory->AddIVCooling(fIVOutCycle);
		fIsStarted = false;		// shut down the Reactor
	}
	else if (fEndOfCycle == false && fShutDown == true) 					//shutdown during Cycle
	{
		fAssociedTreatmentFactory->AddIVCooling(fIVReactor);
		fIVReactor.Clear();
		fIsStarted = false;		// shut down the Reactor
	}
DBGL;
}




