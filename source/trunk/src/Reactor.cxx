#include "CLASSHeaders.hxx"
#include "Reactor.hxx"
#include <omp.h>
//________________________________________________________________________
//
//		Reactor
//
//
//
//
//________________________________________________________________________
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

//________________________________________________________________________
Reactor::~Reactor()
{
	DBGL;
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
		fInternalTime = fCreationTime;
	}

	// Check if the Reactor if started ...
	if(fIsStarted == false) return; 

	long int EvolutionTime = t - fInternalTime; // Calculation of the evolution time (relativ)

	if(EvolutionTime + fInCycleTime < fCycleTime )			//During Cycle
	{
		fIVReactor = fEvolutionDB->GetIsotopicVectorAt(EvolutionTime + fInCycleTime);	// add the new fuel
		fInternalTime += EvolutionTime;							// Update Internal Time
		fInCycleTime += EvolutionTime;							// Update InCycleTime
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
			IsotopicVector IV;
			fIVReactor = IV;
			fShutDown = false;
		}

	}
	else
	{
		// This is so bad!! You will probably unsynchronize all the reactor....
		cout << "!!Warning!! !!!Reactor!!!";
		cout << " Evolution is too log! This is a Bad way to deal the evolution of the reactor...";
		cout << t/365.4/3600/24 << " :" << endl;
		exit (1);
	}
	DBGL;
}

//________________________________________________________________________
void Reactor::Dump()
{
	DBGL;
	if(fEndOfCycle == true  )
	{
		fEndOfCycle = false;

		fAssociedTreatmentFactory->AddIVGodIncome(fIVInCycle - fParc->BuildIsotopicVector(fIVInCycle));

		if(fIsStarted == true)					// A Cycle has already been done
			fAssociedTreatmentFactory->AddIVCooling(fIVOutCycle);
		else fIsStarted = true;					// Just start the first cycle
		if(fShutDown == true) fIsStarted = false;		// shut down the Reactor
	}

	DBGL;
}


//________________________________________________________________________
void Reactor::Write(string Rbasename)
{
	DBGL;
	fIVReactor.Write(Rbasename, fInternalTime);
	DBGL;
}



