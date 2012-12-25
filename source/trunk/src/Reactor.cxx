#include "Reactor.hxx"

#include "EvolutiveProduct.hxx"
#include "EvolutionDataBase.hxx"
#include "TreatmentFactory.hxx"
#include "FabricationPlant.hxx"
#include "Storage.hxx"
#include "CLASS.hxx"

#include "LogFile.hxx"
#include "Defines.hxx"

#include <iostream>
#include <cmath>
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


Reactor::Reactor(EvolutionDataBase<IsotopicVector>* 	fueltypeDB,
		 FabricationPlant* fabricationplant,
 		 TreatmentFactory* treatmentfactory,
 		 double creationtime, double lifetime, double cycletime,
 		 double HMMass, double BurnUp)
{
	DBGL;

	fIsStarted = false;
	fShutDown = false;
	fEndOfCycle = false;

	fFabricationPlant = fabricationplant;
	fFixedFuel = false;
	fBurnUp = BurnUp;
	fHeavyMetalMass = HMMass;

	fAssociedTreatmentFactory = treatmentfactory;

	fFuelTypeDB = fueltypeDB;

	fInternalTime = 0;
	fInCycleTime = 0;
	fCycleTime = cycletime;
	fCreationTime = creationtime;
	fLifeTime = lifetime;
	fPower = BurnUp / (cycletime/3600/24) *1e9 * HMMass; //BU in GWd/t
	DBGL;
}

Reactor::Reactor(EvolutiveProduct evolutivedb, 
 		 TreatmentFactory* treatmentfactory,
 		 double creationtime,
 		 double lifetime,
 		 double cycletime)
{
	DBGL;

	fIsStarted = false;
	fShutDown = false;
	fEndOfCycle = false;

	fEvolutionDB = evolutivedb; 
	fFixedFuel = true;
	fIsStorage = false;

	fAssociedTreatmentFactory = treatmentfactory;

	fInternalTime = 0;
	fInCycleTime = 0;
	fCycleTime = cycletime;
	fCreationTime = creationtime;
	fLifeTime = lifetime;

	fPower = fEvolutionDB.GetPower();

	fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVInCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt(fCycleTime/fEvolutionDB.GetPower()*fPower);
	DBGL;
}


//________________________________________________________________________
Reactor::~Reactor()
{
	DBGL;
	DBGL;
}

//________________________________________________________________________
void Reactor::SetCycleTime(double cycletime)
{
	fCycleTime = cycletime;
	fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt(fCycleTime/fEvolutionDB.GetPower()*fPower);
}

//________________________________________________________________________
void Reactor::SetEvolutionDB(EvolutiveProduct evolutionDB)
{
	DBGL;
	fEvolutionDB = evolutionDB;
	fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt(fCycleTime/fEvolutionDB.GetPower()*fPower);
	fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);

	DBGL;
}

//________________________________________________________________________
void Reactor::SetNewFuel(EvolutiveProduct ivdb)
{
	DBGL;
	SetEvolutionDB(ivdb);
	DBGL;
}

//________________________________________________________________________
void Reactor::Evolution(double t)
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

	double EvolutionTime = t - fInternalTime; // Calculation of the evolution time (relativ)

	if( abs(EvolutionTime + fInCycleTime - fCycleTime) < 3600 )		//End of Cycle
	{
		fEndOfCycle = true;
		fInternalTime += EvolutionTime; 						// Update Internal Time
		fInCycleTime += EvolutionTime;							// Update InCycleTime
		
		if(t >= fCreationTime + fLifeTime)			//if the Next Cycle don't 'Exist...
			fShutDown = true;
	
	}
	else if(EvolutionTime + fInCycleTime < fCycleTime )			//During Cycle
	{
		
		fInternalTime += EvolutionTime;							// Update Internal Time
		fInCycleTime += EvolutionTime;							// Update InCycleTime
	
		fIVReactor = fEvolutionDB.GetIsotopicVectorAt( (fInCycleTime)/fEvolutionDB.GetPower()*fPower );	// update the fuel composition
		
		if(t>=fCreationTime + fLifeTime)	fShutDown = true;
	} 
	else
	{
		// This is so bad!! You will probably unsynchronize all the reactor....
		cout << "!!Warning!! !!!Reactor!!!"
		     << " Evolution is too long! This is a Bad way to deal the evolution of the reactor..."
		     << t/365.25/3600/24 << " :" << endl;
				
		fLog->fLog << "!!Warning!! !!!Reactor!!!"
		           << " Evolution is too long! This is a Bad way to deal the evolution of the reactor..."
		           << t/365.25/3600/24 << " :" << endl;
		exit(1);
	}

	DBGL;
}

//________________________________________________________________________
void Reactor::Dump()
{
DBGL;

	if(fInternalTime < fCreationTime) return;
	if(fShutDown == true && fIsStarted == false) return; // Reactor stopped...

	if(fFixedFuel == true)
	{

		if(fEndOfCycle == true && fShutDown == false )
		{
			fEndOfCycle = false;

			if(fIsStarted == true )					// A Cycle has already been done
			fAssociedTreatmentFactory->AddIVCooling(fIVOutCycle);
			else fIsStarted = true;					// Just start the first cycle

			if(fParc->GetStockManagement() == false && fIsStorage == true)
			{
				IsotopicVector BuildIVtmp ;
				IsotopicVector GodPart;
				
				//Get The Storage Compostion
				BuildIVtmp.Add(fStorage->GetFullStock().GetIsotopicQuantity());
				//Get the rest after IVIn creation
				BuildIVtmp -= fIVInCycle;
				//Get the God part form this rest
				GodPart.Add(BuildIVtmp.GetIsotopicQuantityNeeded()) ;
				//Take what you can from Storage...
				fStorage->TakeFromStock( fIVInCycle - GodPart);
				//And Get the rest from God
				fParc->AddGodIncome(GodPart);

			}
			else	fParc->AddGodIncome(fIVInCycle);
			
			fIVReactor  = fIVBeginCycle;
			fInCycleTime = 0;
		}
		else if (fEndOfCycle == true && fShutDown == true)		//shutdown at end of Cycle
		{

			fAssociedTreatmentFactory->AddIVCooling(fIVOutCycle);
			fIVReactor.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
		else if (fEndOfCycle == false && fShutDown == true) 					//shutdown during Cycle
		{
			fAssociedTreatmentFactory->AddIVCooling(fIVReactor);
			fIVReactor.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
	}
	else
	{
		if(fParc->GetStockManagement() == false)
		{
			cout << "!!Warning!! !!!Reactor!!! Can't have unfixedFuel without stock management'" << endl;
			fLog->fLog << "!!Warning!! !!!Reactor!!! Can't have unfixedFuel without stock management" << endl;
			exit(1);
		}

		
		if(fEndOfCycle == true && fShutDown == false )
		{
			fEndOfCycle = false;

			if(fIsStarted == true )					// A Cycle has already been done
			{
				fAssociedTreatmentFactory->AddIVCooling(fIVOutCycle);
			}
			else fIsStarted = true;					// Just start the first cycle

			SetNewFuel(fFabricationPlant->GetReactorEvolutionDB(fId));
			fFabricationPlant->TakeReactorFuel(fId);


			fIVReactor  = fIVBeginCycle;
			fInCycleTime = 0;

		}
		else if (fEndOfCycle == true && fShutDown == true)		//shutdown at end of Cycle
		{
			fAssociedTreatmentFactory->AddIVCooling(fIVOutCycle);
			fIVReactor.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
		else if (fEndOfCycle == false && fShutDown == true) 					//shutdown during Cycle
		{
			fAssociedTreatmentFactory->AddIVCooling(fIVReactor);
			fIVReactor.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
		


	}

DBGL;
}




