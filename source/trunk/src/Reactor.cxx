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
 		 double creationtime, double lifetime)
{
	DBGL;
	
	fIsStarted = false;
	fShutDown = false;
	fEndOfCycle = false;
	
	fFabricationPlant = fabricationplant;
	fFixedFuel = false;
	fBurnUp = -1.;
	fHeavyMetalMass = -1.;
	
	fAssociedTreatmentFactory = treatmentfactory;
	
	fFuelTypeDB = fueltypeDB;
	
	fInternalTime = 0;
	fInCycleTime = 0;
	fPower = -1.;
	fCycleTime = -1.;	 //BU in GWd/t
	
	fCreationTime = (cSecond)creationtime;
	fLifeTime = (cSecond)lifetime;
	
	cout	<< "!!Info!! !!!Reactor!!! A Reactor has been define :" << endl;
	cout	<< "\t Fuel Composition is not fixed ! "<< endl;
	cout	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	cout	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl << endl;
	cout	<< "!!WARNING!! !!!Reactor!!! You need to set Burn-up/Power/CycleTime (2 of 3) & Heavy Metal Mass Manualy !! " << endl;

	DBGL;
}

Reactor::Reactor(double Power, EvolutionDataBase<IsotopicVector>* 	fueltypeDB,
		 FabricationPlant* fabricationplant,
 		 TreatmentFactory* treatmentfactory,
 		 double creationtime, double lifetime,
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
	fPower = Power;
	fCycleTime = (cSecond) (BurnUp*1e9 / (fPower)  * HMMass  *3600*24);	 //BU in GWd/t

	fCreationTime = (cSecond)creationtime;
	fLifeTime = (cSecond)lifetime;
	
	
	
	cout	<< "!!Info!! !!!Reactor!!! A Reactor has been define :" << endl;
	cout	<< "\t Fuel Composition is not fixed ! "<< endl;
	cout	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	cout	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Effective Thermal Power set at \t " << (double)(fPower *1e-6) << " MW" << endl;
	cout	<< "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	cout	<< "\t The corresponding Cycle Time is\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;
	
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
	fCycleTime = (cSecond)cycletime;
	fCreationTime = (cSecond)creationtime;
	fLifeTime = (cSecond)lifetime;
	fPower = BurnUp*3600.*24. / (fCycleTime) * HMMass *1e9; //BU in GWd/t
	
	
	cout	<< "!!Info!! !!!Reactor!!! A Reactor has been define :" << endl;
	cout	<< "\t Fuel Composition is not fixed ! "<< endl;
	cout	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	cout	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Cycle Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	cout	<< "\t The corresponding Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW" << endl;
	cout	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;
	
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
	fCycleTime = (cSecond)cycletime;
	fCreationTime = (cSecond)creationtime;
	fLifeTime = (cSecond)lifetime;

	fPower = fEvolutionDB.GetPower();

	fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVInCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fCycleTime/fEvolutionDB.GetPower()*fPower) );
	
	
	cout	<< "!!Info!! !!!Reactor!!! A Reactor has been define :" << endl;
	cout	<< "\t Fuel Composition is fixed ! "<< endl;
	cout	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	cout	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Cycle Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW" << endl;
	cout	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;
	
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
	if(fFixedFuel==true)
	{
		fCycleTime = (cSecond)cycletime;
		fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt(fCycleTime/fEvolutionDB.GetPower()*fPower);
	}
	else
	{
		fCycleTime = (cSecond)cycletime;
		fPower = fBurnUp*3600*24 / (fCycleTime) * fHeavyMetalMass *1e9; //BU in GWd/t
	}
}
	//________________________________________________________________________
void Reactor::SetPower(double Power)
{
	if(fFixedFuel==true)
	{
		fPower = Power;
		fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fCycleTime/fEvolutionDB.GetPower()*fPower) );
	}
	else
	{
		fPower = Power;
		fCycleTime = (cSecond)(fBurnUp*1e9 / (fPower)  * fHeavyMetalMass  *3600*24);	 //BU in GWd/t
	}

}

//________________________________________________________________________
void Reactor::SetEvolutionDB(EvolutiveProduct evolutionDB)
{
	DBGL;
	fEvolutionDB = evolutionDB;
	fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fCycleTime/fEvolutionDB.GetPower()*fPower) );
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
void Reactor::Evolution(cSecond t)
{
	DBGL;

	if( fShutDown == true || t < fCreationTime ) return; // Reactor stop or not started...

#pragma omp critical(ParcPowerUpdate)
	{fParc->AddToPower(fPower);}
	
	
	if( t == fInternalTime && t!=0 ) return
	DBGL;
	

	if(fInternalTime == 0 && fIsStarted == false) // Start of the Reactor
	{
		fEndOfCycle = true;
		fIVReactor  = fIVBeginCycle;
		fInternalTime = t;
		
	}

	// Check if the Reactor if started ...
	if(fIsStarted == false) return;			// If the reactor just start don't need to make Fuel evolution 


	cSecond EvolutionTime = t - fInternalTime; // Calculation of the evolution time (relativ)

	if( abs(EvolutionTime + fInCycleTime - fCycleTime) < 3600 )		//End of Cycle
	{
		fEndOfCycle = true;
		fInternalTime += EvolutionTime; 				// Update Internal Time
		fInCycleTime += EvolutionTime;					// Update InCycleTime
		
		if(t >= fCreationTime + fLifeTime)				// if the Next Cycle don't 'Exist...
			fShutDown = true;
	
	}
	else if(EvolutionTime + fInCycleTime < fCycleTime )			// During Cycle
	{
		
		fInternalTime += EvolutionTime;					// Update Internal Time
		fInCycleTime += EvolutionTime;					// Update InCycleTime
	
		fIVReactor = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fInCycleTime/fEvolutionDB.GetPower()*fPower) );	// update the fuel composition
		
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
				fParc->AddGod(GodPart);

			}
			else	fParc->AddGod(fIVInCycle);
			
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




