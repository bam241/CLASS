#include "Reactor.hxx"

#include "EvolutionData.hxx"
#include "PhysicModels.hxx"
#include "Pool.hxx"
#include "FabricationPlant.hxx"
#include "Storage.hxx"
#include "Scenario.hxx"
#include "CLASSConstante.hxx"

#include "CLASSLogger.hxx"

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

Reactor::Reactor():CLASSFacility()
{

	SetFacilityType(4);
	SetName("R_Reactor.");


	fOutBackEndFacility = 0;
	fStorage = 0;
	fFuelTypeDB = 0;
	fFabricationPlant = 0;


	fNextPlan = fLoadingPlan.begin();
}

Reactor::Reactor(CLASSLogger* log):CLASSFacility(log, 4)
{

	fOutBackEndFacility = 0;
	fStorage = 0;
	fFuelTypeDB = 0;
	fFabricationPlant = 0;
	(*this).SetName("R_Reactor.");
	fNextPlan = fLoadingPlan.begin();

}

Reactor::Reactor(CLASSLogger* log, PhysicModels* fueltypeDB,
		 FabricationPlant* fabricationplant,
 		 CLASSBackEnd* Pool,
 		 cSecond creationtime, cSecond lifetime):CLASSFacility(log, creationtime, lifetime, 4)
{

	(*this).SetName("R_Reactor.");


	fIsStarted = false;
	fIsShutDown = false;
	fIsAtEndOfCycle = false;

	fFabricationPlant = fabricationplant;
	fFixedFuel = false;
	fBurnUp = -1.;
	fHeavyMetalMass = -1.;
	fStorage = 0;

	fOutBackEndFacility = Pool;

	fFuelTypeDB = fueltypeDB;

	fPower = -1.;

	fNextPlan = fLoadingPlan.begin();


	INFO << " A Reactor has been define :" << endl;
	INFO << "\t Fuel Composition is not fixed ! "<< endl;
	INFO << "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	INFO << "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl;

	WARNING << " You need to set Burn-up/Power/CycleTime (2 of 3) & Heavy Metal Mass Manualy !! " << endl;

}

Reactor::Reactor(CLASSLogger* log, PhysicModels* fueltypeDB, FabricationPlant* fabricationplant, CLASSBackEnd* Pool,
 		 cSecond creationtime, cSecond lifetime,
 		 double Power, double HMMass, double BurnUp, double ChargeFactor):CLASSFacility(log, creationtime, lifetime, 4)
{
	(*this).SetName("R_Reactor.");


	fStorage = 0;
	fIsStarted = false;
	fIsShutDown = false;
	fIsAtEndOfCycle = false;

	fFabricationPlant = fabricationplant;
	fFixedFuel = false;

	fOutBackEndFacility = Pool;

	fFuelTypeDB = fueltypeDB;


	fBurnUp = BurnUp;
	fHeavyMetalMass = HMMass;
	fPower = Power*ChargeFactor;
	fCycleTime = (cSecond) (fBurnUp*1e9 / (fPower)  * fHeavyMetalMass  *3600*24);	 //BU in GWd/t

	fNextPlan = fLoadingPlan.begin();


	INFO << " A Reactor has been define :" << endl;
	INFO << "\t Fuel Composition is not fixed ! "<< endl;
	INFO << "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	INFO << "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl;
	INFO << "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << Power << " and " << ChargeFactor << " Charge Factor)"<< endl;
	INFO << "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	INFO << "\t The corresponding Cycle Time is\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	INFO << "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;



}

Reactor::Reactor(CLASSLogger* log, PhysicModels* 	fueltypeDB,
		 FabricationPlant* fabricationplant,
 		 CLASSBackEnd* Pool,
 		 cSecond creationtime, cSecond lifetime, cSecond cycletime,
 		 double HMMass, double BurnUp):CLASSFacility(log, creationtime, lifetime, cycletime, 4)
{
	(*this).SetName("R_Reactor.");


	fIsStarted = false;
	fIsShutDown = false;
	fIsAtEndOfCycle = false;

	fStorage = 0;

	fFabricationPlant = fabricationplant;
	fFixedFuel = false;
	fBurnUp = BurnUp;
	fHeavyMetalMass = HMMass;

	fOutBackEndFacility = Pool;

	fFuelTypeDB = fueltypeDB;

	fPower = BurnUp*3600.*24. / (fCycleTime) * HMMass *1e9; //BU in GWd/t

	fNextPlan = fLoadingPlan.begin();

	INFO << " A Reactor has been define :" << endl;
	INFO << "\t Fuel Composition is not fixed ! "<< endl;
	INFO << "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	INFO << "\t Life time (Operating's Duration) set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	INFO << "\t The Cycle Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	INFO << "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	INFO << "\t The corresponding Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW" << endl;
	INFO << "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;



}


Reactor::Reactor(CLASSLogger* log, EvolutionData evolutivedb,
 		 CLASSBackEnd* Pool,
 		 cSecond creationtime,
 		 cSecond lifetime,
 		 double power, double HMMass, double BurnUp, double ChargeFactor ):CLASSFacility(log, creationtime, lifetime, 4)
{
	(*this).SetName("R_Reactor.");


	fIsStarted = false;
	fIsShutDown = false;
	fIsAtEndOfCycle = false;

	fStorage = 0;
	fFuelTypeDB = 0;
	fFabricationPlant = 0;

	fFixedFuel = true;
	fIsStorage = false;

	fOutBackEndFacility = Pool;

	fPower = power * ChargeFactor;

	fHeavyMetalMass = HMMass;

	map<ZAI ,double>::iterator it;
	map<ZAI ,double> isotopicquantity = evolutivedb.GetIsotopicVectorAt(0.).GetActinidesComposition().GetIsotopicQuantity();
	double M0 = 0;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
		M0 += (*it).second*cZAIMass.fZAIMass.find( (*it).first )->second/AVOGADRO*1e-6;

	fEvolutionDB = evolutivedb * (fHeavyMetalMass/M0);

	fBurnUp = BurnUp;
	fCycleTime = (cSecond) (fBurnUp*1e9 / (fPower)  * fHeavyMetalMass  *3600*24);

	fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVInCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fCycleTime/fEvolutionDB.GetPower()*fPower) );



	fNextPlan = fLoadingPlan.begin();

	INFO << " A Reactor has been define :" << endl;
	INFO << "\t Fuel Composition is fixed ! "<< endl;
	INFO << "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	INFO << "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl;
	INFO << "\t The Cycle Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	INFO << "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << power << " and " << ChargeFactor << " Charge Factor)"<< endl;
	INFO << "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;


}


//________________________________________________________________________
Reactor::~Reactor()
{


}

//________________________________________________________________________
void Reactor::SetCycleTime(double cycletime)
{
	fCycleTime = (cSecond)cycletime;

	if(fFixedFuel==true)
	{
		fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt(fCycleTime/fEvolutionDB.GetPower()*fPower);
		fBurnUp = fPower*fCycleTime/3600./24./fHeavyMetalMass;
	}
	else
	{
		fBurnUp = fPower*fCycleTime/3600./24./fHeavyMetalMass;
	}
}
//________________________________________________________________________
void Reactor::SetPower(double Power)
{
	fPower = Power;

	if(fFixedFuel==true)
	{
		fCycleTime = (cSecond) (fBurnUp*1e9 / (fPower)  * fHeavyMetalMass  *3600*24);
		fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fCycleTime/fEvolutionDB.GetPower()*fPower) );
	}
	else
		fCycleTime = (cSecond)(fBurnUp*1e9 / (fPower)  * fHeavyMetalMass  *3600*24);	 //BU in GWd/t


}
//________________________________________________________________________
void Reactor::SetBurnUp(double BU)
{

	fBurnUp = BU;

	if(fFixedFuel==true)
	{
		fCycleTime = (cSecond) (fBurnUp*1e9 / (fPower)  * fHeavyMetalMass  *3600*24);
		fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fCycleTime/fEvolutionDB.GetPower()*fPower) );
	}
	else
		fCycleTime = (cSecond) (fBurnUp*1e9 / (fPower)  * fHeavyMetalMass  *3600*24);

}

//________________________________________________________________________
void Reactor::SetEvolutionDB(EvolutionData evolutionDB)
{
	map<ZAI ,double>::iterator it;
	map<ZAI ,double> isotopicquantity = evolutionDB.GetIsotopicVectorAt(0.).GetActinidesComposition().GetIsotopicQuantity();
	double M0 = 0;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
		M0 += (*it).second*cZAIMass.fZAIMass.find( (*it).first )->second/AVOGADRO*1e-6;
	fEvolutionDB = evolutionDB * (fHeavyMetalMass/M0);

	fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fCycleTime/fEvolutionDB.GetPower()*fPower) );
	fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);


}

//________________________________________________________________________
void Reactor::SetNewFuel(EvolutionData ivdb)
{

	SetEvolutionDB(ivdb);

}

//________________________________________________________________________
void Reactor::Evolution(cSecond t)
{
DBGL


	if( fIsShutDown  || t < GetCreationTime() ) return; // Reactor stop or not started...

	if(Norme(fInsideIV)!=0)
	{
#pragma omp critical(ParcPowerUpdate)
		{GetParc()->AddToPower(fPower);}
	}
	else if(fIsStarted==true)
	{
		WARNING << " Reactor should be working but have no Heavy Nucleus Inside. It's not working so have a zero power..."
		<< " Time : "<< t/365.25/3600/24 << " years" << endl;
	}


	if( t < fInternalTime ) return;
	if( t == fInternalTime && t!=GetCreationTime() ) return;



	if( t == GetCreationTime() && !fIsStarted) // Start of the Reactor
	{
		fIsAtEndOfCycle = true;
		fInsideIV  = fIVBeginCycle;
		fInternalTime = t;
		fInCycleTime = 0;

	}

	// Check if the Reactor if started ...
	if(!fIsStarted) return;			// If the reactor just start don't need to make Fuel evolution


	cSecond EvolutionTime = t - fInternalTime; // Calculation of the evolution time (relativ)


	if( EvolutionTime + fInCycleTime == fCycleTime )		//End of Cycle
	{
		fIsAtEndOfCycle = true;
		fInternalTime += EvolutionTime; 				// Update Internal Time
		fInCycleTime += EvolutionTime;					// Update InCycleTime

		if(t >=  GetCreationTime() + GetLifeTime())				// if the Next Cycle don't 'Exist...
			fIsShutDown = true;


	}
	else if(EvolutionTime + fInCycleTime < fCycleTime )			// During Cycle
	{

		fInternalTime += EvolutionTime;					// Update Internal Time
		fInCycleTime += EvolutionTime;					// Update InCycleTime

		fInsideIV = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fInCycleTime/fEvolutionDB.GetPower()*fPower) );	// update the fuel composition
		if(t>=GetCreationTime() + GetLifeTime())	fIsShutDown = true;
	}
	else
	{
		// This is so bad!! You will probably unsynchronize all the reactor....
		ERROR << " Evolution is too long! There is a problem in Reactor evolution at " << t/365.25/3600/24 << endl;
		exit(1);
	}

DBGL
}

//________________________________________________________________________
void Reactor::Dump()
{
DBGL


	if(fInternalTime < GetCreationTime()) return;
	if(fIsShutDown  && !fIsStarted) return; // Reactor stopped...

	if(fFixedFuel )
	{
		if(fIsAtEndOfCycle  && !fIsShutDown )
		{
			fIsAtEndOfCycle = false;

			if(fIsStarted  )					// A Cycle has already been done
			{
				fOutBackEndFacility->AddIV(fInsideIV);
				AddCumulativeIVOut(fInsideIV);
			}
			else fIsStarted = true;					// Just start the first cycle

			if(!GetParc()->GetStockManagement() && fIsStorage )
			{
				IsotopicVector BuildIVtmp ;
				IsotopicVector GodPart;

				//Get The Storage Compostion
				BuildIVtmp.Add(fStorage->GetInsideIV().GetIsotopicQuantity());
				//Get the rest after IVIn creation
				BuildIVtmp -= fIVInCycle;
				//Get the God part form this rest
				GodPart.Add(BuildIVtmp.GetIsotopicQuantityNeeded()) ;
				//Take what you can from Storage...
				fStorage->TakeFromStock( fIVInCycle - GodPart);
				//And Get the rest from God
				GetParc()->AddGod(GodPart);

			}
			else	GetParc()->AddGod(fIVInCycle);


			if(fNextPlan != fLoadingPlan.end())		// Check if the Fuel change
			{
				if(fInternalTime >= (*fNextPlan).first)
				{
					SetEvolutionDB((*fNextPlan).second.first);
					SetBurnUp((*fNextPlan).second.second);
					fNextPlan++;
				}
			}
			fInsideIV  = fIVBeginCycle;
			AddCumulativeIVIn(fIVBeginCycle);

			fInCycleTime = 0;
		}
		else if (fIsAtEndOfCycle  && fIsShutDown )		//shutdown at end of Cycle
		{

			fOutBackEndFacility->AddIV(fIVOutCycle);
			AddCumulativeIVOut(fIVOutCycle);
			fInsideIV.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
		else if (!fIsAtEndOfCycle && fIsShutDown ) 					//shutdown during Cycle
		{
			fOutBackEndFacility->AddIV(fInsideIV);
			AddCumulativeIVOut(fInsideIV);
			fInsideIV.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
	}
	else
	{
		if(!GetParc()->GetStockManagement())
		{
			ERROR << " Can't have unfixedFuel without stock management" << endl;
			exit(1);
		}


		if(fIsAtEndOfCycle  && !fIsShutDown )
		{
			fIsAtEndOfCycle = false;

			if(fIsStarted  )					// A Cycle has already been done
			{
				fOutBackEndFacility->AddIV(fIVOutCycle);
				AddCumulativeIVOut(fIVOutCycle);
			}
			else fIsStarted = true;					// Just start the first cycle

			SetNewFuel(fFabricationPlant->GetReactorEvolutionDB(GetId()));
			fFabricationPlant->TakeReactorFuel(GetId());


			fInsideIV  = fIVBeginCycle;
			AddCumulativeIVIn(fIVBeginCycle);

			fInCycleTime = 0;

		}
		else if (fIsAtEndOfCycle  && fIsShutDown )		//shutdown at end of Cycle
		{
			fOutBackEndFacility->AddIV(fIVOutCycle);
			AddCumulativeIVOut(fIVOutCycle);
			fInsideIV.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
		else if (!fIsAtEndOfCycle && fIsShutDown ) 					//shutdown during Cycle
		{
			fOutBackEndFacility->AddIV(fInsideIV);
			AddCumulativeIVOut(fInsideIV);
			fInsideIV.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
		
		
		
	}
	
DBGL
}




