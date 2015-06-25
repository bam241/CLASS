#include "Reactor.hxx"

#include "EvolutionData.hxx"
#include "Pool.hxx"
#include "FabricationPlant.hxx"
#include "Storage.hxx"
#include "Scenario.hxx"
#include "CLASSConstante.hxx"

#include <iostream>
#include <cmath>
#include <omp.h>
#include <typeinfo>

//________________________________________________________________________
//
//		Reactor
//
//
//
//
//________________________________________________________________________



ClassImp(Reactor)


Reactor::Reactor():CLASSFacility(4)
{
	
	SetName("R_Reactor.");
	
	fOutBackEndFacility = 0;
	fStorage = 0;
	fFabricationPlant = 0;
	fFuelPlan = 0;
	
}

Reactor::Reactor(CLASSLogger* log):CLASSFacility(log, 4)
{
	DBGL
	
	fOutBackEndFacility = 0;
	fStorage = 0;
	fFabricationPlant = 0;
	fFuelPlan = 0;
	SetName("R_Reactor.");
	
	DBGL
}

Reactor::Reactor(CLASSLogger* log,
		 CLASSBackEnd* Pool,
		 cSecond creationtime,
		 cSecond lifetime,
		 double power, double HMMass, double CapacityFactor ):CLASSFacility(log, creationtime, lifetime, 4)
{
	DBGL
	(*this).SetName("R_Reactor.");
	
	
	fIsStarted = false;
	fIsShutDown = false;
	fIsAtEndOfCycle = false;
	
	fStorage = 0;
	fFabricationPlant = 0;
	
	fFixedFuel = true;
	fIsStorage = false;
	
	fOutBackEndFacility = Pool;
	
	fPower = power * CapacityFactor;
	fEfficiencyFactor = 0.33;
	fElectricPower = fEfficiencyFactor*fPower;
	
	
	fHeavyMetalMass = HMMass;
	
	fBurnUp = -1;
	fCycleTime = (-1);
	
	fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVInCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fCycleTime/fEvolutionDB.GetPower()*fPower) );
	
	
	fFuelPlan = 0;
	
	INFO << " A Reactor has been define :" << endl;
	INFO << "\t Fuel Composition is fixed (for now)! "<< endl;
	INFO << "\t Creation time set at \t " << ((double)GetCreationTime())/((double)cYear) << " year" << endl;
	INFO << "\t Life time (Operating's Duration) set at \t " << ((double)GetLifeTime())/((double)cYear) << " year" << endl;
	INFO << "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << power << " and " << CapacityFactor << " capacity factor)"<< endl;
	INFO << "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;
	
	DBGL
	
}

Reactor::Reactor(CLASSLogger* log,
		 FabricationPlant* fabricationplant, CLASSBackEnd* Pool,
		 cSecond creationtime, cSecond lifetime,
		 double Power, double HMMass, double CapacityFactor):CLASSFacility(log, creationtime, lifetime, 4)
{
	DBGL
	(*this).SetName("R_Reactor.");
	
	
	fStorage = 0;
	fIsStarted = false;
	fIsShutDown = false;
	fIsAtEndOfCycle = false;
	
	fFabricationPlant = fabricationplant;
	fFixedFuel = false;
	
	fOutBackEndFacility = Pool;
	
	fBurnUp = -1;
	fHeavyMetalMass = HMMass;
	fPower = Power*CapacityFactor;
	fEfficiencyFactor = 0.33;
	fElectricPower = fEfficiencyFactor*fPower;
	
	fCycleTime = -1;	 //BU in GWd/t
	
	fFuelPlan = 0;
	
	
	
	INFO << " A Reactor has been define :" << endl;
	INFO << "\t Fuel Composition is not fixed (for now)! "<< endl;
	INFO << "\t Creation time set at \t " << ((double)GetCreationTime())/((double)cYear) << " year" << endl;
	INFO << "\t Life time (Operating's Duration) set at \t " << ((double)GetLifeTime())/((double)cYear) << " year" << endl;
	INFO << "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << Power << " and " << CapacityFactor << " capacity factor)"<< endl;
	INFO << "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;
	
	
	DBGL
	
}


Reactor::Reactor(CLASSLogger* log, PhysicsModels* fueltypeDB, FabricationPlant* fabricationplant, CLASSBackEnd* Pool,
		 cSecond creationtime, cSecond lifetime,
		 double Power, double HMMass, double BurnUp, double CapacityFactor):CLASSFacility(log, creationtime, lifetime, 4)
{
	DBGL
	(*this).SetName("R_Reactor.");
	
	
	fStorage = 0;
	fIsStarted = false;
	fIsShutDown = false;
	fIsAtEndOfCycle = false;
	
	fFabricationPlant = fabricationplant;
	
	fFixedFuel = false;
	
	fOutBackEndFacility = Pool;
	
	fBurnUp = BurnUp;
	fHeavyMetalMass = HMMass;
	fPower = Power*CapacityFactor;
	fEfficiencyFactor = 0.33;
	fElectricPower = fEfficiencyFactor*fPower;
	fCycleTime = (cSecond) (fBurnUp*1e9 / (fPower)  * fHeavyMetalMass  *3600*24);	 //BU in GWd/t
	
	fFuelPlan = new CLASSFuelPlan(log);
	fFuelPlan->AddFuel(creationtime, CLASSFuel(fueltypeDB), fBurnUp);
	
	
	
	INFO << " A Reactor has been define :" << endl;
	INFO << "\t Fuel Composition is not fixed ! "<< endl;
	INFO << "\t Creation time set at \t " << (double)(GetCreationTime()/cYear) << " year" << endl;
	INFO << "\t Life time (Operating's Duration) set at \t " << ((double)GetLifeTime())/((double)cYear) << " year" << endl;
	INFO << "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << Power << " and " << CapacityFactor << " capacity factor)"<< endl;
	INFO << "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	INFO << "\t The corresponding Cycle Time is\t " << ((double)fCycleTime)/((double)cYear) << " year" << endl;
	INFO << "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;
	
	
	DBGL
	
}

Reactor::Reactor(CLASSLogger* log, PhysicsModels* 	fueltypeDB,
		 FabricationPlant* fabricationplant,
		 CLASSBackEnd* Pool,
		 cSecond creationtime, cSecond lifetime, cSecond cycletime,
		 double HMMass, double BurnUp):CLASSFacility(log, creationtime, lifetime, cycletime, 4)
{
	DBGL
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
	fPower = BurnUp*3600.*24. / (fCycleTime) * HMMass *1e9; //BU in GWd/t
	fEfficiencyFactor = 0.33;
	fElectricPower = fEfficiencyFactor*fPower;
	
	fFuelPlan = new CLASSFuelPlan(log);
	fFuelPlan->AddFuel(creationtime, CLASSFuel(fueltypeDB), fBurnUp);
	
	
	INFO << " A Reactor has been define :" << endl;
	INFO << "\t Fuel Composition is not fixed ! "<< endl;
	INFO << "\t Creation time set at \t " << ((double)GetCreationTime())/((double)cYear) << " year" << endl;
	INFO << "\t Life time (Operating's Duration) set at \t " << ((double)GetCreationTime())/((double)cYear) << " year" << endl;
	INFO << "\t The Cycle Time set at\t " << ((double)fCycleTime)/((double)cYear) << " year" << endl;
	INFO << "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	INFO << "\t The corresponding Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW" << endl;
	INFO << "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;
	
	
	DBGL
	
}


Reactor::Reactor(CLASSLogger* log, EvolutionData* evolutivedb,
		 CLASSBackEnd* Pool,
		 cSecond creationtime,
		 cSecond lifetime,
		 double power, double HMMass, double BurnUp, double CapacityFactor ):CLASSFacility(log, creationtime, lifetime, 4)
{
	DBGL
	(*this).SetName("R_Reactor.");
	
	
	fIsStarted = false;
	fIsShutDown = false;
	fIsAtEndOfCycle = false;
	
	fStorage = 0;
	fFabricationPlant = 0;
	
	fFixedFuel = true;
	fIsStorage = false;
	
	fOutBackEndFacility = Pool;
	
	fPower = power * CapacityFactor;
	fEfficiencyFactor = 0.33;
	fElectricPower = fEfficiencyFactor*fPower;
	
	fHeavyMetalMass = HMMass;
	
	double M0 = cZAIMass.GetMass( evolutivedb->GetIsotopicVectorAt(0.).GetActinidesComposition() );
	
	fEvolutionDB = (*evolutivedb) * (fHeavyMetalMass/M0);
	
	fBurnUp = BurnUp;
	fCycleTime = (cSecond) (fBurnUp*1e9 / (fPower)  * fHeavyMetalMass  *3600*24);
	
	fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVInCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fCycleTime/fEvolutionDB.GetPower()*fPower) );
	
	fFuelPlan = new CLASSFuelPlan(log);
	fFuelPlan->AddFuel(creationtime, CLASSFuel(evolutivedb), fBurnUp);
	
	INFO << " A Reactor has been define :" << endl;
	INFO << "\t Fuel Composition is fixed ! "<< endl;
	INFO << "\t Creation time set at \t " << (double)(GetCreationTime()/cYear) << " year" << endl;
	INFO << "\t Life time (Operating's Duration) set at \t " << ((double)GetLifeTime())/((double)cYear) << " year" << endl;
	INFO << "\t The Cycle Time set at\t " << ((double)fCycleTime)/((double)cYear) << " year" << endl;
	INFO << "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << power << " and " << CapacityFactor << " capacity factor)"<< endl;
	INFO << "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;
	
	
	DBGL
}

Reactor::Reactor(CLASSLogger* log, EvolutionData* evolutivedb,
		 CLASSBackEnd* Pool,
		 cSecond creationtime, cSecond lifetime, cSecond cycletime,
		 double HMMass, double BurnUp ):CLASSFacility(log, creationtime, lifetime, cycletime, 4)
{
	DBGL
	(*this).SetName("R_Reactor.");
	
	
	fIsStarted = false;
	fIsShutDown = false;
	fIsAtEndOfCycle = false;
	
	fStorage = 0;
	fFabricationPlant = 0;
	
	fFixedFuel = true;
	fIsStorage = false;
	
	fOutBackEndFacility = Pool;
	
	fPower = BurnUp*3600.*24. / (fCycleTime) * HMMass *1e9; //BU in GWd/t
	fEfficiencyFactor = 0.33;
	fElectricPower = fEfficiencyFactor*fPower;
	
	fHeavyMetalMass = HMMass;
	
	double M0 = cZAIMass.GetMass( evolutivedb->GetIsotopicVectorAt(0.).GetActinidesComposition() );
	
	fEvolutionDB = (*evolutivedb) * (fHeavyMetalMass/M0);
	
	fBurnUp = BurnUp;
	
	fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVInCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fCycleTime/fEvolutionDB.GetPower()*fPower) );
	
	
	fFuelPlan = new CLASSFuelPlan(log);
	fFuelPlan->AddFuel(creationtime, CLASSFuel(evolutivedb), fBurnUp);
	
	INFO << " A Reactor has been define :" << endl;
	INFO << "\t Fuel Composition is fixed ! "<< endl;
	INFO << "\t Creation time set at \t " << ((double)GetCreationTime())/((double)cYear) << " year" << endl;
	INFO << "\t Life time (Operating's Duration) set at \t " << ((double)GetLifeTime())/((double)cYear) << " year" << endl;
	INFO << "\t The Cycle Time set at\t " << ((double)fCycleTime)/((double)cYear) << " year" << endl;
	INFO << "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << fPower << endl;
	INFO << "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;
	
	DBGL
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
	DBGL
	double M0 = cZAIMass.GetMass( evolutionDB.GetIsotopicVectorAt(0.).GetActinidesComposition() );
	fEvolutionDB = evolutionDB * (fHeavyMetalMass/M0);
	
	fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fCycleTime/fEvolutionDB.GetPower()*fPower) );
	fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);

	DBGL
	
}

//________________________________________________________________________
void Reactor::SetNewFuel(EvolutionData ivdb)
{
	DBGL
	SetEvolutionDB(ivdb);
	DBGL
}

//________________________________________________________________________
void Reactor::Evolution(cSecond t)
{
	DBGL
	
	
	if( fIsShutDown  || t < GetCreationTime() ) return; // Reactor stop or not started...
	
	if(Norme(fInsideIV)!=0 && fIsStarted)
	{
#pragma omp critical(ParcPowerUpdate)
		{GetParc()->AddToPower(fPower, fElectricPower);}
	}
	else if(fIsStarted)
	{
		WARNING << " Reactor should be working but have no Heavy Nucleus Inside. It's not working so have a zero power..."
		<< " Time : "<< t/cYear << " years" << endl;
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
		// Evolution goes after the end of cycle.... check it
		ERROR << " " << (*this).GetName() << endl;
		ERROR << " Evolution Time is " << t << endl;
		ERROR << " Evolution is too long! There is a problem in Reactor evolution at " << t/cYear << endl;
		ERROR << " This is too long of : " << EvolutionTime + fInCycleTime - fCycleTime << endl;
		ERROR << " I have spend " << fInCycleTime +EvolutionTime << " and should have been " << fCycleTime << endl;
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
	if(!fIsAtEndOfCycle && !fIsShutDown) return;
	
	// First trash the irradiated fuel
	if(fIsAtEndOfCycle  && !fIsShutDown )
	{
		if(fIsStarted  )					// A Cycle has already been done
		{
			fOutBackEndFacility->AddIV(fInsideIV);
			AddCumulativeIVOut(fInsideIV);
		}
		else fIsStarted = true;					// Just start the first cycle
		
	}
	else if (fIsAtEndOfCycle  && fIsShutDown )			//shutdown at end of Cycle
	{
		
		fOutBackEndFacility->AddIV(fIVOutCycle);
		AddCumulativeIVOut(fIVOutCycle);
		fInsideIV.Clear();
		fInCycleTime = 0;
		fIsStarted = false;					// shut down the Reactor
	}
	else if (!fIsAtEndOfCycle && fIsShutDown ) 			//shutdown during Cycle
	{
		fOutBackEndFacility->AddIV(fInsideIV);
		AddCumulativeIVOut(fInsideIV);
		fInsideIV.Clear();
		fInCycleTime = 0;
		fIsStarted = false;					// shut down the Reactor
	}
	
	
	DBGL
	
	
	// Get the new Fuel !
	pair<CLASSFuel, double> NextFuel = fFuelPlan->GetFuelAt(fInternalTime);
	SetBurnUp((NextFuel).second);
	
	if( NextFuel.first.GetPhysicsModels() )
		fFixedFuel = false;
	else if( NextFuel.first.GetEvolutionData() )
		fFixedFuel = true;
	else
	{
		ERROR << typeid(NextFuel.first).name() << endl;
		ERROR << "WRONG Fuel Format Correct it !! " << endl;
		exit(1);
	}
	DBGL
	
	if(fFixedFuel )
	{
		DBGL
		if(fIsAtEndOfCycle  && !fIsShutDown )
		{
			SetEvolutionDB( *(NextFuel.first.GetEvolutionData()) );
			fIsAtEndOfCycle = false;
			
			
			if(!GetParc()->GetStockManagement() && fIsStorage )
			{
				IsotopicVector BuildIVtmp ;
				IsotopicVector OutIncomePart;
				
				//Get The Storage Compostion
				BuildIVtmp.Add(fStorage->GetInsideIV().GetIsotopicQuantity());
				//Get the rest after IVIn creation
				BuildIVtmp -= fIVInCycle;
				//Get the OutIncome part form this rest
				OutIncomePart.Add(BuildIVtmp.GetIsotopicQuantityNeeded()) ;
				//Take what you can from Storage...
				fStorage->TakeFromStock( fIVInCycle - OutIncomePart);
				//And Get the rest from OutIncome
				GetParc()->AddOutIncome(OutIncomePart);
				
			}
			else	GetParc()->AddOutIncome(fIVInCycle);
			
			
			fInsideIV  = fIVBeginCycle;
			AddCumulativeIVIn(fIVBeginCycle);
			
			fInCycleTime = 0;
		}
		DBGL
		
	}
	else
	{
		DBGL
		if(!GetParc()->GetStockManagement())
		{
			ERROR << " Can't have unfixedFuel without stock management" << endl;
			exit(1);
		}
		
		if(fIsAtEndOfCycle  && !fIsShutDown )
		{
			fIsAtEndOfCycle = false;
			SetNewFuel(fFabricationPlant->GetReactorEvolutionDB(GetId()));
			fFabricationPlant->TakeReactorFuel(GetId());
			
			fInsideIV  = fIVBeginCycle;
			AddCumulativeIVIn(fIVBeginCycle);
			
			fInCycleTime = 0;
		}
		
		DBGL
		
		
		
	}
	
	DBGL
}



//________________________________________________________________________
cSecond Reactor::GetNextCycleTime(cSecond time)
{
	DBGL
	cSecond LastCycle = fInternalTime - fInCycleTime;
	
	while ( LastCycle < time)
	{
		cSecond cycletime = (cSecond)(fFuelPlan->GetFuelAt(LastCycle).second*1e9 / (fPower)  * fHeavyMetalMass  *3600*24);
		LastCycle += cycletime;
	}
	
	DBGL
	return LastCycle;
}

