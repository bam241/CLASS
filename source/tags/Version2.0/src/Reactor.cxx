#include "Reactor.hxx"

#include "EvolutionData.hxx"
#include "DataBank.hxx"
#include "Pool.hxx"
#include "FabricationPlant.hxx"
#include "Storage.hxx"
#include "CLASS.hxx"

#include "LogFile.hxx"

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
	
	
}

Reactor::Reactor(LogFile* log)
{
	
	fLog = log;
	
}

Reactor::Reactor(LogFile* log, DataBank<IsotopicVector>* 	fueltypeDB,
		 FabricationPlant* fabricationplant,
 		 Pool* Pool,
 		 double creationtime, double lifetime)
{
	
	fLog = log;
	
	fIsStarted = false;
	fShutDown = false;
	fEndOfCycle = false;
	
	fFabricationPlant = fabricationplant;
	fFixedFuel = false;
	fBurnUp = -1.;
	fHeavyMetalMass = -1.;
	
	fAssociedPool = Pool;
	
	fFuelTypeDB = fueltypeDB;
	
	fInternalTime = 0;
	fInCycleTime = 0;
	fPower = -1.;
	fCycleTime = -1.;	 //BU in GWd/t
	
	fCreationTime = (cSecond)creationtime;
	fLifeTime = (cSecond)lifetime;
	
	cout	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	cout	<< "\t Fuel Composition is not fixed ! "<< endl;
	cout	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	cout	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(fLifeTime/3600/24/365.25) << " year" << endl << endl;
	cout	<< "!!WARNING!! !!!Reactor!!! You need to set Burn-up/Power/CycleTime (2 of 3) & Heavy Metal Mass Manualy !! " << endl;

	
	fLog->fLog	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	fLog->fLog	<< "\t Fuel Composition is not fixed ! "<< endl;
	fLog->fLog	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	fLog->fLog	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "\t Life time (Operating's Duration) set at \t " << (double)(fLifeTime/3600/24/365.25) << " year" << endl << endl;
	fLog->fLog	<< "!!WARNING!! !!!Reactor!!! You need to set Burn-up/Power/CycleTime (2 of 3) & Heavy Metal Mass Manualy !! " << endl;

	
}

Reactor::Reactor(LogFile* log, DataBank<IsotopicVector>* 	fueltypeDB, FabricationPlant* fabricationplant, Pool* Pool,
 		 double creationtime, double lifetime,
 		 double Power, double HMMass, double BurnUp, double ChargeFactor)
{
	
	fLog = log;
	
	fIsStarted = false;
	fShutDown = false;
	fEndOfCycle = false;
	
	fFabricationPlant = fabricationplant;
	fFixedFuel = false;
	
	fAssociedPool = Pool;
	
	fFuelTypeDB = fueltypeDB;
	
	fInternalTime = 0;
	fInCycleTime = 0;
	
	fBurnUp = BurnUp;
	fHeavyMetalMass = HMMass;
	fPower = Power*ChargeFactor;
	fCycleTime = (cSecond) (fBurnUp*1e9 / (fPower)  * fHeavyMetalMass  *3600*24);	 //BU in GWd/t

	fCreationTime = (cSecond)creationtime;
	fLifeTime = (cSecond)lifetime;
	
	
	
	cout	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	cout	<< "\t Fuel Composition is not fixed ! "<< endl;
	cout	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	cout	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(fLifeTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << Power << " and " << ChargeFactor << " Charge Factor)"<< endl;
	cout	<< "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	cout	<< "\t The corresponding Cycle Time is\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;
	
	fLog->fLog 	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	fLog->fLog 	<< "\t Fuel Composition is not fixed ! "<< endl;
	fLog->fLog 	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	fLog->fLog 	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	fLog->fLog 	<< "\t Life time (Operating's Duration) set at \t " << (double)(fLifeTime/3600/24/365.25) << " year" << endl;
	fLog->fLog 	<< "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << Power << " and " << ChargeFactor << " Charge Factor)"<< endl;
	fLog->fLog 	<< "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	fLog->fLog 	<< "\t The corresponding Cycle Time is\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	fLog->fLog 	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;
	
	
}

Reactor::Reactor(LogFile* log, DataBank<IsotopicVector>* 	fueltypeDB,
		 FabricationPlant* fabricationplant,
 		 Pool* Pool,
 		 double creationtime, double lifetime, double cycletime,
 		 double HMMass, double BurnUp)
{
	
	fLog = log;

	fIsStarted = false;
	fShutDown = false;
	fEndOfCycle = false;

	fFabricationPlant = fabricationplant;
	fFixedFuel = false;
	fBurnUp = BurnUp;
	fHeavyMetalMass = HMMass;

	fAssociedPool = Pool;

	fFuelTypeDB = fueltypeDB;

	fInternalTime = 0;
	fInCycleTime = 0;
	fCycleTime = (cSecond)cycletime;
	fCreationTime = (cSecond)creationtime;
	fLifeTime = (cSecond)lifetime;
	fPower = BurnUp*3600.*24. / (fCycleTime) * HMMass *1e9; //BU in GWd/t
	
	
	cout	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	cout	<< "\t Fuel Composition is not fixed ! "<< endl;
	cout	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	cout	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Cycle Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	cout	<< "\t The corresponding Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW" << endl;
	cout	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;

	fLog->fLog 	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	fLog->fLog	<< "\t Fuel Composition is not fixed ! "<< endl;
	fLog->fLog	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	fLog->fLog	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "\t Life time (Operating's Duration) set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "\t The Cycle Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	fLog->fLog	<< "\t The corresponding Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW" << endl;
	fLog->fLog	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;
	
	
}


Reactor::Reactor(LogFile* log, EvolutionData evolutivedb,
 		 Pool* Pool,
 		 double creationtime,
 		 double lifetime,
 		 double power, double HMMass, double BurnUp, double ChargeFactor )
{
	
	fLog = log;
	
	fIsStarted = false;
	fShutDown = false;
	fEndOfCycle = false;
	
		

	
	


	fFixedFuel = true;
	fIsStorage = false;
	
	fAssociedPool = Pool;
	
	fInternalTime = 0;
	fInCycleTime = 0;
	fCreationTime = (cSecond)creationtime;
	fLifeTime = (cSecond)lifetime;
	
	fPower = power * ChargeFactor;
	
	fHeavyMetalMass = HMMass;
	
	map<ZAI, double> ZAImass;
	ZAImass.insert( pair< ZAI,double >( ZAI(92,238,0), 238050788.247e-6 ) );
	ZAImass.insert( pair< ZAI,double >( ZAI(92,235,0), 235043929.918e-6 ) );
	ZAImass.insert( pair< ZAI,double >( ZAI(94,238,0), 238049559.894e-6 ) );
	ZAImass.insert( pair< ZAI,double >( ZAI(94,239,0), 239052163.381e-6 ) );
	ZAImass.insert( pair< ZAI,double >( ZAI(94,240,0), 240053813.545e-6 ) );
	ZAImass.insert( pair< ZAI,double >( ZAI(94,241,0), 241056851.456e-6 ) );
	ZAImass.insert( pair< ZAI,double >( ZAI(94,242,0), 242058742.611e-6 ) );
	ZAImass.insert( pair< ZAI,double >( ZAI(95,241,0), 241056829.144e-6 ) );
	
	double Na = 6.02214129e23;	//N Avogadro
	map<ZAI ,double>::iterator it;
	map<ZAI ,double> isotopicquantity = evolutivedb.GetIsotopicVectorAt(0.).GetActinidesComposition().GetIsotopicQuantity();
	double M0 = 0;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
		M0 += (*it).second*ZAImass.find( (*it).first )->second/Na*1e-6;
	
	fEvolutionDB = evolutivedb * (fHeavyMetalMass/M0);
	
	fBurnUp = BurnUp;
	fCycleTime = (cSecond) (fBurnUp*1e9 / (fPower)  * fHeavyMetalMass  *3600*24);

	fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVInCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fCycleTime/fEvolutionDB.GetPower()*fPower) );
		
	cout	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	cout	<< "\t Fuel Composition is fixed ! "<< endl;
	cout	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(fLifeTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Cycle Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << power << " and " << ChargeFactor << " Charge Factor)"<< endl;
	cout	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;

	fLog->fLog	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	fLog->fLog	<< "\t Fuel Composition is fixed ! "<< endl;
	fLog->fLog	<< "\t Creation time set at \t " << (double)(fCreationTime/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "\t Life time (Operating's Duration) set at \t " << (double)(fLifeTime/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "\t The Cycle Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << power << " and " << ChargeFactor << " Charge Factor)"<< endl;
	fLog->fLog	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;

	
}


//________________________________________________________________________
Reactor::~Reactor()
{
	
	
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
void Reactor::SetEvolutionDB(EvolutionData evolutionDB)
{
	
	fEvolutionDB = evolutionDB;
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
	

	if( fShutDown == true || t < fCreationTime ) return; // Reactor stop or not started...

	if(Norme(fIVReactor)!=0){
#pragma omp critical(ParcPowerUpdate)
		{fParc->AddToPower(fPower);}
	}
	else if(fIsStarted==true){
	fLog->fLog << "!!Warning!! !!!Reactor!!!"
		   << " Reactor should be working but there is no Heavy Nucleus Inside. It's not working so have a zero power..."
		   << " Time : "<< t/365.25/3600/24 << " years" << endl;	
	}

	
	if( t == fInternalTime && t!=0 ) return; 
	
	

	if(fInternalTime == 0 && fIsStarted == false) // Start of the Reactor
	{
		fEndOfCycle = true;
		fIVReactor  = fIVBeginCycle;
		fInternalTime = t;
		
	}

	// Check if the Reactor if started ...
	if(fIsStarted == false) return;			// If the reactor just start don't need to make Fuel evolution 


	cSecond EvolutionTime = t - fInternalTime; // Calculation of the evolution time (relativ)

	if( EvolutionTime + fInCycleTime == fCycleTime )		//End of Cycle
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

	
}

//________________________________________________________________________
void Reactor::Dump()
{


	if(fInternalTime < fCreationTime) return;
	if(fShutDown == true && fIsStarted == false) return; // Reactor stopped...

	if(fFixedFuel == true)
	{

		if(fEndOfCycle == true && fShutDown == false )
		{
			fEndOfCycle = false;

			if(fIsStarted == true )					// A Cycle has already been done
				fAssociedPool->AddIVCooling(fIVReactor);
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

			fAssociedPool->AddIVCooling(fIVOutCycle);
			fIVReactor.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
		else if (fEndOfCycle == false && fShutDown == true) 					//shutdown during Cycle
		{
			fAssociedPool->AddIVCooling(fIVReactor);
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
				fAssociedPool->AddIVCooling(fIVOutCycle);
			}
			else fIsStarted = true;					// Just start the first cycle

			SetNewFuel(fFabricationPlant->GetReactorEvolutionDB(fId));
			fFabricationPlant->TakeReactorFuel(fId);


			fIVReactor  = fIVBeginCycle;
			fInCycleTime = 0;

		}
		else if (fEndOfCycle == true && fShutDown == true)		//shutdown at end of Cycle
		{
			fAssociedPool->AddIVCooling(fIVOutCycle);
			fIVReactor.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
		else if (fEndOfCycle == false && fShutDown == true) 					//shutdown during Cycle
		{
			fAssociedPool->AddIVCooling(fIVReactor);
			fIVReactor.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
		


	}


}




