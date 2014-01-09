#include "Reactor.hxx"

#include "EvolutionData.hxx"
#include "DataBank.hxx"
#include "Pool.hxx"
#include "FabricationPlant.hxx"
#include "Storage.hxx"
#include "CLASS.hxx"
#include "CLASSHeaders.hxx"

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

	SetLog(log);

}

Reactor::Reactor(LogFile* log, DataBank<IsotopicVector>* 	fueltypeDB,
		 FabricationPlant* fabricationplant,
 		 Pool* Pool,
 		 double creationtime, double lifetime)
{

	SetLog(log);

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

	SetCreationTime( (cSecond)creationtime );
	SetLifeTime( (cSecond)lifetime );

	cout	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	cout	<< "\t Fuel Composition is not fixed ! "<< endl;
	cout	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	cout	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl << endl;
	cout	<< "!!WARNING!! !!!Reactor!!! You need to set Burn-up/Power/CycleTime (2 of 3) & Heavy Metal Mass Manualy !! " << endl;


	GetLog()->fLog	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	GetLog()->fLog	<< "\t Fuel Composition is not fixed ! "<< endl;
	GetLog()->fLog	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	GetLog()->fLog	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl << endl;
	GetLog()->fLog	<< "!!WARNING!! !!!Reactor!!! You need to set Burn-up/Power/CycleTime (2 of 3) & Heavy Metal Mass Manualy !! " << endl;

}

Reactor::Reactor(LogFile* log, DataBank<IsotopicVector>* 	fueltypeDB, FabricationPlant* fabricationplant, Pool* Pool,
 		 double creationtime, double lifetime,
 		 double Power, double HMMass, double BurnUp, double ChargeFactor)
{

	SetLog(log);

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

	SetCreationTime( (cSecond)creationtime );
	SetLifeTime( (cSecond)lifetime );



	cout	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	cout	<< "\t Fuel Composition is not fixed ! "<< endl;
	cout	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	cout	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << Power << " and " << ChargeFactor << " Charge Factor)"<< endl;
	cout	<< "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	cout	<< "\t The corresponding Cycle Time is\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;

	GetLog()->fLog 	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	GetLog()->fLog 	<< "\t Fuel Composition is not fixed ! "<< endl;
	GetLog()->fLog 	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	GetLog()->fLog 	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	GetLog()->fLog 	<< "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl;
	GetLog()->fLog 	<< "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << Power << " and " << ChargeFactor << " Charge Factor)"<< endl;
	GetLog()->fLog 	<< "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	GetLog()->fLog 	<< "\t The corresponding Cycle Time is\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	GetLog()->fLog 	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;



}

Reactor::Reactor(LogFile* log, DataBank<IsotopicVector>* 	fueltypeDB,
		 FabricationPlant* fabricationplant,
 		 Pool* Pool,
 		 double creationtime, double lifetime, double cycletime,
 		 double HMMass, double BurnUp)
{

	SetLog(log);

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
	SetCreationTime( (cSecond)creationtime );
	SetLifeTime( (cSecond)lifetime );
	fPower = BurnUp*3600.*24. / (fCycleTime) * HMMass *1e9; //BU in GWd/t


	cout	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	cout	<< "\t Fuel Composition is not fixed ! "<< endl;
	cout	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	cout	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Cycle Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	cout	<< "\t The corresponding Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW" << endl;
	cout	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;

	GetLog()->fLog 	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	GetLog()->fLog	<< "\t Fuel Composition is not fixed ! "<< endl;
	GetLog()->fLog	<< "\t Fuel Type set to : \t "<<  fFuelTypeDB->GetFuelType() << endl;
	GetLog()->fLog	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t Life time (Operating's Duration) set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t The Cycle Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t Burn-Up at end of Cycle set at \t " << (double)(fBurnUp) << " GWj/t" << endl;
	GetLog()->fLog	<< "\t The corresponding Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW" << endl;
	GetLog()->fLog	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;



}


Reactor::Reactor(LogFile* log, EvolutionData evolutivedb,
 		 Pool* Pool,
 		 double creationtime,
 		 double lifetime,
 		 double power, double HMMass, double BurnUp, double ChargeFactor )
{

	SetLog(log);

	fIsStarted = false;
	fShutDown = false;
	fEndOfCycle = false;



	fFixedFuel = true;
	fIsStorage = false;

	fAssociedPool = Pool;

	fInternalTime = 0;
	fInCycleTime = 0;
	SetCreationTime( (cSecond)creationtime );
	SetLifeTime( (cSecond)lifetime );

	fPower = power * ChargeFactor;

	fHeavyMetalMass = HMMass;

	double Na = 6.02214129e23;	//N Avogadro
	map<ZAI ,double>::iterator it;
	map<ZAI ,double> isotopicquantity = evolutivedb.GetIsotopicVectorAt(0.).GetActinidesComposition().GetIsotopicQuantity();
	double M0 = 0;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
		M0 += (*it).second*cZAIMass.fZAIMass.find( (*it).first )->second/Na*1e-6;

	fEvolutionDB = evolutivedb * (fHeavyMetalMass/M0);

	fBurnUp = BurnUp;
	fCycleTime = (cSecond) (fBurnUp*1e9 / (fPower)  * fHeavyMetalMass  *3600*24);

	fIVBeginCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVInCycle = fEvolutionDB.GetIsotopicVectorAt(0);
	fIVOutCycle = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fCycleTime/fEvolutionDB.GetPower()*fPower) );

	cout	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	cout	<< "\t Fuel Composition is fixed ! "<< endl;
	cout	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	cout	<< "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Cycle Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << power << " and " << ChargeFactor << " Charge Factor)"<< endl;
	cout	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;

	GetLog()->fLog	<< "!!INFO!! !!!Reactor!!! A Reactor has been define :" << endl;
	GetLog()->fLog	<< "\t Fuel Composition is fixed ! "<< endl;
	GetLog()->fLog	<< "\t Creation time set at \t " << (double)(GetCreationTime()/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t Life time (Operating's Duration) set at \t " << (double)(GetLifeTime()/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t The Cycle Time set at\t " << (double)(fCycleTime/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t The Effective Thermal Power is \t " << (double)(fPower *1e-6) << " MW (with Full Power " << power << " and " << ChargeFactor << " Charge Factor)"<< endl;
	GetLog()->fLog	<< "\t The Heavy Metal Mass in the Core set at " << (double)(fHeavyMetalMass) << " tons" << endl << endl;


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
		fBurnUp = fPower*fCycleTime/3600./24./fHeavyMetalMass;
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
	double Na = 6.02214129e23;	//N Avogadro
	map<ZAI ,double>::iterator it;
	map<ZAI ,double> isotopicquantity = evolutionDB.GetIsotopicVectorAt(0.).GetActinidesComposition().GetIsotopicQuantity();
	double M0 = 0;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
		M0 += (*it).second*cZAIMass.fZAIMass.find( (*it).first )->second/Na*1e-6;
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


	if( fShutDown == true || t < GetCreationTime() ) return; // Reactor stop or not started...

	if(Norme(fInsideIV)!=0)
	{
#pragma omp critical(ParcPowerUpdate)
		{GetParc()->AddToPower(fPower);}
	}
	else if(fIsStarted==true)
	{
		GetLog()->fLog << "!!Warning!! !!!Reactor!!!"
		<< " Reactor should be working but there is no Heavy Nucleus Inside. It's not working so have a zero power..."
		<< " Time : "<< t/365.25/3600/24 << " years" << endl;
	}


	if( t == fInternalTime && t!=0 ) return;
	if( t < fInternalTime ) return;



	if(fInternalTime == 0 && fIsStarted == false) // Start of the Reactor
	{
		fEndOfCycle = true;
		fInsideIV  = fIVBeginCycle;
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

		if(t >=  GetCreationTime() + GetLifeTime())				// if the Next Cycle don't 'Exist...
			fShutDown = true;

	}
	else if(EvolutionTime + fInCycleTime < fCycleTime )			// During Cycle
	{

		fInternalTime += EvolutionTime;					// Update Internal Time
		fInCycleTime += EvolutionTime;					// Update InCycleTime

		fInsideIV = fEvolutionDB.GetIsotopicVectorAt( (cSecond)(fInCycleTime/fEvolutionDB.GetPower()*fPower) );	// update the fuel composition
		if(t>=GetCreationTime() + GetLifeTime())	fShutDown = true;
	}
	else
	{
		// This is so bad!! You will probably unsynchronize all the reactor....
		cout << "!!Warning!! !!!Reactor!!!"
		<< " Evolution is too long! This is a Bad way to deal the evolution of the reactor..."
		<< t/365.25/3600/24 << " :" << endl;

		GetLog()->fLog << "!!Warning!! !!!Reactor!!!"
		<< " Evolution is too long! This is a Bad way to deal the evolution of the reactor..."
		<< t/365.25/3600/24 << " :" << endl;
		exit(1);
	}


}

//________________________________________________________________________
void Reactor::Dump()
{


	if(fInternalTime < GetCreationTime()) return;
	if(fShutDown == true && fIsStarted == false) return; // Reactor stopped...

	if(fFixedFuel == true)
	{

		if(fEndOfCycle == true && fShutDown == false )
		{
			fEndOfCycle = false;

			if(fIsStarted == true )					// A Cycle has already been done
				fAssociedPool->AddIVCooling(fInsideIV);
			else fIsStarted = true;					// Just start the first cycle

			if(GetParc()->GetStockManagement() == false && fIsStorage == true)
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
				GetParc()->AddGod(GodPart);

			}
			else	GetParc()->AddGod(fIVInCycle);

			fInsideIV  = fIVBeginCycle;
			fInCycleTime = 0;
		}
		else if (fEndOfCycle == true && fShutDown == true)		//shutdown at end of Cycle
		{

			fAssociedPool->AddIVCooling(fIVOutCycle);
			fInsideIV.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
		else if (fEndOfCycle == false && fShutDown == true) 					//shutdown during Cycle
		{
			fAssociedPool->AddIVCooling(fInsideIV);
			fInsideIV.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
	}
	else
	{
		if(GetParc()->GetStockManagement() == false)
		{
			cout << "!!Warning!! !!!Reactor!!! Can't have unfixedFuel without stock management'" << endl;
			GetLog()->fLog << "!!Warning!! !!!Reactor!!! Can't have unfixedFuel without stock management" << endl;
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

			SetNewFuel(fFabricationPlant->GetReactorEvolutionDB(GetId()));
			fFabricationPlant->TakeReactorFuel(GetId());


			fInsideIV  = fIVBeginCycle;
			fInCycleTime = 0;

		}
		else if (fEndOfCycle == true && fShutDown == true)		//shutdown at end of Cycle
		{
			fAssociedPool->AddIVCooling(fIVOutCycle);
			fInsideIV.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
		else if (fEndOfCycle == false && fShutDown == true) 					//shutdown during Cycle
		{
			fAssociedPool->AddIVCooling(fInsideIV);
			fInsideIV.Clear();
			fInCycleTime = 0;
			fIsStarted = false;		// shut down the Reactor
		}
		
		
		
	}
	
	
}




