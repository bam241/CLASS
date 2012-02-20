#include "CLASSHeaders.hxx"
#include "Reactor.hxx"
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
	IsStarted = false;
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
	// Check if the TF has been created ...
	if(t<fCreationTime) return;
	DBGL;
	
	
	if(fInternalTime == 0 && IsStarted == false)
	{
		IsStarted = true;
		
		fAssociedTreatmentFactory->AddIVGodIncome(fIVBeginCycle - fParc->BuildIsotopicVector(fIVBeginCycle));
		fIVReactor  = fIVBeginCycle;
		
		
		fInternalTime = fCreationTime;
	}
	
	long int EvolutionTime = t - fInternalTime;


	if(EvolutionTime + fInCycleTime < fCycleTime )
	{
		fIVReactor = fEvolutionDB->GetIsotopicVectorAt(EvolutionTime + fInCycleTime);
		fInternalTime += EvolutionTime;
		fInCycleTime += EvolutionTime;

	}
	else if(EvolutionTime + fInCycleTime == fCycleTime )
	{

		IsotopicVector IVtoCool = fIVOutCycle;

		fInternalTime += EvolutionTime;
		fAssociedTreatmentFactory->Evolution(fInternalTime); //Just for Safety... 

		fAssociedTreatmentFactory->AddIVCooling(IVtoCool);
		if(t < fCreationTime + fLifeTime)
		{

		fAssociedTreatmentFactory->AddIVGodIncome(fIVInCycle - fParc->BuildIsotopicVector(fIVInCycle));
		fIVReactor  = fIVBeginCycle;
		fInCycleTime = 0;
		}
		else 
		{
			IsotopicVector IV;
			fIVReactor = IV;
		}

	}
	else
	{

		// This is so bad!! You will probably unsynchronize all the reactor....
		cout << "!!Warning!! !!!Reactor!!! Evolution is too log! This is a Bad way to deal the evolution of the reactor..."<< t/365.4/3600/24 << " :" << endl;
		exit (1);

//		cout << "!!Warning!! !!!Reactor!!! Will do it anyway, but you should change something in your evolution method...." << endl;
//		
//		int EndToLastCycle = fCycleTime - fInCycleTime;
//		IsotopicVector IVtoCool = fEvolutionDB->GetIsotopicVectorAt(fCycleTime);
//		fInternalTime += EndToLastCycle;
//		EvolutionTime -= EndToLastCycle;
//		fAssociedTreatmentFactory->Evolution(fInternalTime); 
//		fAssociedTreatmentFactory->AddIVCooling(IVtoCool);
//		fAssociedTreatmentFactory->TakeFromStock(fIVBeginCycle,fInternalTime);
//				
//		int Nb_Cycle = EvolutionTime / fCycleTime;
//		int RemainingTime = EvolutionTime % fCycleTime;
//		for(int i = 0; i< Nb_Cycle; i++)
//		{
//			IsotopicVector IVtoCool = fEvolutionDB->GetIsotopicVectorAt(fCycleTime);
//			fInternalTime += fCycleTime;
//			fAssociedTreatmentFactory->Evolution(fInternalTime); //Ohhh this is bad... don't do it !!' 
//			fAssociedTreatmentFactory->AddIVCooling(IVtoCool);
//			fAssociedTreatmentFactory->TakeFromStock(fIVBeginCycle,fInternalTime);
//		}
//		fIVReactor = fEvolutionDB->GetIsotopicVectorAt(RemainingTime);
//		fInternalTime += RemainingTime;
//		fInCycleTime += RemainingTime;
		
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



