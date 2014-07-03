#include "Scenario.hxx"

#include "Storage.hxx"
#include "Reactor.hxx"
#include "CLASSBackEnd.hxx"
#include "Pool.hxx"
#include "FabricationPlant.hxx"
#include "LogFile.hxx"


#include <ctime>
#include "time.h"
#include <cmath>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <omp.h>

//________________________________________________________________________
//
//		CLASS
//
//
//
//
//________________________________________________________________________



float random(float a, float b) //peak random numebr between a and b
{
	float range = pow(2., 31);
	srand(time(NULL)); //initialize the srand
	return (float)a + (float)(b-a)*rand()/range;
}

string dtoa(double num)
{
	ostringstream os(ostringstream::out);
	os<<setprecision(3)<<num;
	return os.str();
}



//________________________________________________________________________
Scenario::Scenario()
{


	fNewTtree = true;
	fPrintStep = (cSecond)(3600*24*365.25);  // One Step per Year
	fAbsoluteTime = 0;
	fStartingTime = 0;

	fStockManagement = true;

	fOutputFileName = "CLASS_Default.root";
	fOutputTreeName = "Data";
	fOutFile = 0;
	fOutT = 0;

	SetLog(new LogFile("CLASS.log"));
	fParcPower = 0;


	// Warning

	cout	<< "!!INFO!! !!!Scenario!!! A Parc has been define :" << endl;
	cout	<< "\t Print step at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	cout	<< "\t Absolute Time set at " << (double)(fAbsoluteTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t StockManagement set at : true" << endl;
	cout	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	cout	<< "\t Log will be in " << GetLog()->GetLogFileName() << endl << endl;

	GetLog()->fLog 	<< "!!INFO!! !!!Scenario!!! Parc has been define :" << endl;
	GetLog()->fLog	<< "\t Print  set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t StockManagement set at : true" << endl;
	GetLog()->fLog	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	GetLog()->fLog	<< "\t Log will be in " << GetLog()->GetLogFileName() << endl << endl;



}
//________________________________________________________________________
Scenario::Scenario(LogFile* Log):CLASSObject(Log)
{


	fNewTtree = true;
	fPrintStep = (cSecond)(3600*24*365.25);  // One Step per Year
	fAbsoluteTime = 0;
	fStartingTime = 0;

	fStockManagement = true;

	fOutputFileName = "CLASS_Default.root";
	fOutputTreeName = "Data";
	fOutFile = 0;
	fOutT = 0;

	fParcPower = 0;


	// Warning

	cout	<< "!!INFO!! !!!Scenario!!! A Parc has been define :" << endl;
	cout	<< "\t Print set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	cout	<< "\t Absolute Time set at " << (double)(fAbsoluteTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t StockManagement set at : true" << endl;
	cout	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	cout	<< "\t Log will be in " << GetLog()->GetLogFileName() << endl << endl;

	GetLog()->fLog 	<< "!!INFO!! !!!Scenario!!! Parc has been define :" << endl;
	GetLog()->fLog	<< "\t Print  set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t StockManagement set at : true" << endl;
	GetLog()->fLog	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	GetLog()->fLog	<< "\t Log will be in " << GetLog()->GetLogFileName() << endl << endl;



}
//________________________________________________________________________
Scenario::Scenario(double abstime)
{

	fNewTtree = true;
	fPrintStep = (cSecond)(3600*24*365.25);  // One Step per Year
	fAbsoluteTime = (cSecond)abstime;
	fStartingTime = fAbsoluteTime;

	fStockManagement = true;

	fOutputFileName = "CLASS_Default.root";
	fOutputTreeName = "Data";
	fOutFile = 0;
	fOutT = 0;

	SetLog(new LogFile("CLASS.log"));
	fParcPower = 0;




	// Warning

	cout	<< "!!INFO!! !!!Scenario!!! A Parc has been define :" << endl;
	cout	<< "\t Print set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	cout	<< "\t StockManagement set at : true" << endl;
	cout	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	cout	<< "\t Log will be in " << GetLog()->GetLogFileName() << endl;

	GetLog()->fLog 	<< "!!INFO!! !!!Scenario!!! Parc has been define :" << endl;
	GetLog()->fLog	<< "\t Print  set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t StockManagement set at : true" << endl;
	GetLog()->fLog	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	GetLog()->fLog	<< "\t Log will be in " << GetLog()->GetLogFileName() << endl << endl;



}


//________________________________________________________________________
Scenario::~Scenario()
{

#pragma omp single
	{CloseOutputTree();}


}

//________________________________________________________________________
void Scenario::AddPool(Pool* Pool)
{


	fPool.push_back(Pool);
	fPool.back()->SetParc(this);
	fPool.back()->SetDecayDataBank( (*this).GetDecayDataBase() );
	fPool.back()->SetLog(GetLog());
	fPool.back()->SetId((int)fPool.size()-1);


	string Pool_name = fPool.back()->GetName();
	if(Pool_name == "P_Pool.")
	{
		Pool_name = "P_Pool";
		Pool_name += dtoa(fPool.back()->GetId());
		Pool_name += ".";
		fPool.back()->SetName(Pool_name.c_str());
	}
	else
	{
		string name_tmp = Pool_name;
		Pool_name = "P_";
		Pool_name += name_tmp;
		Pool_name += ".";
		fPool.back()->SetName(Pool_name.c_str());
	}

	if(!fNewTtree)
		fOutT->Branch(fPool.back()->GetName(), "Pool", &fPool.back());


}

//________________________________________________________________________
void Scenario::AddReactor(Reactor* reactor)
{

	fReactor.push_back(reactor);
	fReactor.back()->SetParc(this);
	fReactor.back()->SetLog(GetLog());
	fReactor.back()->SetId((int)fReactor.size()-1);
	if(!fReactor.back()->IsFuelFixed())
		fReactor.back()->GetFabricationPlant()->AddReactor( (int)fReactor.size()-1,fReactor.back()->GetCreationTime() );


	string Reactor_name = fReactor.back()->GetName();
	if(Reactor_name == "R_Reactor.")
	{
		Reactor_name = "R_Reactor";
		Reactor_name += dtoa(fReactor.back()->GetId());
		Reactor_name += ".";
		fReactor.back()->SetName(Reactor_name.c_str());
	}
	else
	{
		string name_tmp = Reactor_name;
		Reactor_name = "R_";
		Reactor_name += name_tmp;
		Reactor_name += ".";
		fReactor.back()->SetName(Reactor_name.c_str());
	}

	if(!fNewTtree)
		fOutT->Branch(fReactor.back()->GetName(), "Reactor", &fReactor.back());


}

//________________________________________________________________________
void Scenario::AddStorage(Storage* storage)
{

	fStorage.push_back(storage);
	fStorage.back()->SetParc(this);
	fStorage.back()->SetDecayDataBank( (*this).GetDecayDataBase() );
	fStorage.back()->SetLog(GetLog());
	fStorage.back()->SetId((int)fStorage.size()-1);

	string Storage_name = fStorage.back()->GetName();

	if(Storage_name == "S_Storage.")
	{
		Storage_name = "S_Storage";
		Storage_name += dtoa(fStorage.back()->GetId());
		Storage_name += ".";
		fStorage.back()->SetName(Storage_name.c_str());
	}
	else
	{
		string name_tmp = Storage_name;
		Storage_name = "S_";
		Storage_name += name_tmp;
		Storage_name += ".";
		fStorage.back()->SetName(Storage_name.c_str());
	}

	if(!fNewTtree)
		fOutT->Branch(fStorage.back()->GetName(), "Storage", &fStorage.back());


}
//________________________________________________________________________
void Scenario::AddFabricationPlant(FabricationPlant* fabricationplant)
{

	fFabricationPlant.push_back(fabricationplant);
	fFabricationPlant.back()->SetParc(this);
	fFabricationPlant.back()->SetDecayDataBank( (*this).GetDecayDataBase() );
	fFabricationPlant.back()->SetLog(GetLog());
	fFabricationPlant.back()->SetId((int)fStorage.size()-1);


	string FP_name = fFabricationPlant.back()->GetName();
	if(FP_name == "F_FabricationPlant.")
	{
		FP_name = "F_FabricationPlant";
		FP_name += dtoa(fFabricationPlant.back()->GetId());
		FP_name += ".";
		fFabricationPlant.back()->SetName(FP_name.c_str());
	}
	else
	{
		string name_tmp = FP_name;
		FP_name = "F_";
		FP_name += name_tmp;
		FP_name += ".";
		fFabricationPlant.back()->SetName(FP_name.c_str());
	}

	if(!fNewTtree)
		fOutT->Branch(fFabricationPlant.back()->GetName(), "FabricationPlant", &fFabricationPlant.back());
}
//________________________________________________________________________
map<cSecond,int> Scenario::GetTheBackEndTimePath(Reactor* reactor)
{
	cSecond step = 0;
	map<cSecond, int> TheBackEndTimePath;

	{
		pair< map<cSecond, int>::iterator, bool > IResult;
//		IResult = TheBackEndTimePath.insert(pair<cSecond, double> ( step,reactor->GetFacilityType() ) );
//		if( !IResult.second ) IResult.first->second |= reactor->GetFacilityType();

	}



	vector<CLASSBackEnd*> BackEndPath;
	BackEndPath.push_back(reactor->GetOutBackEndFacility());
	while (!BackEndPath.back()->GetStorageType())
	{
		step += BackEndPath.back()->GetCycleTime();
		int FacilityType = BackEndPath.back()->GetFacilityType();
		pair< map<cSecond, int>::iterator, bool > IResult  = TheBackEndTimePath.insert(pair<cSecond,int> (step, FacilityType));
		if( !IResult.second ) IResult.first->second |= FacilityType;

		BackEndPath.push_back(BackEndPath[BackEndPath.size()-1]->GetOutBackEndFacility());

	}

	return TheBackEndTimePath;
}

//________________________________________________________________________
void Scenario::BuildTimeVector(cSecond t)
{
	fTimeStep.clear();
	fTimeStep.insert( pair<cSecond ,int>(t,1) );
	//********* Printing Step *********//
	{
		cSecond step = fStartingTime;

		if(step >= fAbsoluteTime )
			fTimeStep.insert( pair<cSecond ,int>(step,1) );
		step += fPrintStep;
		do
		{

			if(step >= fAbsoluteTime )
				fTimeStep.insert( pair<cSecond ,int>(step,1) );
			step += fPrintStep;
		}
		while( step < t );
	}


	for(int i = 0; i < (int)fReactor.size();i++)
	{
		cSecond ReactorStaringTime = fReactor[i]->GetCreationTime();
		cSecond ReactorShutDownTime = fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime();
		cSecond	ReactorCycleTime = fReactor[i]->GetCycleTime();
		cSecond FabricationCycleTime = 0;
		int ReactorFacilityType = fReactor[i]->GetFacilityType();

		cSecond step = ReactorStaringTime;

		map< cSecond, int > BackEndTimePath = GetTheBackEndTimePath(fReactor[i]);
		if(!fReactor[i]->IsFuelFixed())
			FabricationCycleTime = fReactor[i]->GetFabricationPlant()->GetCycleTime();


		//********* Reactor Evolution Step *********//
		// ShutDown of a reactor

		// Test if the sutdown of the reactor is after the actual time (AbsolutreTime) and before the end of the evolution (t)
		if( ReactorShutDownTime < t )
		{
			//********* Reactor Shutdown *********//
			if( ReactorShutDownTime > fAbsoluteTime)
			{
				pair< map<cSecond, int>::iterator, bool > IResult;
				IResult = fTimeStep.insert( pair<cSecond ,int>(ReactorShutDownTime, 2) );
				if( !IResult.second )
					IResult.first->second |= 2;
			}

			map< cSecond, int >::iterator TV_it; // the time vector iterator
			//********* BackEnd fuel Cycle after reactor Shutdown *********//
			for(TV_it = BackEndTimePath.begin(); TV_it != BackEndTimePath.end(); TV_it++)	// Loop on the BackEnd fuel Cycle Time path
			{
				// Test if each step of the Fuel Cycle BackEnd is after the actual time (AbsolutreTime) and before the end of the evolution (t)
				if( ReactorShutDownTime + (*TV_it).first >= fAbsoluteTime && ReactorShutDownTime + (*TV_it).first <= t)
				{
					pair< map<cSecond, int>::iterator, bool > IResult;
					IResult = fTimeStep.insert( pair<cSecond ,int>(ReactorShutDownTime + (*TV_it).first, (*TV_it).second) );
					if( !IResult.second )
						IResult.first->second |= (*TV_it).second;
				}
			}


		}

		// Start the reactor and the Fuel Fabrication
		if(step >= fAbsoluteTime &&  step <= t && step < ReactorShutDownTime)
		{
			pair< map<cSecond, int>::iterator, bool > IResult;
			IResult = fTimeStep.insert( pair<cSecond ,int>(ReactorStaringTime, ReactorFacilityType) );
			if( !IResult.second )
				IResult.first->second |= ReactorFacilityType;
		}


		//********* FabricationPlant Evolution Step *********//
		if(!fReactor[i]->IsFuelFixed())
		{
			if( (step - FabricationCycleTime) >= fAbsoluteTime && (step - FabricationCycleTime) <= t )
			{
				pair< map<cSecond, int>::iterator, bool > IResult;
				IResult = fTimeStep.insert( pair<cSecond ,int>(step - FabricationCycleTime,16) );
				if( !IResult.second )
					IResult.first->second |= 16;
			}
			else if( step - FabricationCycleTime < fStartingTime )
			{
				cout	   << "!!Warning!! !!!CLASS!!! Can't Build Fuel before Scenario's start\"\n" << endl;
				GetLog()->fLog   << "!!Warning!! !!!CLASS!!! Can't Build Fuel before Scenario's start\"\n" << endl;
				exit(1);
			}
		}

		map<cSecond, pair<EvolutionData, double> >	ReactorLoadingPlan = fReactor[i]->GetLoadingPlan();
		map<cSecond, pair<EvolutionData, double> >::iterator	ReactorNextPlan = ReactorLoadingPlan.begin();




		if (ReactorCycleTime !=0)
		{
			step += ReactorCycleTime;
			do
			{
				if(ReactorNextPlan != ReactorLoadingPlan.end())		// Check if the Fuel change
				{
					if(step >= (*ReactorNextPlan).first)
					{
						ReactorCycleTime = (cSecond) ((*ReactorNextPlan).second.second * 1e9
									      / (fReactor[i]->GetPower())
									      * fReactor[i]->GetHeavyMetalMass()  *3600*24);
						ReactorNextPlan++;

					}
				}

				//********* FabricationPlant Evolution Step *********//
				if(!fReactor[i]->IsFuelFixed())
					if(step - FabricationCycleTime >= fAbsoluteTime && step - FabricationCycleTime <= t && step < ReactorShutDownTime)
					{						// Set End of reactor cycle
						pair< map<cSecond, int>::iterator, bool > IResult;
						IResult = fTimeStep.insert( pair<cSecond ,int>(step - FabricationCycleTime,16) );
						if( !IResult.second ) IResult.first->second  |= 16;
					}

				if(step > fAbsoluteTime && step <= t && step < ReactorShutDownTime)
				{						// Set End of reactor cycle
					pair< map<cSecond, int>::iterator, bool > IResult = fTimeStep.insert( pair<cSecond ,int>(step,4) );
					if( !IResult.second ) IResult.first->second  |= 4;
				}

				//********* End/Start Of Reactor Cycle Step *********//
				map< cSecond, int >::iterator TV_it; // the time vector iterator
				//********* BackEnd fuel Cycle *********//
				for(TV_it = BackEndTimePath.begin(); TV_it != BackEndTimePath.end(); TV_it++)	// Loop on the BackEnd fuel Cycle Time path
				{
					// Test if each step of the Fuel Cycle BackEnd is after the actual time (AbsolutreTime) and before the end of the evolution (t)
					if( step + (*TV_it).first >= fAbsoluteTime && step + (*TV_it).first <= t)
					{
						pair< map<cSecond, int>::iterator, bool > IResult;
						IResult = fTimeStep.insert( pair<cSecond ,int>(step + (*TV_it).first, (*TV_it).second) );
						if( !IResult.second )
							IResult.first->second |= (*TV_it).second;
					}
				}



				step += ReactorCycleTime;
			}
			while(step <= t && step <= ReactorShutDownTime );
		}
		else
		{
			cout	       << "!!Warning!! !!!Reactor!!! Be carefull a reactor cycletime is set to 0 second....\"\n" << endl;
			GetLog()->fLog << "!!Warning!! !!!Reactor!!! Be carefull a reactor cycletime is set to 0 second....\"\n" << endl;
		}

	}




	//****** Print the Time Index ******//
	ofstream TimeStepfile("CLASS_TimeStep", ios_base::app);		// Open the File

	if(!TimeStepfile)
	{
		cout		<< "!!Warning!! !!!CLASS!!! Can't open \" CLASS_TimeStep \"\n" << endl;
		GetLog()->fLog 	<< "!!Warning!! !!!CLASS!!! Can't open \" CLASS_TimeStep \"\n" << endl;
	}
	map<cSecond ,int >::iterator it;
	for( it = fTimeStep.begin(); it != fTimeStep.end(); it++)
		TimeStepfile << (*it).first << " " << (*it).second << endl;

}
//________________________________________________________________________
void Scenario::OldBuildTimeVector(cSecond t)
{
	fTimeStep.clear();
	fTimeStep.insert( pair<cSecond ,int>(t,1) );
	//********* Printing Step *********//
	{
		cSecond step = 0;
		if(fAbsoluteTime == fStartingTime)
			step = fStartingTime;
		else
		{
			step = fAbsoluteTime;
		}
		if(step >= fAbsoluteTime )
			fTimeStep.insert( pair<cSecond ,int>(step,1) );
		step += fPrintStep;
		do
		{

			if(step >= fAbsoluteTime )
				fTimeStep.insert( pair<cSecond ,int>(step,1) );
			step += fPrintStep;
		}
		while( step < t );
	}

	for(int i = 0; i < (int)fReactor.size();i++)
	{
		double step = fReactor[i]->GetCreationTime();
		double coolingstep = fReactor[i]->GetOutBackEndFacility()->GetCycleTime();
		double fabricationstep = 0;

		if(!fReactor[i]->IsFuelFixed())
			fabricationstep = fReactor[i]->GetFabricationPlant()->GetCycleTime();


		//********* Reactor Evolution Step *********//
		// set destruction of a reactor
		if( (fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime() > fAbsoluteTime) &&
		   (fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime() < t) )
		{
			//********* Reactor Shutdown *********//
			pair< map<cSecond, int>::iterator, bool > IResult  = fTimeStep.insert( pair<cSecond ,int>(fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime(),2) );
			if( !IResult.second ) IResult.first->second |= 2;

			//********* End of Cooling after reactor Shutdown *********//
			if(fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime()+coolingstep >= fAbsoluteTime
			   && fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime()+coolingstep <= t)
			{
				pair< map<cSecond, int>::iterator, bool > IResult = fTimeStep.insert( pair<cSecond ,int>(fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime()+coolingstep,8) );
				if( !IResult.second ) IResult.first->second |= 8;
			}
		}


		//********* Start Of Reactor First Cycle *********//
		if(step >= fAbsoluteTime &&  step <= t && step < fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime())
		{
			pair< map<cSecond, int>::iterator, bool > IResult = fTimeStep.insert( pair<cSecond ,int>(step,4) );
			if( !IResult.second ) IResult.first->second |= 4;
		}

		//********* FabricationPlant Evolution Step *********//
		if(!fReactor[i]->IsFuelFixed())
		{
			if(step > fAbsoluteTime && step - fabricationstep <= t && step < fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime() )
			{
				pair< map<cSecond, int>::iterator, bool > IResult = fTimeStep.insert( pair<cSecond ,int>(step -fabricationstep,16) );
				if( !IResult.second ) IResult.first->second  |= 16;
			}
			else if(step - fabricationstep < fStartingTime)
			{
				cout		<< "!!Warning!! !!!CLASS!!! Can't Build Fuel before Scenario's start\"\n" << endl;
				GetLog()->fLog 	<< "!!Warning!! !!!CLASS!!! Can't Build Fuel before Scenario's start\"\n" << endl;
				exit(1);
			}
		}


		//********* Reactor related Step *********//
		step += fReactor[i]->GetCycleTime();
		if (fReactor[i]->GetCycleTime() !=0)
		{
			do
			{


				//********* FabricationPlant Evolution Step *********//
				if(!fReactor[i]->IsFuelFixed())
					if(step > fAbsoluteTime && step - fabricationstep <= t && step < fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime())
					{						// Set End of reactor cycle
						pair< map<cSecond, int>::iterator, bool > IResult = fTimeStep.insert( pair<cSecond ,int>(step -fabricationstep,16) );
						if( !IResult.second ) IResult.first->second  |= 16;
					}


				//********* End/Start Of Reactor Cycle Step *********//
				if(step > fAbsoluteTime && step <= t && step < fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime())
				{						// Set End of reactor cycle
					pair< map<cSecond, int>::iterator, bool > IResult = fTimeStep.insert( pair<cSecond ,int>(step,4) );
					if( !IResult.second ) IResult.first->second  |= 4;
				}

				//********* End of Cooling Step *********//
				if(step >= fAbsoluteTime && step + coolingstep <= t)			// Set End of Cooling
				{
					pair< map<cSecond, int>::iterator, bool > IResult = fTimeStep.insert( pair<cSecond ,int>(step+coolingstep,8) );
					if( !IResult.second ) IResult.first->second |= 8;
				}
				step += fReactor[i]->GetCycleTime();
			}
			while(step <= t && step <= fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime() );
		}
	}


	//*** In Case of Evolution Restart ****//
	for(int i =0; i < (int)fPool.size(); i++)
	{
		//********* End of Cooling Step *********//
		for(int j = 0; j<(int)fPool[i]->GetIVArray().size(); j++ )// Set End of Cooling
		{
			if(fPool[i]->GetCoolingStartingTime()[j] +  fPool[i]->GetCycleTime() > fAbsoluteTime )
			{
				pair< map<cSecond, int>::iterator, bool > IResult;
				IResult = fTimeStep.insert( pair<cSecond ,int>(fPool[i]->GetCoolingStartingTime()[j] + fPool[i]->GetCycleTime(),8) );
				if( !IResult.second ) IResult.first->second |= 8;
			}
		}
	}



	//****** Print the Time Index ******//
	ofstream TimeStepfile("CLASS_TimeStep", ios_base::app);		// Open the File

	if(!TimeStepfile)
	{
		cout		<< "!!Warning!! !!!CLASS!!! Can't open \" CLASS_TimeStep \"\n" << endl;
		GetLog()->fLog 	<< "!!Warning!! !!!CLASS!!! Can't open \" CLASS_TimeStep \"\n" << endl;
	}
	map<cSecond ,int >::iterator it;
	for( it = fTimeStep.begin(); it != fTimeStep.end(); it++)
		TimeStepfile << (*it).first << " " << (*it).second << endl;

}



//________________________________________________________________________
//___________________________ Evolution Method ___________________________
//________________________________________________________________________

void Scenario::PoolEvolution()
{

#pragma omp parallel for
	for(int i = 0; i < (int) fPool.size();i++)
		fPool[i]->Evolution(fAbsoluteTime);

	for(int i = 0; i < (int) fPool.size();i++)
		fPool[i]->Dump();

}

void Scenario::StorageEvolution()
{

#pragma omp parallel for
	for(int i = 0; i < (int) fStorage.size();i++)
		fStorage[i]->Evolution(fAbsoluteTime);


}

void Scenario::FabricationPlantEvolution()
{

	//#pragma omp parallel for
	for(int i = 0; i < (int) fFabricationPlant.size();i++)
		fFabricationPlant[i]->Evolution(fAbsoluteTime);


}

//________________________________________________________________________
void Scenario::ReactorEvolution()
{

	fParcPower = 0;
#pragma omp parallel for
	for(int i = 0; i < (int)fReactor.size(); i++)
		fReactor[i]->Evolution(fAbsoluteTime);


	for(int i = 0; i < (int)fReactor.size(); i++)
		fReactor[i]->Dump();


}

//________________________________________________________________________
void Scenario::Evolution(double t)
{


	BuildTimeVector( (cSecond)t );

	if(fNewTtree )
	{
		OpenOutputTree();
		OutAttach();
		UpdateParc();
		fOutT->Fill();
	}

	map<cSecond ,int >::iterator it;
	for(it = fTimeStep.begin(); it != fTimeStep.end(); it++)
	{

		fAbsoluteTime = (*it).first;


		if( (*it).second & 1 || (*it).second & 4 || (*it).second & 8 || (*it).second & 16 )
			StorageEvolution();

		if( (*it).second & 1 || (*it).second & 2 || (*it).second & 4 || (*it).second & 8 || (*it).second & 16 )
			PoolEvolution();

		if( (*it).second & 1 || (*it).second & 2 || (*it).second & 4 || (*it).second & 16 )
			FabricationPlantEvolution();

		if( (*it).second & 1 || (*it).second & 2 || (*it).second & 4 )
			ReactorEvolution();

		if( (*it).second & 1 || it == fTimeStep.begin() )
		{
#pragma omp single
			{
				UpdateParc();
				fOutT->Fill();
				ProgressPrintout( (cSecond)t);
			}
		}

	}
	cout << endl;


}

void Scenario::ProgressPrintout(cSecond t)
{

	double Time = (fAbsoluteTime-fStartingTime)/3600/24/365.25 ;
	double Total = (t-fStartingTime)/3600/24/365.25;

	cout << "                                                                                                       " << flush ;
	cout << "\r[";
	for(int i = 0; i < (int)(Time/Total*100.0); i++)
		cout << "|";
	for(int i = 100; i >= (int)(Time/Total*100.0); i--)
		cout << "-";
	cout << "] ";

	cout << " Processed ";
	if (Time < 10) cout << " ";
	if (Time < 100) cout << " ";
	cout << (int)Time << " / " << (int)Total << " Years \r" << flush;


	GetLog()->fLog 	<< "Proccessed";
	if (Time < 10) GetLog()->fLog << " ";
	if (Time < 100) GetLog()->fLog << " ";
	GetLog()->fLog << (int)Time << " / " << (int)Total << " Years \r" << endl;

}

//________________________________________________________________________
//______________________________ Out Method ______________________________
//________________________________________________________________________
void Scenario::UpdateParc()
{

	ResetQuantity();

	for(int i =0; i < (int)fFabricationPlant.size(); i++)
		fFuelFabrication += fFabricationPlant[i]->GetInsideIV();

	for(int i = 0; i < (int) fPool.size();i++)
		fTotalCooling += fPool[i]->GetInsideIV();

	for(int i = 0; i < (int)fStorage.size(); i++)
		fTotalStorage += fStorage[i]->GetInsideIV();

	for(int i = 0; i < (int)fReactor.size(); i++)
		fTotalInReactor += fReactor[i]->GetIVReactor();

	fIVTotal = fWaste + fTotalStorage + fTotalCooling + fFuelFabrication + fTotalInReactor;
	fIVInCycleTotal = fTotalStorage + fTotalCooling + fFuelFabrication + fTotalInReactor;


}



void Scenario::ResetQuantity()
{


	fTotalInReactor.Clear();
	fTotalStorage.Clear();
	fTotalCooling.Clear();
	fFuelFabrication.Clear();
	fIVInCycleTotal.Clear();
	fIVTotal.Clear();

}

//________________________________________________________________________
void Scenario::OpenOutputTree()
{


	cout << "Opening OutPut File ...\t";
	GetLog()->fLog << "Opening : " << fOutputFileName << " ...\t";
	fOutFile = new TFile(fOutputFileName.c_str(),"UPDATE");

	if(!fOutFile)
	{
		cout << "\nCould not open " << fOutputFileName <<endl;
		GetLog()->fLog << "\nCould not open " << fOutputFileName <<endl;
		exit(-1);
	}
	cout << "\t ...OK!" << endl;


	fOutT = new TTree(fOutputTreeName.c_str(), "Data Tree");
	cout << "Creating Data Tree ...\t";
	GetLog()->fLog << "Creating Data Tree ...\t";
	if(!fOutT)
	{
		cout << "\nCould not create Data Tree in " << fOutputFileName << endl;
		GetLog()->fLog << "\nCould not create Data Tree in " << fOutputFileName << endl;
		exit(-1);
	}
	fNewTtree = false;
	cout << "\t ...OK!" << endl;
	GetLog()->fLog <<  "\t ...OK!" << endl;

}
void Scenario::CloseOutputTree()
{


	fOutFile->ls();
	cout << "Writing outTree " << fOutputFileName << endl;
	GetLog()->fLog << "Writing outTree " << fOutputFileName << endl;
	fOutFile->Write();

	if(fOutFile->IsOpen()) {
		cout << "Deleting outTree : " << endl;
		GetLog()->fLog << "Deleting outTree : " << endl;
		delete fOutT;
		cout << "Closing file : " << fOutputFileName <<endl;
		GetLog()->fLog << "Closing file : " << fOutputFileName <<endl;
		fOutFile-> Close();
		delete fOutFile;
	} else {
		cout << "File was not opened " << fOutputFileName << endl;
		GetLog()->fLog << "File was not opened " << fOutputFileName << endl;
		exit(-1);
	}
}
//________________________________________________________________________
void Scenario::OutAttach()
{

	ResetQuantity();
	//Branch Absolut Time
	fOutT->Branch("AbsTime",&fAbsoluteTime,"AbsoluteTime/L");
	//Branch The Power installed in the Parc
	fOutT->Branch("ParcPower",&fParcPower,"ParcPower/D");



	// Branch the Sum IV


	fOutT->Branch("STOCK.", "IsotopicVector", &fTotalStorage);
	fOutT->Branch("FUELFABRICATION.", "IsotopicVector", &fFuelFabrication);
	fOutT->Branch("COOLING.", "IsotopicVector", &fTotalCooling);
	fOutT->Branch("REACTOR.", "IsotopicVector", &fTotalInReactor);
	fOutT->Branch("INCYCLE.", "IsotopicVector", &fIVInCycleTotal);
	fOutT->Branch("TOTAL.", "IsotopicVector", &fIVTotal);

	fOutT->Branch("GOD.", "IsotopicVector", &fGod);
	fOutT->Branch("WASTE.", "IsotopicVector", &fWaste);

	// Branch the separate object

	for(int i = 0; i < (int)fStorage.size(); i++)
		fOutT->Branch(fStorage[i]->GetName(), "Storage", &fStorage[i]);

	for(int i = 0; i < (int)fPool.size(); i++)

		fOutT->Branch(fPool[i]->GetName(), "Pool", &fPool[i]);

	for(int i = 0; i < (int)fReactor.size(); i++)

		fOutT->Branch(fReactor[i]->GetName(), "Reactor", &fReactor[i]);

	for(int i = 0; i < (int)fFabricationPlant.size(); i++)
		fOutT->Branch(fFabricationPlant[i]->GetName(), "FabricationPlant", &fFabricationPlant[i]);


}

//________________________________________________________________________
void Scenario::Write()
{




}

//________________________________________________________________________
void Scenario::Print()
{

	for(int i = 0; i < (int) fPool.size();i++)
	{
		cout << "!!!!!!!!!STEP : " << fAbsoluteTime/(int)(3600*24*365.25) << endl;
		cout << "Pool : " << endl;
		cout << "Cooling ";
		cout << fPool[i]->GetIVArray().size()<< endl;
	}

	for(int i = 0; i < (int)fReactor.size(); i++)
	{
		cout << "Reactor" << endl;
		fReactor[i]->GetIVReactor().Print();
	}

}
