#include "CLASS.hxx"

#include "Storage.hxx"
#include "Reactor.hxx"
#include "TreatmentFactory.hxx"
#include "FabricationPlant.hxx"
#include "Defines.hxx"
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
CLASS::CLASS()
{
DBGL;
	fPrintStep = (cSecond)(3600*24*365.25);  // One Step per Year
	fAbsoluteTime = 0;
	fStartingTime = 0;

	fStockManagement = true;	
	
	fOutputFileName = "CLASS_Default.root";
	fOutputTreeName = "Data";
	
	fOutputLogName = "CLASS.log";
	fLog = new LogFile(fOutputLogName.c_str());
	fParcPower = 0;


	// Warning
	
	cout	<< "!!Info!! !!!CLASS!!! A Parc has been define :" << endl;
	cout	<< "\t Print set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	cout	<< "\t Absolute Time set at " << (double)(fAbsoluteTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t StockManagement set at : true" << endl;
	cout	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	cout	<< "\t Log will be in " << fOutputLogName << endl << endl;
	
	fLog->fLog 	<< "!!Info!! !!!CLASS!!! Parc has been define :" << endl;
	fLog->fLog	<< "\t Print  set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "\t StockManagement set at : true" << endl;
	fLog->fLog	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	fLog->fLog	<< "\t Log will be in " << fOutputLogName << endl << endl;


DBGL;
}

//________________________________________________________________________
CLASS::CLASS(double abstime)
{
DBGL;
	fPrintStep = (cSecond)(3600*24*365.25);  // One Step per Year
	fAbsoluteTime = (cSecond)abstime;
	fStartingTime = fAbsoluteTime;
	
	fStockManagement = true;
	
	fOutputFileName = "CLASS_Default.root";
	fOutputTreeName = "Data";
	
	fOutputLogName = "CLASS.log" ;
	fLog = new LogFile(fOutputLogName.c_str());
	fParcPower = 0;

	
	
	
	// Warning
	
	cout	<< "!!Info!! !!!CLASS!!! A Parc has been define :" << endl;
	cout	<< "\t Print set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	cout	<< "\t StockManagement set at : true" << endl;
	cout	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	cout	<< "\t Log will be in " << fOutputLogName << endl;
	
	fLog->fLog 	<< "!!Info!! !!!CLASS!!! Parc has been define :" << endl;
	fLog->fLog	<< "\t Print  set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	fLog->fLog	<< "\t StockManagement set at : true" << endl;
	fLog->fLog	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	fLog->fLog	<< "\t Log will be in " << fOutputLogName << endl << endl;


DBGL;
}


//________________________________________________________________________
CLASS::~CLASS()
{
DBGL;
DBGL;
}

//________________________________________________________________________
void CLASS::AddTreatmentFactory(TreatmentFactory* treatmentfactory)
{
DBGL;
	
	fTreatmentFactory.push_back(treatmentfactory);
	fTreatmentFactory.back()->SetParc(this);
	fTreatmentFactory.back()->SetDecayDataBase( (*this).GetDecayDataBase() );
	
	fTreatmentFactory.back()->SetLog(fLog);
	fTreatmentFactory.back()->SetId((int)fTreatmentFactory.size()-1);
DBGL;
}

//________________________________________________________________________
void CLASS::AddReactor(Reactor* reactor)
{
DBGL;
	fReactor.push_back(reactor);
	fReactor.back()->SetParc(this);
	fReactor.back()->SetLog(fLog);
	fReactor.back()->SetId((int)fReactor.size()-1);
	if(fReactor.back()->IsFuelFixed() == false)
		fReactor.back()->GetFabricationPlant()->AddReactor( (int)fReactor.size()-1,fReactor.back()->GetCreationTime() );
DBGL;
}

//________________________________________________________________________
void CLASS::AddStorage(Storage* storage)
{
DBGL;
	fStorage.push_back(storage);
	fStorage.back()->SetParc(this);
	fStorage.back()->SetDecayDataBase( (*this).GetDecayDataBase() );
	fStorage.back()->SetLog(fLog);
	fStorage.back()->SetId((int)fStorage.size()-1);

DBGL;
}
//________________________________________________________________________
void CLASS::AddFabricationPlant(FabricationPlant* fabricationplant)
{
DBGL;
	fFabricationPlant.push_back(fabricationplant);
	fFabricationPlant.back()->SetParc(this);
	fFabricationPlant.back()->SetDecayDataBase( (*this).GetDecayDataBase() );
	fFabricationPlant.back()->SetLog(fLog);
	fFabricationPlant.back()->SetId((int)fStorage.size()-1);

DBGL;
}


//________________________________________________________________________
void CLASS::BuildTimeVector(cSecond t)
{
DBGL;
	
	fTimeStep.clear();
	fTimeStep.insert( pair<double ,int>(t,1) );

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
		double step = fReactor[i]->GetCreationTime();
		double coolingstep = fReactor[i]->GetAssociedTreatmentFactory()->GetCoolingTime();
		double fabricationstep = 0;

		if(fReactor[i]->IsFuelFixed() == false)
			fabricationstep = fReactor[i]->GetFabricationPlant()->GetFabricationTime();


//********* Reactor Evolution Step *********//
		// set destruction of a reactor
		if( (fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime() > fAbsoluteTime) &&
		    (fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime() < t) )
		{
			//********* Reactor Shutdown *********//
			pair< map<cSecond, int>::iterator, bool > IResult  = fTimeStep.insert( pair<cSecond ,int>(fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime(),2) );
			if( IResult.second == false ) IResult.first->second |= 2;
			
			//********* End of Cooling after reactor Shutdown *********//
			if(fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime()+coolingstep >= fAbsoluteTime 
			&& fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime()+coolingstep <= t)
			{
				pair< map<cSecond, int>::iterator, bool > IResult = fTimeStep.insert( pair<cSecond ,int>(fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime()+coolingstep,8) );
				if( IResult.second == false ) IResult.first->second |= 8;
			}
		}

		
		//********* Start Of Reactor First Cycle *********//
		if(step >= fAbsoluteTime &&  step <= t && step < fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime())		
		{
			pair< map<cSecond, int>::iterator, bool > IResult = fTimeStep.insert( pair<cSecond ,int>(step,4) );
			if( IResult.second == false ) IResult.first->second |= 4;
		}

		//********* FabricationPlant Evolution Step *********//
		if(fReactor[i]->IsFuelFixed() == false)
		{
			if(step > fAbsoluteTime && step - fabricationstep <= t && step < fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime() )
			{
				pair< map<cSecond, int>::iterator, bool > IResult = fTimeStep.insert( pair<cSecond ,int>(step -fabricationstep,16) );
				if( IResult.second == false ) IResult.first->second  |= 16;
			}
			else if(step - fabricationstep < fStartingTime)
			{
				cout		<< "!!Warning!! !!!CLASS!!! Can't Build Fuel before Scenario's start\"\n" << endl;
				fLog->fLog 	<< "!!Warning!! !!!CLASS!!! Can't Build Fuel before Scenario's start\"\n" << endl;
				exit(1);
			}
		}


//********* Reactor related Step *********//
		step += fReactor[i]->GetCycleTime();
		do 
		{


			//********* FabricationPlant Evolution Step *********//
			if(fReactor[i]->IsFuelFixed() == false)
				if(step > fAbsoluteTime && step - fabricationstep <= t && step < fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime())
				{						// Set End of reactor cycle
					pair< map<cSecond, int>::iterator, bool > IResult = fTimeStep.insert( pair<cSecond ,int>(step -fabricationstep,16) );
					if( IResult.second == false ) IResult.first->second  |= 16;
				}


			//********* End/Start Of Reactor Cycle Step *********//
			if(step > fAbsoluteTime && step <= t && step < fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime())
			{						// Set End of reactor cycle
				pair< map<cSecond, int>::iterator, bool > IResult = fTimeStep.insert( pair<cSecond ,int>(step,4) );
				if( IResult.second == false ) IResult.first->second  |= 4;
			}

			//********* End of Cooling Step *********//
			if(step >= fAbsoluteTime && step + coolingstep <= t)			// Set End of Cooling
			{
				pair< map<cSecond, int>::iterator, bool > IResult = fTimeStep.insert( pair<cSecond ,int>(step+coolingstep,8) );
				if( IResult.second == false ) IResult.first->second |= 8;
			}

			step += fReactor[i]->GetCycleTime();
		}
		while(step <= t && step <= fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime() );
	}


//*** In Case of Evolution Restart ****//
	for(int i =0; i < (int)fTreatmentFactory.size(); i++) 
	{
		
		
		//********* End of Cooling Step *********//
		for(int j = 0; j<(int)fTreatmentFactory[i]->GetIVCooling().size(); j++ )// Set End of Cooling
		{
			if(fTreatmentFactory[i]->GetCoolingStartingTime()[j] +  fTreatmentFactory[i]->GetCoolingTime() > fAbsoluteTime )
			{
				pair< map<cSecond, int>::iterator, bool > IResult;
				IResult = fTimeStep.insert( pair<cSecond ,int>(fTreatmentFactory[i]->GetCoolingStartingTime()[j] +  fTreatmentFactory[i]->GetCoolingTime(),8) );
				if( IResult.second == false ) IResult.first->second |= 8;
			}
		}
	}

	
	
//****** Print the Time Index ******//
	ofstream TimeStepfile("CLASS_TimeStep", ios_base::app);		// Open the File
	
	if(!TimeStepfile)
	{
		cout		<< "!!Warning!! !!!CLASS!!! Can't open \" CLASS_TimeStep \"\n" << endl;
		fLog->fLog 	<< "!!Warning!! !!!CLASS!!! Can't open \" CLASS_TimeStep \"\n" << endl;
	}
	map<cSecond ,int >::iterator it;
	for( it = fTimeStep.begin(); it != fTimeStep.end(); it++)
		TimeStepfile << (double)((*it).first/3600/24./365.25) << " " << (*it).second << endl;
	
	DBGL;	
}



//________________________________________________________________________
//___________________________ Evolution Method ___________________________
//________________________________________________________________________

void CLASS::TreatmentEvolution()
{
DBGL;
	for(int i = 0; i < (int) fTreatmentFactory.size();i++)
		fTreatmentFactory[i]->Evolution(fAbsoluteTime);
	
	for(int i = 0; i < (int) fTreatmentFactory.size();i++)
		fTreatmentFactory[i]->Dump();
DBGL;
}

void CLASS::StorageEvolution()
{
DBGL;
	for(int i = 0; i < (int) fStorage.size();i++)
		fStorage[i]->Evolution(fAbsoluteTime);

DBGL;
}

void CLASS::FabricationPlantEvolution()
{
DBGL;
	for(int i = 0; i < (int) fFabricationPlant.size();i++)
		fFabricationPlant[i]->Evolution(fAbsoluteTime);

DBGL;
}

//________________________________________________________________________
void CLASS::ReactorEvolution()
{
DBGL;
	fParcPower = 0;
#pragma omp parallel for 
	for(int i = 0; i < (int)fReactor.size(); i++)
		fReactor[i]->Evolution(fAbsoluteTime);
	
	
	for(int i = 0; i < (int)fReactor.size(); i++)
		fReactor[i]->Dump();

DBGL;
}

//________________________________________________________________________
void CLASS::Evolution(double t)
{
DBGL;

	BuildTimeVector( (cSecond)t );


	
	OpenOutputTree();

	OutAttach();

	fOutT->Fill();
	map<cSecond ,int >::iterator it;
	for(it = fTimeStep.begin(); it != fTimeStep.end(); it++)
	{

		fAbsoluteTime = (*it).first;
		if( (*it).second & 2 || (*it).second & 1 )
		{
			if( (*it).second & 1 )
			StorageEvolution();
			TreatmentEvolution();
			FabricationPlantEvolution();
			ReactorEvolution();
		
			if((*it).second & 2 )
				(*it).second ^= 2;
			if((*it).second & 4 )
				(*it).second ^= 4;
			if((*it).second & 8 )
				(*it).second ^= 8;
			if((*it).second & 16 )
				(*it).second ^= 16;
		}

		
		if( (*it).second & 4 )
		{
			StorageEvolution();
			TreatmentEvolution();
			ReactorEvolution();

			(*it).second ^= 4;
			if((*it).second & 8 )
				(*it).second ^= 8;
		}

		if( (*it).second & 16 )
		{
			StorageEvolution();
			TreatmentEvolution();
			FabricationPlantEvolution();
		
			(*it).second ^= 16;
			if((*it).second & 8 )
				(*it).second ^= 8;
		}


		if( (*it).second & 8 )
		{
			StorageEvolution();
			TreatmentEvolution();

			(*it).second ^= 8;
		}
	

	
		if( (*it).second & 1 || (*it).first == (*fTimeStep.begin()).first )
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
#pragma omp single
	{CloseOutputTree();}
DBGL;
}

void CLASS::ProgressPrintout(cSecond t)
{
DBGL;
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
DBGL;
}

//________________________________________________________________________
//______________________________ Out Method ______________________________
//________________________________________________________________________
void CLASS::UpdateParc()
{
DBGL;
	ResetQuantity();
	
	for (int i =0; i < (int)fFabricationPlant.size(); i++)
	{
		map<int, IsotopicVector > reactorNextStep = fFabricationPlant[i]->GetReactorFuturIncome();
		map<int, IsotopicVector >::iterator it;
		for ( it = reactorNextStep.begin(); it != reactorNextStep.end(); it++)
			fFuelFabrication += (*it).second;
	}

	for(int i = 0; i < (int) fTreatmentFactory.size();i++)
	{
		for(int j=0; j < (int)fTreatmentFactory[i]->GetIVCooling().size(); j++)
			fTotalCooling += fTreatmentFactory[i]->GetIVCooling()[j];
	}
	
	for(int i = 0; i < (int)fStorage.size(); i++)
	{
		fTotalStorage += fStorage[i]->GetFullStock();
	}
	for(int i = 0; i < (int)fReactor.size(); i++)
		fTotalInReactor += fReactor[i]->GetIVReactor();

	fIVTotal = fWaste + fTotalStorage + fTotalCooling + fFuelFabrication + fTotalInReactor;
	fIVInCycleTotal = fTotalStorage + fTotalCooling + fFuelFabrication + fTotalInReactor;
	
DBGL;	
}



void CLASS::ResetQuantity()
{
DBGL;

	fTotalInReactor.Clear();
	fTotalStorage.Clear();
	fTotalCooling.Clear();
	fFuelFabrication.Clear();
	fIVInCycleTotal.Clear();
	fIVTotal.Clear();
DBGL;
}

//________________________________________________________________________
void CLASS::OpenOutputTree()
{
DBGL;

	cout << "Opening OutPut File ...\t";
	fLog->fLog << "Opening : " << fOutputFileName << " ...\t";
	fOutFile = new TFile(fOutputFileName.c_str(),"UPDATE");

	if(!fOutFile)
	{
		cout << "\nCould not open " << fOutputFileName <<endl;
		fLog->fLog << "\nCould not open " << fOutputFileName <<endl;
		exit(-1);
	}
	cout << "\t ...OK!" << endl;

	fOutT = new TTree(fOutputTreeName.c_str(), "Data Tree");
	cout << "Creating Data Tree ...\t";
	fLog->fLog << "Creating Data Tree ...\t";
	if(!fOutT)
	{
		cout << "\nCould not create Data Tree in " << fOutputFileName << endl;
		fLog->fLog << "\nCould not create Data Tree in " << fOutputFileName << endl;
		exit(-1);
	}
	cout << "\t ...OK!" << endl;
	fLog->fLog <<  "\t ...OK!" << endl;

}
void CLASS::CloseOutputTree()
{
DBGL;

	fOutFile->ls();
	cout << "Writing outTree " << fOutputFileName << endl;
	fLog->fLog << "Writing outTree " << fOutputFileName << endl;
	fOutFile->Write();

	if(fOutFile->IsOpen()) {
		cout << "Deleting outTree : " << endl;
		fLog->fLog << "Deleting outTree : " << endl;
		delete fOutT;
		cout << "Closing file : " << fOutputFileName <<endl;
		fLog->fLog << "Closing file : " << fOutputFileName <<endl;
		fOutFile-> Close();
		delete fOutFile;
	} else {
		cout << "File was not opened " << fOutputFileName << endl;
		fLog->fLog << "File was not opened " << fOutputFileName << endl;
		exit(-1);
	}
}
//________________________________________________________________________
void CLASS::OutAttach()
{
DBGL;
	ResetQuantity();
	//Branch Absolut Time
	fOutT->Branch("AbsTime",&fAbsoluteTime,"AbsoluteTime/l");
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
	
	for(int i = 0; i < (int)fStorage.size(); i++)
	{
		string Storage_name = "Storage";
		Storage_name += dtoa(i);
		Storage_name += ".";
		
		fOutT->Branch(Storage_name.c_str(), "Storage", &fStorage[i]);
	}

	for(int i = 0; i < (int)fTreatmentFactory.size(); i++)
	{
		string TF_name = "TreatmentFactory";
		TF_name += dtoa(i);
		TF_name += ".";		
		fOutT->Branch(TF_name.c_str(), "TreatmentFactory", &fTreatmentFactory[i]);
	}
	
	for(int i = 0; i < (int)fReactor.size(); i++)
	{
		string R_name = "Reactor";
		R_name += dtoa(i);
		R_name += ".";		
		fOutT->Branch(R_name.c_str(), "Reactor", &fReactor[i]);
	}
DBGL;
}

//________________________________________________________________________
void CLASS::Write()
{	
DBGL;


DBGL;
}

//________________________________________________________________________
void CLASS::Print()
{
DBGL;
	for(int i = 0; i < (int) fTreatmentFactory.size();i++)
	{
		cout << "!!!!!!!!!STEP : " << fAbsoluteTime/(int)(3600*24*365.25) << endl;
		cout << "TreatmentFactory : " << endl;
		cout << "Cooling ";
		cout << fTreatmentFactory[i]->GetIVCooling().size()<< endl;
	}

	for(int i = 0; i < (int)fReactor.size(); i++)
	{	
		cout << "Reactor" << endl;
		fReactor[i]->GetIVReactor().Print();
	}
DBGL;
}
