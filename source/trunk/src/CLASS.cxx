#include "CLASS.hxx"

#include "Storage.hxx"
#include "Reactor.hxx"
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
CLASS::CLASS()
{
	
	
	fNewTtree = true;
	fPrintStep = (cSecond)(3600*24*365.25);  // One Step per Year
	fAbsoluteTime = 0;
	fStartingTime = 0;
	
	fStockManagement = true;
	
	fOutputFileName = "CLASS_Default.root";
	fOutputTreeName = "Data";
	
	SetLog(new LogFile("CLASS.log"));
	fParcPower = 0;
	
	
		// Warning
	
	cout	<< "!!INFO!! !!!CLASS!!! A Parc has been define :" << endl;
	cout	<< "\t Print set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	cout	<< "\t Absolute Time set at " << (double)(fAbsoluteTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t StockManagement set at : true" << endl;
	cout	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	cout	<< "\t Log will be in " << GetLog()->GetLogFileName() << endl << endl;
	
	GetLog()->fLog 	<< "!!INFO!! !!!CLASS!!! Parc has been define :" << endl;
	GetLog()->fLog	<< "\t Print  set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t StockManagement set at : true" << endl;
	GetLog()->fLog	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	GetLog()->fLog	<< "\t Log will be in " << GetLog()->GetLogFileName() << endl << endl;
	
	
	
}
	//________________________________________________________________________
CLASS::CLASS(LogFile* Log)
{
	
	
	fNewTtree = true;
	fPrintStep = (cSecond)(3600*24*365.25);  // One Step per Year
	fAbsoluteTime = 0;
	fStartingTime = 0;
	
	fStockManagement = true;
	
	fOutputFileName = "CLASS_Default.root";
	fOutputTreeName = "Data";
	
	SetLog(Log);
	fParcPower = 0;
	
	
		// Warning
	
	cout	<< "!!INFO!! !!!CLASS!!! A Parc has been define :" << endl;
	cout	<< "\t Print set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	cout	<< "\t Absolute Time set at " << (double)(fAbsoluteTime/3600/24/365.25) << " year" << endl;
	cout	<< "\t StockManagement set at : true" << endl;
	cout	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	cout	<< "\t Log will be in " << GetLog()->GetLogFileName() << endl << endl;
	
	GetLog()->fLog 	<< "!!INFO!! !!!CLASS!!! Parc has been define :" << endl;
	GetLog()->fLog	<< "\t Print  set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t StockManagement set at : true" << endl;
	GetLog()->fLog	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	GetLog()->fLog	<< "\t Log will be in " << GetLog()->GetLogFileName() << endl << endl;
	
	
	
}
	//________________________________________________________________________
CLASS::CLASS(double abstime)
{
	
	fNewTtree = true;
	fPrintStep = (cSecond)(3600*24*365.25);  // One Step per Year
	fAbsoluteTime = (cSecond)abstime;
	fStartingTime = fAbsoluteTime;
	
	fStockManagement = true;
	
	fOutputFileName = "CLASS_Default.root";
	fOutputTreeName = "Data";
	
	SetLog(new LogFile("CLASS.log"));
	fParcPower = 0;
	
	
	
	
		// Warning
	
	cout	<< "!!INFO!! !!!CLASS!!! A Parc has been define :" << endl;
	cout	<< "\t Print set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	cout	<< "\t StockManagement set at : true" << endl;
	cout	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	cout	<< "\t Log will be in " << GetLog()->GetLogFileName() << endl;
	
	GetLog()->fLog 	<< "!!INFO!! !!!CLASS!!! Parc has been define :" << endl;
	GetLog()->fLog	<< "\t Print  set at : " << (double)(fPrintStep/3600/24/365.25) << " year" << endl;
	GetLog()->fLog	<< "\t StockManagement set at : true" << endl;
	GetLog()->fLog	<< "\t OutPut will be in \"" << fOutputFileName << "\" File and \"" << fOutputTreeName << "\" TTree" << endl;
	GetLog()->fLog	<< "\t Log will be in " << GetLog()->GetLogFileName() << endl << endl;
	
	
	
}


	//________________________________________________________________________
CLASS::~CLASS()
{
	
#pragma omp single
	{CloseOutputTree();}
	
}

	//________________________________________________________________________
void CLASS::AddPool(Pool* Pool)
{
	
	
	fPool.push_back(Pool);
	fPool.back()->SetParc(this);
	fPool.back()->SetDecayDataBase( (*this).GetDecayDataBase() );
	fPool.back()->SetLog(GetLog());
	fPool.back()->SetId((int)fPool.size()-1);
	
	
	if(fNewTtree == false)
	{
		string Pool_name = "Pool";
		Pool_name += dtoa(fPool.back()->GetId());
		Pool_name += ".";		
		fOutT->Branch(Pool_name.c_str(), "Pool", &fPool.back());
	}
	
}

	//________________________________________________________________________
void CLASS::AddReactor(Reactor* reactor)
{
	
	fReactor.push_back(reactor);
	fReactor.back()->SetParc(this);
	fReactor.back()->SetLog(GetLog());
	fReactor.back()->SetId((int)fReactor.size()-1);
	if(fReactor.back()->IsFuelFixed() == false)
		fReactor.back()->GetFabricationPlant()->AddReactor( (int)fReactor.size()-1,fReactor.back()->GetCreationTime() );

	
	
	if(fNewTtree == false)
	{
		string Reactor_name = "Reactor";
		Reactor_name += dtoa(fReactor.back()->GetId());
		Reactor_name += ".";
		fOutT->Branch(Reactor_name.c_str(), "Reactor", &fReactor.back());
	}

	
}

	//________________________________________________________________________
void CLASS::AddStorage(Storage* storage)
{
	
	fStorage.push_back(storage);
	fStorage.back()->SetParc(this);
	fStorage.back()->SetDecayDataBase( (*this).GetDecayDataBase() );
	fStorage.back()->SetLog(GetLog());
	fStorage.back()->SetId((int)fStorage.size()-1);

	if(fNewTtree == false)
	{
		string Storage_name = "Storage";
		Storage_name += dtoa(fStorage.back()->GetId());
		Storage_name += ".";
		fOutT->Branch(Storage_name.c_str(), "Storage", &fStorage.back());
	}

	
}
	//________________________________________________________________________
void CLASS::AddFabricationPlant(FabricationPlant* fabricationplant)
{
	
	fFabricationPlant.push_back(fabricationplant);
	fFabricationPlant.back()->SetParc(this);
	fFabricationPlant.back()->SetDecayDataBase( (*this).GetDecayDataBase() );
	fFabricationPlant.back()->SetLog(GetLog());
	fFabricationPlant.back()->SetId((int)fStorage.size()-1);
	
	if(fNewTtree == false)
	{
		string FabricationPlant_name = "FabricationPlant";
		FabricationPlant_name += dtoa(fFabricationPlant.back()->GetId());
		FabricationPlant_name += ".";
		fOutT->Branch(FabricationPlant_name.c_str(), "FabricationPlant", &fFabricationPlant.back());
	}

	
}

	//________________________________________________________________________
void CLASS::BuildTimeVector(cSecond t)
{
	
	fTimeStep.clear();
	fTimeStep.insert( pair<double ,int>(t,1) );
	
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
		double coolingstep = fReactor[i]->GetAssociedPool()->GetCoolingTime();
		double fabricationstep = 0;
		
		if(fReactor[i]->IsFuelFixed() == false)
			fabricationstep = fReactor[i]->GetFabricationPlant()->GetCycleTime();
		
		
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
				GetLog()->fLog 	<< "!!Warning!! !!!CLASS!!! Can't Build Fuel before Scenario's start\"\n" << endl;
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
	for(int i =0; i < (int)fPool.size(); i++)
	{
		
		
			//********* End of Cooling Step *********//
		for(int j = 0; j<(int)fPool[i]->GetIVCooling().size(); j++ )// Set End of Cooling
		{
			if(fPool[i]->GetCoolingStartingTime()[j] +  fPool[i]->GetCoolingTime() > fAbsoluteTime )
			{
				pair< map<cSecond, int>::iterator, bool > IResult;
				IResult = fTimeStep.insert( pair<cSecond ,int>(fPool[i]->GetCoolingStartingTime()[j] +  fPool[i]->GetCoolingTime(),8) );
				if( IResult.second == false ) IResult.first->second |= 8;
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
		TimeStepfile << (double)((*it).first/3600/24./365.25) << " " << (*it).second << endl;
	
	
}



	//________________________________________________________________________
	//___________________________ Evolution Method ___________________________
	//________________________________________________________________________

void CLASS::PoolEvolution()
{
	
	for(int i = 0; i < (int) fPool.size();i++)
		fPool[i]->Evolution(fAbsoluteTime);
	
	for(int i = 0; i < (int) fPool.size();i++)
		fPool[i]->Dump();
	
}

void CLASS::StorageEvolution()
{
	
	for(int i = 0; i < (int) fStorage.size();i++)
		fStorage[i]->Evolution(fAbsoluteTime);
	
	
}

void CLASS::FabricationPlantEvolution()
{
	
	for(int i = 0; i < (int) fFabricationPlant.size();i++)
		fFabricationPlant[i]->Evolution(fAbsoluteTime);
	
	
}

	//________________________________________________________________________
void CLASS::ReactorEvolution()
{
	
	fParcPower = 0;
#pragma omp parallel for
	for(int i = 0; i < (int)fReactor.size(); i++)
		fReactor[i]->Evolution(fAbsoluteTime);
	
	
	for(int i = 0; i < (int)fReactor.size(); i++)
		fReactor[i]->Dump();
	
	
}

	//________________________________________________________________________
void CLASS::Evolution(double t)
{
	
	
	BuildTimeVector( (cSecond)t );

	if(fNewTtree == true)
	{
		OpenOutputTree();
		OutAttach();
		fOutT->Fill();
	}

	map<cSecond ,int >::iterator it;
	for(it = fTimeStep.begin(); it != fTimeStep.end(); it++)
	{
		
		fAbsoluteTime = (*it).first;
		if( (*it).second & 2 || (*it).second & 1 )
		{
			if( (*it).second & 1 )
				StorageEvolution();
			PoolEvolution();
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
			PoolEvolution();
			ReactorEvolution();
			
			(*it).second ^= 4;
			if((*it).second & 8 )
				(*it).second ^= 8;
		}
		
		if( (*it).second & 16 )
		{
			StorageEvolution();
			PoolEvolution();
			FabricationPlantEvolution();
			
			(*it).second ^= 16;
			if((*it).second & 8 )
				(*it).second ^= 8;
		}
		
		
		if( (*it).second & 8 )
		{
			StorageEvolution();
			PoolEvolution();
			
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
	
	
}

void CLASS::ProgressPrintout(cSecond t)
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
	
}

	//________________________________________________________________________
	//______________________________ Out Method ______________________________
	//________________________________________________________________________
void CLASS::UpdateParc()
{
	
	ResetQuantity();
	
	for (int i =0; i < (int)fFabricationPlant.size(); i++)
		fFuelFabrication += fFabricationPlant[i]->GetFullFabrication();
	
	for(int i = 0; i < (int) fPool.size();i++)
		fTotalCooling += fPool[i]->GetFullCooling();
	
	for(int i = 0; i < (int)fStorage.size(); i++)
		fTotalStorage += fStorage[i]->GetFullStock();
	
	for(int i = 0; i < (int)fReactor.size(); i++)
		fTotalInReactor += fReactor[i]->GetIVReactor();
	
	fIVTotal = fWaste + fTotalStorage + fTotalCooling + fFuelFabrication + fTotalInReactor;
	fIVInCycleTotal = fTotalStorage + fTotalCooling + fFuelFabrication + fTotalInReactor;
	
	
}



void CLASS::ResetQuantity()
{
	
	
	fTotalInReactor.Clear();
	fTotalStorage.Clear();
	fTotalCooling.Clear();
	fFuelFabrication.Clear();
	fIVInCycleTotal.Clear();
	fIVTotal.Clear();
	
}

	//________________________________________________________________________
void CLASS::OpenOutputTree()
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
void CLASS::CloseOutputTree()
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
void CLASS::OutAttach()
{
	
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
	
		// Branch the separate object
	
	for(int i = 0; i < (int)fStorage.size(); i++)
	{
		string Storage_name = "Storage";
		Storage_name += dtoa(i);
		Storage_name += ".";
		
		fOutT->Branch(Storage_name.c_str(), "Storage", &fStorage[i]);
	}
	
	for(int i = 0; i < (int)fPool.size(); i++)
	{
		string TF_name = "Pool";
		TF_name += dtoa(i);
		TF_name += ".";
		fOutT->Branch(TF_name.c_str(), "Pool", &fPool[i]);
	}
	
	for(int i = 0; i < (int)fReactor.size(); i++)
	{
		string R_name = "Reactor";
		R_name += dtoa(i);
		R_name += ".";
		fOutT->Branch(R_name.c_str(), "Reactor", &fReactor[i]);
	}
	for(int i = 0; i < (int)fFabricationPlant.size(); i++)
	{
		string FP_name = "FabricationPlant";
		FP_name += dtoa(i);
		FP_name += ".";
		fOutT->Branch(FP_name.c_str(), "FabricationPlant", &fFabricationPlant[i]);
	}
	
}

	//________________________________________________________________________
void CLASS::Write()
{
	
	
	
	
}

	//________________________________________________________________________
void CLASS::Print()
{
	
	for(int i = 0; i < (int) fPool.size();i++)
	{
		cout << "!!!!!!!!!STEP : " << fAbsoluteTime/(int)(3600*24*365.25) << endl;
		cout << "Pool : " << endl;
		cout << "Cooling ";
		cout << fPool[i]->GetIVCooling().size()<< endl;
	}
	
	for(int i = 0; i < (int)fReactor.size(); i++)
	{
		cout << "Reactor" << endl;
		fReactor[i]->GetIVReactor().Print();
	}
	
}
