#include "CLASS.hxx"
#include "Reactor.hxx"
#include "TreatmentFactory.hxx"
#include "Defines.hxx"
#include "IsotopicVector.hxx"

#include <ctime>
#include "time.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <TFile.h>
#include <TTree.h>


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

//________________________________________________________________________
CLASS::CLASS()
{
	DBGL;
	fPrintStep = (long int)(3600*24*365.4);  // One Step per Month
	fAbsoluteTime = 0;
	fTreatmentFactoryLastIndex = 0;
	fReactorLastIndex = 0;
	fStockManagement = true;
	fOutputName = "CLASS_Default.root";
	DBGL;
}

//________________________________________________________________________
CLASS::CLASS(long int abstime)
{
	DBGL;
	fPrintStep = (long int)(3600*24*365.4);  // One Step per Month
	fAbsoluteTime = abstime;
	fTreatmentFactoryLastIndex = 0;
	fReactorLastIndex = 0;
	fStockManagement = true;
	fOutputName = "CLASS_Default.root";
	DBGL;
}

//________________________________________________________________________
CLASS::CLASS(string name, long int abstime)
{
	DBGL;
	fPrintStep = (long int)(3600*24*365.4);  // One Step per Month
	fAbsoluteTime = abstime;
	fTreatmentFactoryLastIndex = 0;
	fReactorLastIndex = 0;
	fStockManagement = true;
	fOutputName = name;
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
	fTreatmentFactory.back()->SetStockManagement(fStockManagement);
	fTreatmentFactoryLastIndex++;
	fTreatmentFactoryIndex.push_back(fTreatmentFactoryLastIndex);
	DBGL;
}

//________________________________________________________________________
void CLASS::AddReactor(Reactor* reactor)
{
	DBGL;
	fReactor.push_back(reactor);
	fReactor.back()->SetParc(this);
	fReactorLastIndex++;
	fReactorIndex.push_back(fReactorLastIndex);
	DBGL;
}

//________________________________________________________________________
void CLASS::RemoveReactor()
{
	DBGL;
	for(int i = (int)fReactor.size()-1; i >= 0; i--)
	{
		if(fAbsoluteTime == fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime())
		{
			fReactor[i]->GetAssociedTreatmentFactory()->AddIVCooling(fReactor[i]->GetIVReactor());
			delete (fReactor[i]);
			fReactor.erase(fReactor.begin()+i);
			fReactorIndex.erase(fReactorIndex.begin()+i);
		}
	}
	DBGL;
}

//________________________________________________________________________
//____________________________ Building Method ___________________________
//________________________________________________________________________

IsotopicVector CLASS::BuildIsotopicVector(IsotopicVector isotopicvector)
{
	DBGL;
	IsotopicVector BuildIV;
	
	if (fStockManagement == false )
		return BuildIV;

	vector< map<pair<int,int>, double> > Zcontribution;
	vector<int>  Zlist;
	int k;
#pragma omp parallel
	{
		IsotopicVector BuildIVtmp;
		IsotopicVector BuildIVz;
		IsotopicVector BuildIVztmp;
		int z = 0;
		
		map<pair<int,int>, double> contribution;
		#pragma omp for nowait
		for(k = 0; k < (int)isotopicvector.GetAtomicSpecies().size(); k++ ) // Loop on the Atomic Species
		{
			z = isotopicvector.GetAtomicSpecies()[k]; 		// Get the Atomic Species Number

			IsotopicVector IVAtomicCompo = isotopicvector.GetAtomicComposition(z);
			double delta = 1e-3*Norme(IVAtomicCompo);
			vector<IsotopicVector>	IVStock;
			vector<double>		IVMaxQuantity;
			vector<int>		FactoryIndex;
			vector<int>		StockIndex;
		
			for(int i=0; i< (int)fTreatmentFactory.size(); i++ )
			{
				for(int j = 0; j < (int)fTreatmentFactory[i]->GetIVStock().size(); j++ )
				{
					FactoryIndex.push_back(i);
					StockIndex.push_back(j);
					double Quantitytmp = Norme(fTreatmentFactory[i]->GetIVStock()[j].GetAtomicComposition(z));
					IVMaxQuantity.push_back(Quantitytmp);
					IVStock.push_back( fTreatmentFactory[i]->GetIVStock()[j].GetAtomicComposition(z) / Quantitytmp );
				}
			}
		
		
			vector<IsotopicVector>	IVStocktmp 	= IVStock;
			vector<int>		FactoryIndextmp	= FactoryIndex;
			vector<int>		StockIndextmp 	= StockIndex;

			// Choose the Contribution
//			IsotopicVector BuildIVz;
//			IsotopicVector BuildIVztmp;
			double DistanceN = Norme(IVAtomicCompo);	// Distance at the Nst iteration
			double DistanceN_1 = Norme(IVAtomicCompo);	// Distance at the N-1 st iteration

			while(IVStocktmp.size() !=0)
			{
				int IVIndex = IVStocktmp.size();
				double IVIndexQuantityneed = 0;
			
			
				while(IVIndex == (int)IVStocktmp.size()) IVIndex = (int)random(0,IVStocktmp.size());
		
				BuildIVztmp = BuildIVz + (IVStocktmp[IVIndex] * delta);
				DistanceN = Distance( BuildIVztmp, isotopicvector.GetAtomicComposition(z) );
			
				// Get the contribution of the IVStock IVIndex
				map< pair<int,int>, double >::iterator it = 
						contribution.find(pair<int,int>(FactoryIndextmp[IVIndex], StockIndextmp[IVIndex] ));
				if(it != contribution.end() )
					IVIndexQuantityneed = delta + (*it).second;
		
				if(DistanceN > DistanceN_1 || IVIndexQuantityneed > IVMaxQuantity[IVIndex] )
				{
					IVStocktmp.erase(IVStocktmp.begin()+IVIndex);
					FactoryIndextmp.erase(FactoryIndextmp.begin()+IVIndex);
					StockIndextmp.erase(StockIndextmp.begin()+IVIndex);
				}
				else
				{
					// Add the delta contribution from the IVIndex Stock IV
					pair< map< pair<int,int>, double >::iterator, bool > IResult;
				
					IResult = contribution.insert( pair< pair< int ,int > , double > 
						( pair< int ,int >(FactoryIndextmp[IVIndex], StockIndextmp[IVIndex] ), delta) ) ;
					if(IResult.second == false)
						IResult.first->second += delta;

					// Reset value
					DistanceN_1 = DistanceN;
					BuildIVz = BuildIVztmp;
					IVStocktmp 	= IVStock;
					FactoryIndextmp = FactoryIndex;
					StockIndextmp 	= StockIndex;
				}
			}
			BuildIVtmp += BuildIVz;
//			BuildIV += BuildIVz;
//			Zcontribution.push_back( contribution );
//			Zlist.push_back(z);
		
		}
		#pragma omp critical(BuildIV)
		{
			BuildIV += BuildIVtmp;
			Zcontribution.push_back( contribution );
			Zlist.push_back(z);
		}
	
	}
//	#pragma omp critical(BuildIVStockupdate)
	{
	for(int k = 0; k < (int)Zlist.size(); k++ ) //Loop on the Atomic Species
	{
		// Take it and build the IsotopicVector
		map<pair<int,int>, double> contribution = Zcontribution[k];
		for(map<pair<int,int>, double>::iterator it = contribution.begin(); it != contribution.end(); it++)
		{
			int z = Zlist[k];
			double IVNorme = Norme(fTreatmentFactory[(*it).first.first]->GetIVStock()[(*it).first.second].GetAtomicComposition(z));
			IsotopicVector IV = fTreatmentFactory[(*it).first.first]->GetIVStock()[(*it).first.second].GetAtomicComposition(z);
			fTreatmentFactory[(*it).first.first]->TakeFromStock( IV / IVNorme * (*it).second , (*it).first.second );
		}
		
		
	 }
	}
	
	RemoveTotalStock(BuildIV);
	DBGL;
	return BuildIV;
}

//________________________________________________________________________
void CLASS::BuildTimeVector(long int t)
{
	DBGL;
	fTimeStep.clear();
	fTimeStep.insert( pair<long int ,int>(t,1) );

//********* Printing Step *********//
	{
		long int step = 0;
		if(step >= fAbsoluteTime	)
			fTimeStep.insert( pair<long int ,int>(step,1) );

		while( step < t )
		{
			step += fPrintStep;
			if(step >= fAbsoluteTime	)
				fTimeStep.insert( pair<long int ,int>(step,1) );
		}
	}

	for(int i = 0; i < (int)fReactor.size();i++)
	{
		long int step = fReactor[i]->GetCreationTime();
		long int coolingstep = fReactor[i]->GetAssociedTreatmentFactory()->GetCoolingTime();
		long int seprationstep = fReactor[i]->GetAssociedTreatmentFactory()->GetSeparationTime();

//********* Reactor Evolution Step *********//
		// set destruction of a reactor
		if( (fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime() > fAbsoluteTime) &&
		    (fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime()<t) )
		{
			pair< map<long int, int>::iterator, bool > IResult ;
			IResult = fTimeStep.insert( pair<long int ,int>(fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime(),2) );
			if( IResult.second == false ) IResult.first->second |= 2;
		}
		
		// Set End of reactor cycle
		if(step >= fAbsoluteTime &&  step <= t && step < fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime())		
		{
			pair< map<long int, int>::iterator, bool > IResult = fTimeStep.insert( pair<long int ,int>(step,4) );
			if( IResult.second == false ) IResult.first->second |= 4;
		}
		

//********* TF Evolution Step *********//
		while(step <= t && step <= fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime() )
		{
			step += fReactor[i]->GetCycleTime();
			if(step > fAbsoluteTime && step <= t && step < fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime())
			{						// Set End of reactor cycle
				pair< map<long int, int>::iterator, bool > IResult = fTimeStep.insert( pair<long int ,int>(step,4) );
				if( IResult.second == false ) IResult.first->second  |= 4;
			}

			if(step >= fAbsoluteTime && step + coolingstep <= t)			// Set End of Cooling
			{
				pair< map<long int, int>::iterator, bool > IResult = fTimeStep.insert( pair<long int ,int>(step+coolingstep,8) );
				if( IResult.second == false ) IResult.first->second |= 8;
			}

			if(step >= fAbsoluteTime && step + coolingstep + seprationstep <= t)	// Set end of Separation
			{
				pair< map<long int, int>::iterator, bool > IResult;
				IResult = fTimeStep.insert( pair<long int ,int>(step+coolingstep+seprationstep,16) );
				if( IResult.second == false ) IResult.first->second |= 16;
			}
		}
	}
	
//********* TF Evolution Step *********//
//*** In Case of Evolution Restart ****//
	for(int i =0; i < (int)fTreatmentFactory.size(); i++) 
	{
		
		for(int j = 0; j<(int)fTreatmentFactory[i]->GetIVCooling().size(); j++ )// Set End of Cooling
		{
			if(fTreatmentFactory[i]->GetCoolingStartingTime()[j] +  fTreatmentFactory[i]->GetCoolingTime() > fAbsoluteTime )
			{
				pair< map<long int, int>::iterator, bool > IResult;
				IResult = fTimeStep.insert( pair<long int ,int>(fTreatmentFactory[i]->GetCoolingStartingTime()[j] +  fTreatmentFactory[i]->GetCoolingTime(),8) );
				if( IResult.second == false ) IResult.first->second |= 8;
			}
		}
		for(int j = 0; j< (int)fTreatmentFactory[i]->GetIVSeparating().size(); j++ )// Set end of Separation
		{
			if(fTreatmentFactory[i]->GetSeparatingStartingTime()[j] +  fTreatmentFactory[i]->GetSeparationTime() > fAbsoluteTime )
			{
				pair< map<long int, int>::iterator, bool > IResult = fTimeStep.insert( pair<long int ,int>(fTreatmentFactory[i]->GetSeparatingStartingTime()[j] +  fTreatmentFactory[i]->GetSeparationTime(),16) );
				if( IResult.second == false ) IResult.first->second |= 16;
			}
		}
	}

	
	
//****** Print the Time Index ******//
	ofstream TimeStepfile("CLASS_TimeStep", ios_base::app);		// Open the File
	
	if(!TimeStepfile)
		cout << "!!Warning!! !!!CLASS!!! Can't open \" CLASS_TimeStep \"\n" << endl;

	for(map<long int ,int >::iterator it = fTimeStep.begin(); it != fTimeStep.end(); it++)
		TimeStepfile << (*it).first/365.4/3600/24 << " " << (*it).second << endl;
	
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
	DBGL;
}

//________________________________________________________________________
void CLASS::ReactorEvolution()
{
	DBGL;
	

 #pragma omp sections
	{
	#pragma omp section
		{TreatmentEvolution();}

	#pragma omp section
	{
		#pragma omp parallel for
			for(int i = 0; i < (int)fReactor.size(); i++)
				fReactor[i]->Evolution(fAbsoluteTime);
	}
	}
	

	for(int i = 0; i < (int) fTreatmentFactory.size();i++)
		fTreatmentFactory[i]->Dump();

	for(int i = 0; i < (int)fReactor.size(); i++)
		fReactor[i]->Dump();

	DBGL;
}

//________________________________________________________________________
void CLASS::Evolution(long int t)
{
	DBGL;
	OpenOutputTree();
	OutAttach();
	BuildTimeVector(t);

	int Start = time(NULL);
	for(map<long int ,int >::iterator it = fTimeStep.begin(); it != fTimeStep.end(); it++)
	{

		ResetQuantity();
		fAbsoluteTime = (*it).first;
		for(int i=0; i<4;i++)
			fIVTotal.addtest(fAbsoluteTime+i);
		if(fAbsoluteTime > t) return;
		
		if( (*it).second & 2 || (*it).second & 1 )
		{
			ReactorEvolution();
			if((*it).second & 2 )
				RemoveReactor();
			
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
			ReactorEvolution();
			if((*it).second & 4 )
				(*it).second ^= 4;
			if((*it).second & 8 )
				(*it).second ^= 8;
			if((*it).second & 16 )
				(*it).second ^= 16;
		}
		
		if( (*it).second & 8 || (*it).second & 16 )
		{
			TreatmentEvolution();
			for(int i = 0; i < (int) fTreatmentFactory.size();i++)
				fTreatmentFactory[i]->Dump();

			if((*it).second & 8 )
				(*it).second ^= 8;
			if((*it).second & 16 )
				(*it).second ^= 16;
		}
		

		
		if( (*it).second & 1 )
			#pragma omp critical(FillTTree)
			{
//				Print();
//				Write();
				fOutT->Fill();
				PorgressionPrintOut(Start, t);
			}
	
	}
	CloseOutputTree();
	DBGL;
}
void CLASS::PorgressionPrintOut(int starttime, long int t)
{
	DBGL;
	int TimeNow = time(NULL);
	int Spent = difftime(TimeNow, starttime);
	double Time = fAbsoluteTime/3600/24/365.4;
	double Total = t/3600/24/365.4;
	double Remain =  (Spent/Time * Total - Spent);
	cout << "                                                                                             " << flush ;
	cout << "\rProcessed " << Time << " / " << Total << " STEP "
	     << (Time/Total*100.0) << "%) ---  I still need : "
	     << (long int)Remain/60 << " min "<<  (long int)(Remain)%60 << " sec to finish ! \r"  << flush;
	DBGL;
}

//________________________________________________________________________
//______________________________ Out Method ______________________________
//________________________________________________________________________
void CLASS::Print()
{
	DBGL;
	for(int i = 0; i < (int) fTreatmentFactory.size();i++)
	{
		cout << "!!!!!!!!!STEP : " << fAbsoluteTime/(int)(3600*24*365.4) << endl;
		cout << "TreatmentFactory : " << endl;
		cout << "Cooling ";
		cout << fTreatmentFactory[i]->GetIVCooling().size()<< endl;
		cout << "Waste " << endl;
		fTreatmentFactory[i]->GetIVWaste().Print();
		cout << "Stock :" << endl;
		for(int j=0; j < (int)fTreatmentFactory[i]->GetIVStock().size(); j++)
		{
		cout << j << endl;
		fTreatmentFactory[i]->GetIVStock()[j].Print();
		}
	}

	for(int i = 0; i < (int)fReactor.size(); i++)
	{	
		cout << "Reactor" << endl;
		fReactor[i]->GetIVReactor().Print();
	}
	DBGL;
}
void CLASS::ResetQuantity()
{
	DBGL;
	IsotopicVector EmptyIV;
	fTotalInReactor		= EmptyIV;
	fTotalWaste		= EmptyIV;
	fTotalStock		= EmptyIV;
	fTotalGodIncome		= EmptyIV;
	fTotalCooling		= EmptyIV;
	fTotalSeparating	= EmptyIV;
	fIVInCycleTotal		= EmptyIV;
	fIVTotal		= EmptyIV;
	fIVTotal.cleartest();
	DBGL;
}

//________________________________________________________________________
void CLASS::OpenOutputTree()
{
	DBGL;

	cout << "Opening : " << fOutputName << " ...";

	fOutTree = new TFile(fOutputName.c_str(),"recreate");

	if(!fOutTree)
	{
		cout << "\nCould not open " << fOutputName <<endl;
		exit(-1);
	}
	cout << "\t ...O.K." << endl;

	fOutT = new TTree("Data","Data Tree");
	cout << "Creating Data Tree...";
	if(!fOutTree)
	{
		cout << "\nCould not create Data Tree in " << fOutputName << endl;
		exit(-1);
	}
	cout << "\t ...O.K!" << endl;
	// fOutT->SetAutoFlush(1000000);
	// outT->SetAutoSave(100000);
}
void CLASS::CloseOutputTree()
{
	DBGL;

	fOutTree->ls();
	cout << "Writing outTree " << fOutputName << endl;
	fOutTree->Write();

	if(fOutTree->IsOpen()) {
		cout << "Deleting outTree : " << endl;
		delete fOutT;
		cout << "Closing file : " << fOutputName <<endl;
		fOutTree-> Close();
		delete fOutTree;
	} else {
		cout << "File was not opened " << fOutputName << endl;
		exit(-1);
	}
}
//________________________________________________________________________
void CLASS::OutAttach()
{
	DBGL;
	ResetQuantity();
	//Branch Absolut Time
	fOutT->Branch("AbsTime",&fAbsoluteTime,"AbsoluteTime/L");
	
	// Branch the Sum IV

	
	fOutT->Branch("TF_WASTE.", "IsotopicVector", &fTotalWaste);
	fOutT->Branch("TF_STOCK.", "IsotopicVector", &fTotalStock);
	fOutT->Branch("TF_SEPARATING.", "IsotopicVector", &fTotalSeparating);
	fOutT->Branch("TF_COOLING.", "IsotopicVector", &fTotalCooling);
	fOutT->Branch("GODINCOME.", "IsotopicVector", &fTotalGodIncome);
	fOutT->Branch("INCYCLE.", "IsotopicVector", &fIVInCycleTotal);
	fOutT->Branch("TOTAL.", "IsotopicVector", &fIVTotal);
	
	fOutT->Branch("TreatmentFactory", "TreatmentFactory", &fTreatmentFactory.at(0));
	fOutT->Branch("Reactor", "Reactor", &fReactor.at(0));
	
	DBGL;
}



//________________________________________________________________________
void CLASS::Write()
{	
	DBGL;
	string basename = "CLASS";


//--- Reactor Writing ---/

	//*** Individual ***//
//	for(int i = 0; i < (int)fReactor.size(); i++)
//	{
//		ostringstream ostr;
//		ostr << fReactorIndex[i] ;
//		string Rbasename = basename + "_R" + ostr.str();
//		fReactor[i]->Write(Rbasename);
//	}

	//****** Total *****//
	IsotopicVector TotalInReactor;
	for(int i = 0; i < (int)fReactor.size(); i++)
	{	
		TotalInReactor += fReactor[i]->GetIVReactor();
	}
	string RTotal = basename + "_R_TOTAL";
	TotalInReactor.Write(RTotal, fAbsoluteTime);



//------ TF Writing -----/

	//*** Individual ***//
//	for(int i = 0; i < (int) fTreatmentFactory.size();i++)
//	{
//		ostringstream ostr;
//		ostr << fTreatmentFactoryIndex[i] ;
//		string TFbasename = basename + "_TF" + ostr.str();
//		fTreatmentFactory[i]->Write(TFbasename);
//	}

	//****** Total *****//
	IsotopicVector TotalWaste;
	IsotopicVector TotalGodIncome;
	IsotopicVector TotalStock;
	IsotopicVector TotalCooling;
	IsotopicVector TotalSeparating;


	for(int i = 0; i < (int) fTreatmentFactory.size();i++)
	{
		TotalWaste += fTreatmentFactory[i]->GetIVWaste();
		TotalStock += fTreatmentFactory[i]->GetIVFullStock();
		TotalGodIncome += fTreatmentFactory[i]->GetIVGodIncome();
		
		for(int j=0; j < (int)fTreatmentFactory[i]->GetIVCooling().size(); j++)
			TotalCooling += fTreatmentFactory[i]->GetIVCooling()[j];
		for(int j=0; j < (int)fTreatmentFactory[i]->GetIVSeparating().size(); j++)
			TotalSeparating += fTreatmentFactory[i]->GetIVSeparating()[j];
	}

	string TFWaste = basename + "_TF_WASTE";
	string TFStock = basename + "_TF_STOCK";
	string TFGodIncome  = basename + "_TF_GODINCOME";
	string TFSeparating = basename + "_TF_SPERATING";

	TotalWaste.Write(TFWaste, fAbsoluteTime);
	TotalStock.Write(TFStock, fAbsoluteTime);
	TotalGodIncome.Write(TFGodIncome, fAbsoluteTime);
	TotalSeparating.Write(TFSeparating, fAbsoluteTime);

//---- Total Writing ----/

	string InCycleTotal 	= basename + "_IN_CYCLE_TOTAL";
	string Total 		= basename + "_TOTAL";

	IsotopicVector IVInCycleTotal;
	IsotopicVector IVTotal;

	IVInCycleTotal 	= TotalCooling + TotalInReactor + TotalSeparating + TotalStock;
	IVTotal 	= IVInCycleTotal + TotalWaste;

	IVInCycleTotal.Write(InCycleTotal, fAbsoluteTime);
	IVTotal.Write(Total, fAbsoluteTime);


	DBGL;
}

