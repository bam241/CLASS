#include "CLASSHeaders.hxx"
#include "CLASS.hxx"
#include <ctime>
#include <cmath>
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
		if(fAbsoluteTime == fReactor.at(i)->GetCreationTime() + fReactor.at(i)->GetLifeTime())
		{
			fReactor.at(i)->GetAssociedTreatmentFactory()->AddIVCooling(fReactor.at(i)->GetIVReactor());
			delete (fReactor.at(i));
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
	vector< map<pair<int,int>, double> > Zcontribution;
	

	for(int k = 0; k < (int)isotopicvector.GetAtomicSpecies().size(); k++ ) // Loop on the Atomic Species
	{
		int z = isotopicvector.GetAtomicSpecies().at(k); 		// Get the Atomic Species Number

		IsotopicVector IVAtomicCompo = isotopicvector.GetAtomicComposition(z);
		double delta = 1e-5*Norme(IVAtomicCompo);
		vector<IsotopicVector>	IVStock;
		vector<double>		IVMaxQuantity;
		vector<int>		FactoryIndex;
		vector<int>		StockIndex;
		
		for(int i=0; i< (int)fTreatmentFactory.size(); i++ )
		{
			for(int j = 0; j < (int)fTreatmentFactory.at(i)->GetIVStock().size(); j++ )
			{
				FactoryIndex.push_back(i);
				StockIndex.push_back(j);
				double Quantitytmp = Norme(fTreatmentFactory.at(i)->GetIVStock().at(j).GetAtomicComposition(z));
				IVMaxQuantity.push_back(Quantitytmp);
				IVStock.push_back( fTreatmentFactory.at(i)->GetIVStock().at(j).GetAtomicComposition(z) / Quantitytmp );
			}
		}
		
		
		vector<IsotopicVector>	IVStocktmp 	= IVStock;
		vector<int>		FactoryIndextmp	= FactoryIndex;
		vector<int>		StockIndextmp 	= StockIndex;

		map<pair<int,int>, double> contribution;
		
		
		
		
		// Choose the Contribution
		IsotopicVector BuildIVz;
		IsotopicVector BuildIVztmp;
		double DistanceN = Norme(IVAtomicCompo);	// Distance at the Nst iteration
		double DistanceN_1 = Norme(IVAtomicCompo);	// Distance at the N-1 st iteration

		while(IVStocktmp.size() !=0)
		{
			int IVIndex = IVStocktmp.size();
			double IVIndexQuantityneed = 0;
			
			
			while(IVIndex == (int)IVStocktmp.size()) IVIndex = (int)random(0,IVStocktmp.size());
		
			BuildIVztmp = BuildIVz + (IVStocktmp.at(IVIndex) * delta);
			DistanceN = Distance( BuildIVztmp, isotopicvector.GetAtomicComposition(z) );
			
			// Get the contribution of the IVStock IVIndex
			map< pair<int,int>, double >::iterator it = 
					contribution.find(pair<int,int>(FactoryIndextmp.at(IVIndex), StockIndextmp.at(IVIndex) ));
			if(it != contribution.end() )
				IVIndexQuantityneed = delta + (*it).second;
			
			

		
			if(DistanceN > DistanceN_1 || IVIndexQuantityneed > IVMaxQuantity.at(IVIndex) )
			{
				IVStocktmp.erase(IVStocktmp.begin()+IVIndex);
				FactoryIndextmp.erase(FactoryIndextmp.begin()+IVIndex);
				StockIndextmp.erase(StockIndextmp.begin()+IVIndex);
			}
			else
			{
				// Add the delta contribution from the IVIndex Stock IV²
				pair< map< pair<int,int>, double >::iterator, bool > IResult;
				
				IResult = contribution.insert( pair< pair< int ,int > , double > 
					( pair< int ,int >(FactoryIndextmp.at(IVIndex), StockIndextmp.at(IVIndex) ), delta) ) ;

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
		
		BuildIV += BuildIVz;
		Zcontribution.push_back( contribution );
	}

	for(int k = 0; k < (int)isotopicvector.GetAtomicSpecies().size(); k++ ) //Loop on the Atomic Species
	{
		// Take it and build the IsotopicVector
		map<pair<int,int>, double> contribution = Zcontribution.at(k);
		for(map<pair<int,int>, double>::iterator it = contribution.begin(); it != contribution.end(); it++)
		{
			int z = isotopicvector.GetAtomicSpecies().at(k);
			double IVNorme = Norme(fTreatmentFactory.at((*it).first.first)->GetIVStock().at((*it).first.second).GetAtomicComposition(z));
			IsotopicVector IV = fTreatmentFactory.at((*it).first.first)->GetIVStock().at((*it).first.second).GetAtomicComposition(z);
			fTreatmentFactory.at((*it).first.first)->TakeFromStock( IV / IVNorme * (*it).second , (*it).first.second );
		}
		
		
	 }

	
	DBGL;
	return BuildIV;
}

//________________________________________________________________________
void CLASS::BuildTimeVector(long int t)
{
	DBGL;
	fTimeStep.clear();
	fTimeStep.insert( pair<long int ,int>(t,1) );
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
	
		long int step = fReactor.at(i)->GetCreationTime();
		long int coolingstep = fReactor.at(i)->GetAssociedTreatmentFactory()->GetCoolingTime();
		long int seprationstep = fReactor.at(i)->GetAssociedTreatmentFactory()->GetSeparationTime();
		if(fReactor.at(i)->GetCreationTime() + fReactor.at(i)->GetLifeTime() > fAbsoluteTime && fReactor.at(i)->GetCreationTime() + fReactor.at(i)->GetLifeTime()<t)
		{						// set destruction of a reactor
			pair< map<long int, int>::iterator, bool > IResult ;
			IResult = fTimeStep.insert( pair<long int ,int>(fReactor.at(i)->GetCreationTime() + fReactor.at(i)->GetLifeTime(),2) );
			if( IResult.second == false ) IResult.first->second |= 2;
		}

		if(step >= fAbsoluteTime &&  step <= t && step < fReactor.at(i)->GetCreationTime() + fReactor.at(i)->GetLifeTime())		// Set End of reactor cycle
		{
			pair< map<long int, int>::iterator, bool > IResult = fTimeStep.insert( pair<long int ,int>(step,4) );
			if( IResult.second == false ) IResult.first->second |= 4;
		}
		


		while(step <= t && step <= fReactor.at(i)->GetCreationTime() + fReactor.at(i)->GetLifeTime() )
		{
			step += fReactor.at(i)->GetCycleTime();
			if(step > fAbsoluteTime && step <= t && step < fReactor.at(i)->GetCreationTime() + fReactor.at(i)->GetLifeTime())
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
	
	
	for(int i =0; i < (int)fTreatmentFactory.size(); i++) 
	{
		
		for(int j = 0; j<(int)fTreatmentFactory.at(i)->GetIVCooling().size(); j++ )// Set End of Cooling
		{
			if(fTreatmentFactory.at(i)->GetCoolingStartingTime().at(j) +  fTreatmentFactory.at(i)->GetCoolingTime() > fAbsoluteTime )
			{
				pair< map<long int, int>::iterator, bool > IResult;
				IResult = fTimeStep.insert( pair<long int ,int>(fTreatmentFactory.at(i)->GetCoolingStartingTime().at(j) +  fTreatmentFactory.at(i)->GetCoolingTime(),8) );
				if( IResult.second == false ) IResult.first->second |= 8;
			}
		}
		for(int j = 0; j< (int)fTreatmentFactory.at(i)->GetIVSeparating().size(); j++ )// Set end of Separation
		{
			if(fTreatmentFactory.at(i)->GetSeparatingStartingTime().at(j) +  fTreatmentFactory.at(i)->GetSeparationTime() > fAbsoluteTime )
			{
				pair< map<long int, int>::iterator, bool > IResult = fTimeStep.insert( pair<long int ,int>(fTreatmentFactory.at(i)->GetSeparatingStartingTime().at(j) +  fTreatmentFactory.at(i)->GetSeparationTime(),16) );
				if( IResult.second == false ) IResult.first->second |= 16;
			}
		}
	}

	
	
// //Print the Time Index
	ofstream TimeStepfile("CLASS_TimeStep", ios_base::app);		// Open the File
	if(!TimeStepfile)
	cout << "!!Warning!! !!!CLASS!!! Can't open \" CLASS_TimeStep \"\n" << endl;

	for(map<long int ,int >::iterator it = fTimeStep.begin(); it != fTimeStep.end(); it++)
	{
		TimeStepfile << (*it).first/365.4/3600/24 << " " << (*it).second << endl;
	}
	
	DBGL;	
}



//________________________________________________________________________
//___________________________ Evolution Method ___________________________
//________________________________________________________________________

void CLASS::TreatmentEvolution()
{
	DBGL;
	for(int i = 0; i < (int) fTreatmentFactory.size();i++)
		fTreatmentFactory.at(i)->Evolution(fAbsoluteTime);
	DBGL;
}

//________________________________________________________________________
void CLASS::ReactorEvolution()
{
	DBGL;
	for(int i = 0; i < (int)fReactor.size(); i++)
		fReactor.at(i)->Evolution(fAbsoluteTime);
	DBGL;
}

//________________________________________________________________________
void CLASS::Evolution(long int t)
{
	DBGL;
	BuildTimeVector(t);

	for(map<long int ,int >::iterator it = fTimeStep.begin(); it != fTimeStep.end(); it++)
	{
		fAbsoluteTime = (*it).first;
		if(fAbsoluteTime > t) return;
		
		if( (*it).second & 2 || (*it).second & 1 )
		{
			TreatmentEvolution();
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
			TreatmentEvolution();
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
			if((*it).second & 8 )
				(*it).second ^= 8;
			if((*it).second & 16 )
				(*it).second ^= 16;
		}
		
		if( (*it).second & 1 )
		{
			
//			Print();
			Write();
		}
			
	}
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
		cout << fTreatmentFactory.at(i)->GetIVCooling().size()<< endl;
		cout << "Waste " << endl;
		fTreatmentFactory.at(i)->GetIVWaste().Print();
		cout << "Stock :" << endl;
		for(int j=0; j < (int)fTreatmentFactory.at(i)->GetIVStock().size(); j++)
			{
			cout << j << endl;
			fTreatmentFactory.at(i)->GetIVStock().at(j).Print();
			}
	}

	for(int i = 0; i < (int)fReactor.size(); i++)
	{	
		cout << "Reactor" << endl;
		fReactor.at(i)->GetIVReactor().Print();
	}
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
//		ostr << fReactorIndex.at(i) ;
//		string Rbasename = basename + "_R" + ostr.str();
//		fReactor.at(i)->Write(Rbasename);
//	}

	//****** Total *****//
	IsotopicVector TotalInReactor;
	for(int i = 0; i < (int)fReactor.size(); i++)
	{	
		TotalInReactor += fReactor.at(i)->GetIVReactor();
	}
	string RTotal = basename + "_R_TOTAL";
	TotalInReactor.Write(RTotal, fAbsoluteTime);



//------ TF Writing -----/

	//*** Individual ***//
//	for(int i = 0; i < (int) fTreatmentFactory.size();i++)
//	{
//		ostringstream ostr;
//		ostr << fTreatmentFactoryIndex.at(i) ;
//		string TFbasename = basename + "_TF" + ostr.str();
//		fTreatmentFactory.at(i)->Write(TFbasename);
//	}

	//****** Total *****//
	IsotopicVector TotalWaste;
	IsotopicVector TotalGodIncome;
	IsotopicVector TotalStock;
	IsotopicVector TotalCooling;
	IsotopicVector TotalSeparating;


	for(int i = 0; i < (int) fTreatmentFactory.size();i++)
	{
		TotalWaste += fTreatmentFactory.at(i)->GetIVWaste();
		TotalStock += fTreatmentFactory.at(i)->GetIVFullStock();
		TotalGodIncome += fTreatmentFactory.at(i)->GetIVGodIncome();
		
		for(int j=0; j < (int)fTreatmentFactory.at(i)->GetIVCooling().size(); j++)
			TotalCooling += fTreatmentFactory.at(i)->GetIVCooling().at(j);
		for(int j=0; j < (int)fTreatmentFactory.at(i)->GetIVSeparating().size(); j++)
			TotalSeparating += fTreatmentFactory.at(i)->GetIVSeparating().at(j);
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

	cout << "CLASS : " << (long int)fAbsoluteTime/3600/24/365.4 << " STEP DONE !" << endl;
	DBGL;
}

