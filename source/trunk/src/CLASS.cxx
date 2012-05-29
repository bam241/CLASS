#include "CLASS.hxx"
#include "Reactor.hxx"
#include "TreatmentFactory.hxx"
#include "Defines.hxx"
#include "IsotopicVector.hxx"
#include "LogFile.hxx"

#include <ctime>
#include "time.h"
#include <cmath>
#include <iomanip>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <omp.h>
#include <TFile.h>
#include <TTree.h>
//#include <gsl/gsl_vector.h>
//#include <gsl/gsl_multimin.h>


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
	fBuildingMethod = 0;
	
	fOutputName = "CLASS_Default.root";
	string logname = "CLASS.log";
	fLog = new LogFile("CLASS.log");

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
	fBuildingMethod = 0;
	
	fOutputName = "CLASS_Default.root";
	string logname = "CLASS.log";
	fLog = new LogFile("CLASS.log");

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
	fBuildingMethod = 0;
	fOutputName = name;
	string logname = "CLASS.log";
	fLog = new LogFile("CLASS.log");

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
	fTreatmentFactory[fTreatmentFactory.size()-1]->SetParc(this);
	fTreatmentFactory[fTreatmentFactory.size()-1]->SetLog(fLog);
	fTreatmentFactory[fTreatmentFactory.size()-1]->SetStockManagement(fStockManagement);
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
	fReactor.back()->SetLog(fLog);
	fReactorLastIndex++;
	fReactorIndex.push_back(fReactorLastIndex);
DBGL;
}

//________________________________________________________________________
//void CLASS::RemoveReactor()
//{
//DBGL;
//	for(int i = (int)fReactor.size()-1; i >= 0; i--)
//	{
//		if(fAbsoluteTime == fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime())
//		{
//			fReactor[i]->GetAssociedTreatmentFactory()->AddIVCooling(fReactor[i]->GetIVReactor());
//			delete (fReactor[i]);
//			fReactor.erase(fReactor.begin()+i);
//			fReactorIndex.erase(fReactorIndex.begin()+i);
//		}
//	}
//DBGL;
//}

//________________________________________________________________________
//____________________________ Building Method ___________________________
//________________________________________________________________________

void CLASS::UpdateParcStock()
{
DBGL;
	
	fParcStock.clear();
	
	for(int i=0; i< (int)fTreatmentFactory.size(); i++ )
	{
		for(int j = 0; j < (int)fTreatmentFactory[i]->GetIVStock().size(); j++ )
		{
			fParcStock.insert( pair< pair< int ,int > , IsotopicVector > 
				( pair< int ,int >(i, j), fTreatmentFactory[i]->GetIVStock()[j]) ) ;
		}
	}
	
DBGL;
}

//________________________________________________________________________
IsotopicVector CLASS::BuildIsotopicVector(IsotopicVector isotopicvector)
{
DBGL;
	UpdateParcStock();
	
	IsotopicVector BuildIV;
	if (fStockManagement == false )
		return BuildIV;

	if(fBuildingMethod == 0)
		return MonteCarloBuild(isotopicvector);
	
//	else if(fBuildingMethod == 1)
//		return MinimizationBuild(isotopicvector);
	else 
	{
		cout		<< "!!Warning!! !!!CLASS!!! Wrong BUilding Method !\n" << endl;
		fLog->fLog 	<< "!!Warning!! !!!CLASS!!! Wrong BUilding Method !\n" << endl;
		exit(1);
	}
	
DBGL;
}


//________________________________________________________________________
void CLASS::DumpParcStock()
{
DBGL;
	if(fStockManagement == false) return;
	for(int i = 0; i < (int)fTreatmentFactory.size(); i++)
		fTreatmentFactory[i]->ClearIVStock();
	
	
	map< pair<int,int>, IsotopicVector >::iterator it;
	for(it = fParcStock.begin(); it != fParcStock.end(); it++)
		fTreatmentFactory[(*it).first.first]->AddIVStock((*it).second);
	
DBGL;	
}




//________________________________________________________________________

//________________________________________________________________________

/*struct MinimizationParameter {
	int z;
	vector<int> FactoryIndex;
	vector<int> StockIndex;
	map< pair<int,int>, IsotopicVector > ParcStock;
	IsotopicVector isotopicvector;
	IsotopicVector TotalStock;

};


double	my_f (const gsl_vector *v, void *params)
{
	
	MinimizationParameter *par = (MinimizationParameter *)params;
	int z = (*par).z;
	vector<int> FactoryIndex = (*par).FactoryIndex;
	vector<int> StockIndex	 = (*par).StockIndex;
	map< pair<int,int>, IsotopicVector > ParcStock = (*par).ParcStock;
	IsotopicVector isotopicvector = (*par).isotopicvector;
	IsotopicVector TotalStock = (*par).TotalStock;

	double x[FactoryIndex.size()];
	
		IsotopicVector IVSum;
	for (int i =0; i < (int)FactoryIndex.size(); i++)
	{
		x[i] = gsl_vector_get(v, i);
		if(x[i]<0) x[i] = -x[i];
	}

	
	for (int i = 0; i < (int)FactoryIndex.size(); i++)
	{
		map< pair<int,int>, IsotopicVector >::iterator it;
		it = ParcStock.find(pair<int,int>(FactoryIndex[i], StockIndex[i] ));
		IVSum += (*it).second.GetAtomicComposition(z)*x[i]/Norme((*it).second.GetAtomicComposition(z));
	}
	return Distance(IVSum,isotopicvector.GetAtomicComposition(z)); 
}

// The gradient of f, df = (df/dx, df/dy). 
void my_df (const gsl_vector *v, void *params, gsl_vector *df)
{
	MinimizationParameter *par = (MinimizationParameter *)params;
	int z = (*par).z;
	vector<int> FactoryIndex = (*par).FactoryIndex;
	vector<int> StockIndex	 = (*par).StockIndex;
	map< pair<int,int>, IsotopicVector > ParcStock = (*par).ParcStock;
	IsotopicVector isotopicvector = (*par).isotopicvector;
	IsotopicVector TotalStock = (*par).TotalStock;


	double x[FactoryIndex.size()];	
	
	for (int i =0; i < (int)FactoryIndex.size(); i++)
	{
		x[i] = gsl_vector_get(v, i);
		if(x[i]<0) x[i] = -x[i];
	}
	for (int i = 0; i < (int)FactoryIndex.size(); i++)
	{
		double Sum = 0;
		map< pair<int,int>, IsotopicVector >::iterator it_i;
		it_i = ParcStock.find(pair<int,int>(FactoryIndex[i], StockIndex[i] ));
		map<ZAI ,double > isotopicquantity = TotalStock.GetAtomicComposition(z).GetIsotopicQuantity();
		for(map<ZAI ,double >::iterator it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		{	
			Sum += 2*(*it_i).second.GetAtomicComposition(z).GetZAIIsotopicQuantity((*it).first)
			*isotopicvector.GetAtomicComposition(z).GetZAIIsotopicQuantity((*it).first);
			
			for (int j = 0; j < (int)FactoryIndex.size(); j++) 
			{
				//if(i == j) j++;

				map< pair<int,int>, IsotopicVector >::iterator it_j;
				it_j = ParcStock.find(pair<int,int>(FactoryIndex[j], StockIndex[j] ));
				
				Sum -= 2*x[j]*(*it_i).second.GetZAIIsotopicQuantity((*it).first)
				*(*it_j).second.GetZAIIsotopicQuantity((*it).first);
			}
		}
		gsl_vector_set(df, i, Sum);
		
	}
}

// Compute both f and df together. 
void my_fdf (const gsl_vector *x, void *params, double *f, gsl_vector *df) 
{
	*f = my_f(x, params); 
	my_df(x, params, df);
}



IsotopicVector CLASS::MinimizationBuild(IsotopicVector isotopicvector )
{
DBGL;
	IsotopicVector BuildIV;
	vector<int>		FactoryIndex;
	vector<int>		StockIndex;
	
	map< pair<int,int>, IsotopicVector >::iterator	it;

	for(it = fParcStock.begin(); it != fParcStock.end(); it++ )
	{
			FactoryIndex.push_back((*it).first.first);
			StockIndex.push_back((*it).first.second);
	}
	
	
	if( (int)FactoryIndex.size() <= 3) return BuildIV;
	
	for(int k = 0; k < (int)isotopicvector.GetAtomicSpecies().size(); k++ ) // Loop on the Atomic Species
	{
		int z = isotopicvector.GetAtomicSpecies()[k]; 		// Get the Atomic Species Number
		IsotopicVector BuildIVz;
		IsotopicVector IVAtomicCompo = isotopicvector.GetAtomicComposition(z);
		//double delta = 1e-3*Norme(IVAtomicCompo);

		
		
		
				
		size_t iter = 0;
		int status;
		
		const gsl_multimin_fdfminimizer_type *T;
		gsl_multimin_fdfminimizer *s;
		
		// Index of IV 
		MinimizationParameter par ;
		par.z = z;
		par.FactoryIndex = FactoryIndex;
		par.StockIndex	 = StockIndex;
		par.ParcStock	 = fParcStock;
		par.isotopicvector = isotopicvector;
		par.TotalStock = fTotalStock;
		//double par[5] = { 1.0, 2.0, 10.0, 20.0, 30.0 };

		gsl_vector *x;
		gsl_multimin_function_fdf my_func;

		my_func.n = (int)FactoryIndex.size();  // number of function components 
		my_func.f = my_f;
		my_func.df = my_df;
		my_func.fdf = my_fdf;
		my_func.params = &par;
		
		// Starting point, x = (5,7) 
		x = gsl_vector_alloc (FactoryIndex.size());
		for (int i = 0; i < (int)FactoryIndex.size(); i++) 
		{
			map< pair<int,int>, IsotopicVector >::iterator	it;
			it = fParcStock.find(pair<int,int>(FactoryIndex[i], StockIndex[i] ));
			gsl_vector_set (x, i, Norme((*it).second.GetAtomicComposition(z))/2. );
		}
		T = gsl_multimin_fdfminimizer_conjugate_fr;
		s = gsl_multimin_fdfminimizer_alloc (T, FactoryIndex.size());
		
		gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-1);
		
		do
		{
			iter++;
			status = gsl_multimin_fdfminimizer_iterate (s);
			
			if (status)
				break;
			
			status = gsl_multimin_test_gradient (s->gradient, 1e-3);
			//status = gsl_multimin_test_size (0, 1e-3);
			if (status == GSL_SUCCESS)
				printf ("Minimum found !\n");
			
		}
		while (status == GSL_CONTINUE && iter < 100);
		
		
		for (int i = 0; i < (int)FactoryIndex.size(); i++)
		{
			map< pair<int,int>, IsotopicVector >::iterator it;
			it = fParcStock.find(pair<int,int>(FactoryIndex[i], StockIndex[i] ));
			double Contribution_i = gsl_vector_get (s->x, i);
			BuildIVz += (*it).second.GetAtomicComposition(z)*Contribution_i/Norme((*it).second.GetAtomicComposition(z));
			(*it).second -= (*it).second.GetAtomicComposition(z)*Contribution_i/Norme((*it).second.GetAtomicComposition(z));
		}
		
		
		gsl_multimin_fdfminimizer_free (s);
		gsl_vector_free (x);
		
		BuildIV += BuildIVz;
	}
	return BuildIV;
DBGL;	
}
*/

//________________________________________________________________________
IsotopicVector CLASS::MonteCarloBuild(IsotopicVector isotopicvector )
{

	IsotopicVector BuildIV;
	
	
	int z = 0;
	
	vector<int>		FactoryIndex;
	vector<int>		StockIndex;
	map< pair<int,int>, IsotopicVector >::iterator	it;

	for(it = fParcStock.begin(); it != fParcStock.end(); it++ )
	{
			FactoryIndex.push_back((*it).first.first);
			StockIndex.push_back((*it).first.second);
	}
	
	if( (int)FactoryIndex.size() == 0) return BuildIV;

	for(int k = 0; k < (int)isotopicvector.GetAtomicSpecies().size(); k++ ) // Loop on the Atomic Species
	{
		z = isotopicvector.GetAtomicSpecies()[k]; 		// Get the Atomic Species Number
		
		IsotopicVector IVAtomicCompo = isotopicvector.GetAtomicComposition(z);
		IsotopicVector BuildIVz;
		double delta = 1e-3*Norme(IVAtomicCompo);
		
		
		
		vector<int>		FactoryIndextmp	= FactoryIndex;
		vector<int>		StockIndextmp 	= StockIndex;
		
		double DistanceN = Norme(IVAtomicCompo);	// Distance at the Nst iteration
		double DistanceN_1 = Norme(IVAtomicCompo);	// Distance at the N-1 st iteration
		
		while(FactoryIndextmp.size() !=0)
		{
			int IVIndex = FactoryIndextmp.size();
			
			// Find a random IV contribution
			while(IVIndex == (int)FactoryIndextmp.size()) IVIndex = (int)random(0,FactoryIndextmp.size());
			
			map< pair<int,int>, IsotopicVector >::iterator it;
			it = fParcStock.find(pair<int,int>(FactoryIndextmp[IVIndex], StockIndextmp[IVIndex] ));	
			IsotopicVector IVToAdd = (*it).second.GetAtomicComposition(z) ; 
			
			
			IsotopicVector BuildIVzTmp = BuildIVz + IVToAdd/Norme(IVToAdd) * delta; // Build the new IVBuild 
			
			DistanceN = Distance( BuildIVzTmp, isotopicvector.GetAtomicComposition(z) ); //Progress determination
			
			if(DistanceN > DistanceN_1 || delta > Norme(IVToAdd) ) // Not Good Peak remove it from Index
			{
				FactoryIndextmp.erase(FactoryIndextmp.begin()+IVIndex);
				StockIndextmp.erase(StockIndextmp.begin()+IVIndex);
			}
			else	// Good Peak, Add the contribution to IVBuild and remove from stock...
			{
				// Validation of Values....
				DistanceN_1 = DistanceN;
				BuildIVz = BuildIVzTmp;
				
				// remove contribution from stock
				map< pair<int,int>, IsotopicVector >::iterator it;
				it = fParcStock.find(pair<int,int>(FactoryIndextmp[IVIndex], StockIndextmp[IVIndex] ));
				
				(*it).second -= IVToAdd/Norme(IVToAdd) * delta;
			}
		}
		BuildIV += BuildIVz;
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

//********* Printing Step *********//
	{
		long int step = 0;
		if(step >= fAbsoluteTime	)
			fTimeStep.insert( pair<long int ,int>(step,1) );
		step += fPrintStep;
		do 
		{
			
			if(step >= fAbsoluteTime	)
				fTimeStep.insert( pair<long int ,int>(step,1) );
			step += fPrintStep;
		}
		while( step < t );
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

			if(fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime()+coolingstep >= fAbsoluteTime 
			&& fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime()+coolingstep <= t)			// Set End of Cooling after end of reactor
			{
				pair< map<long int, int>::iterator, bool > IResult;
				IResult = fTimeStep.insert( pair<long int ,int>(fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime()+coolingstep,8) );
				if( IResult.second == false ) IResult.first->second |= 8;
			}

			if(fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime()+coolingstep+seprationstep >= fAbsoluteTime 
			&& fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime()+coolingstep+seprationstep <= t)	// Set end of Separation after end of reactor
			{
				pair< map<long int, int>::iterator, bool > IResult;
				IResult = fTimeStep.insert( pair<long int ,int>(fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime()+coolingstep+seprationstep,16) );
				if( IResult.second == false ) IResult.first->second |= 16;
			}

		}
		
		// Set End of reactor cycle
		if(step >= fAbsoluteTime &&  step <= t && step < fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime())		
		{
			pair< map<long int, int>::iterator, bool > IResult = fTimeStep.insert( pair<long int ,int>(step,4) );
			if( IResult.second == false ) IResult.first->second |= 4;
		
		}
		

//********* TF Evolution Step *********//
		step += fReactor[i]->GetCycleTime();
		do 
		{
			
			if(step > fAbsoluteTime && step <= t && step < fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime())
			{						// Set End of reactor cycle
				pair< map<long int, int>::iterator, bool > IResult = fTimeStep.insert( pair<long int ,int>(step,4) );
				if( IResult.second == false ) IResult.first->second  |= 4;
			}

			if(step >= fAbsoluteTime && step + coolingstep <= t)			// Set End of Cooling
			{
				pair< map<long int, int>::iterator, bool > IResult;
				IResult = fTimeStep.insert( pair<long int ,int>(step+coolingstep,8) );
				if( IResult.second == false ) IResult.first->second |= 8;
			}

			if(step >= fAbsoluteTime && step + coolingstep + seprationstep <= t)	// Set end of Separation
			{
				pair< map<long int, int>::iterator, bool > IResult;
				IResult = fTimeStep.insert( pair<long int ,int>(step+coolingstep+seprationstep,16) );
				if( IResult.second == false ) IResult.first->second |= 16;
			}
			step += fReactor[i]->GetCycleTime();
		}
		while(step <= t && step <= fReactor[i]->GetCreationTime() + fReactor[i]->GetLifeTime() );
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
	{
		cout		<< "!!Warning!! !!!CLASS!!! Can't open \" CLASS_TimeStep \"\n" << endl;
		fLog->fLog 	<< "!!Warning!! !!!CLASS!!! Can't open \" CLASS_TimeStep \"\n" << endl;
	}
	for(map<long int ,int >::iterator it = fTimeStep.begin(); it != fTimeStep.end(); it++)
		TimeStepfile << (*it).first << " " << (*it).second << endl;
	
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
	{
	TreatmentEvolution();

	for(int i = 0; i < (int) fTreatmentFactory.size();i++)
		fTreatmentFactory[i]->Dump();
	}

	#pragma omp section
	{
		#pragma omp parallel for
			for(int i = 0; i < (int)fReactor.size(); i++)
				fReactor[i]->Evolution(fAbsoluteTime);
	}
	}
	


	UpdateParcStock();
	for(int i = 0; i < (int)fReactor.size(); i++)
		fReactor[i]->Dump();

	DumpParcStock();
	
DBGL;
}

//________________________________________________________________________
void CLASS::Evolution(long int t)
{
DBGL;
	OpenOutputTree();
	OutAttach();
	BuildTimeVector(t);
	fOutT->Fill();
	int Start = time(NULL);
	for(map<long int ,int >::iterator it = fTimeStep.begin(); it != fTimeStep.end(); it++)
	{
		ResetQuantity();
		fAbsoluteTime = (*it).first;
		if( (*it).second & 2 || (*it).second & 1 )
		{
			ReactorEvolution();
			if((*it).second & 2 )
//				RemoveReactor();
		
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
		{
#pragma omp single
				{
				UpdateParc();
				fOutT->Fill();
				ProgressPrintout(Start, t);
				}
		}
	}
	cout << endl;
#pragma omp single
	{CloseOutputTree();}
DBGL;
}

void CLASS::ProgressPrintout(int starttime, long int t)
{
DBGL;
	int TimeNow = time(NULL);
	int Spent = (int)difftime(TimeNow, starttime);
	double Time = fAbsoluteTime/3600/24/365.4;
	double Total = t/3600/24/365.4;
	double Remain =  (Spent/Time * Total - Spent);
	cout << "                                                                                             " << flush ;
	cout << "\rProcessed " << setprecision(4) << Time << " / " << setprecision(4) << Total << " STEP "
	     << setprecision(4) << (Time/Total*100.0) << "%) ---  I still need : "
	<< (long int)Remain/60 << " min "<<  (long int)(Remain)%60 << " sec to finish ! \r" << flush;
DBGL;
}

//________________________________________________________________________
//______________________________ Out Method ______________________________
//________________________________________________________________________
void CLASS::UpdateParc()
{
DBGL;

	for(int i = 0; i < (int) fTreatmentFactory.size();i++)
	{
		fTotalWaste += fTreatmentFactory[i]->GetIVWaste();
		fTotalStock += fTreatmentFactory[i]->GetIVFullStock();
		fTotalGodIncome += fTreatmentFactory[i]->GetIVGodIncome();
		
		for(int j=0; j < (int)fTreatmentFactory[i]->GetIVCooling().size(); j++)
			fTotalCooling += fTreatmentFactory[i]->GetIVCooling()[j];
		for(int j=0; j < (int)fTreatmentFactory[i]->GetIVSeparating().size(); j++)
			fTotalSeparating += fTreatmentFactory[i]->GetIVSeparating()[j];
	}
	
	for(int i = 0; i < (int)fReactor.size(); i++)
		fTotalInReactor += fReactor[i]->GetIVReactor();

	fIVTotal = fTotalWaste + fTotalStock + fTotalCooling + fTotalSeparating + fTotalInReactor;
	fIVInCycleTotal = fTotalStock + fTotalCooling + fTotalSeparating + fTotalInReactor;
	
DBGL;	
}



void CLASS::ResetQuantity()
{
DBGL;
	fTotalInReactor.Clear();
	fTotalWaste.Clear();
	fTotalStock.Clear();
	fTotalGodIncome.Clear();
	fTotalCooling.Clear();
	fTotalSeparating.Clear();
	fIVInCycleTotal.Clear();
	fIVTotal.Clear();
DBGL;
}

//________________________________________________________________________
void CLASS::OpenOutputTree()
{
DBGL;

	cout << "Opening : " << fOutputName << " ...";
	fLog->fLog << "Opening : " << fOutputName << " ...";
	fOutTree = new TFile(fOutputName.c_str(),"UPDATE");

	if(!fOutTree)
	{
		cout << "\nCould not open " << fOutputName <<endl;
		fLog->fLog << "\nCould not open " << fOutputName <<endl;
		exit(-1);
	}
	cout << "\t ...O.K." << endl;

	fOutT = new TTree("Data","Data Tree");
	cout << "Creating Data Tree...";
	fLog->fLog << "Creating Data Tree...";
	if(!fOutTree)
	{
		cout << "\nCould not create Data Tree in " << fOutputName << endl;
		fLog->fLog << "\nCould not create Data Tree in " << fOutputName << endl;
		exit(-1);
	}
	cout << "\t ...O.K!" << endl;
	fLog->fLog <<  "\t ...O.K!" << endl;
	// fOutT->SetAutoFlush(1000000);
	// fOutT->SetAutoSave(100000);
}
void CLASS::CloseOutputTree()
{
DBGL;

	fOutTree->ls();
	cout << "Writing outTree " << fOutputName << endl;
	fLog->fLog << "Writing outTree " << fOutputName << endl;
	fOutTree->Write();

	if(fOutTree->IsOpen()) {
		cout << "Deleting outTree : " << endl;
		fLog->fLog << "Deleting outTree : " << endl;
		delete fOutT;
		cout << "Closing file : " << fOutputName <<endl;
		fLog->fLog << "Closing file : " << fOutputName <<endl;
		fOutTree-> Close();
		delete fOutTree;
	} else {
		cout << "File was not opened " << fOutputName << endl;
		fLog->fLog << "File was not opened " << fOutputName << endl;
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
	fOutT->Branch("INTOTALINREACTOR.", "IsotopicVector", &fTotalInReactor);
	fOutT->Branch("INCYCLE.", "IsotopicVector", &fIVInCycleTotal);
	fOutT->Branch("TOTAL.", "IsotopicVector", &fIVTotal);
	
	fOutT->Branch("TreatmentFactory", "TreatmentFactory", &fTreatmentFactory.at(0));
//	fOutT->Branch("Reactor", "Reactor", &fReactor.at(0));
	
DBGL;
}

//________________________________________________________________________
void CLASS::Write()
{	
DBGL;
	string basename = "CLASS";

	//****** Total *****//
	IsotopicVector TotalInReactor;
	for(int i = 0; i < (int)fReactor.size(); i++)
	{	
		TotalInReactor += fReactor[i]->GetIVReactor();
	}
	string RTotal = basename + "_R_TOTAL";
	TotalInReactor.Write(RTotal, fAbsoluteTime);

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
