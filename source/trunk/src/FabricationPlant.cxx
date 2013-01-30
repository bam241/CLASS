#include "FabricationPlant.hxx"

#include "Storage.hxx"
#include "FabricationPlant.hxx"
#include "Reactor.hxx"
#include "EvolutiveProduct.hxx"
#include "EvolutionDataBase.hxx"
#include "IsotopicVector.hxx"
#include "CLASS.hxx"
#include "Defines.hxx"
#include "LogFile.hxx"




//#include "Minuit2/FunctionMinimum.h"
//#include "Minuit2/MnMigrad.h"
//#include "Minuit2/MnMinimize.h"
//#include "Minuit2/MnMinos.h"
//#include "Minuit2/MnSimplex.h"
//#include "Minuit2/MnUserParameters.h"
//#include "Minuit2/FCNBase.h"
//#include "Minuit2/MnPrint.h"
#include "TMatrixT.h"

#include <sstream>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

//________________________________________________________________________
//
//		FabricationPlant
//
//
//
//
//________________________________________________________________________
ClassImp(FabricationPlant)

template <class T>  T random(T a, T b) //peak random numebr between a and b
{
	double range = pow(2., 31);
	srand(time(NULL)); //initialize the srand
	return (T)a + (T)(b-a)*rand()/range;
}


FabricationPlant::FabricationPlant()
{
DBGL;
	fChronologicalTimePriority = false;
DBGL;
}

FabricationPlant::FabricationPlant(Storage* storage, Storage* reusable, double fabircationtime)
{
DBGL;
	fFabricationTime = fabircationtime;
	fChronologicalTimePriority = false;
	fStorage = storage;
	fReUsable = reusable;
DBGL;
}



//________________________________________________________________________
FabricationPlant::~FabricationPlant()
{
DBGL;
DBGL;
}


//________________________________________________________________________
void	FabricationPlant::AddValorisableIV(ZAI zai, double factor)
{
DBGL;
	pair<map<ZAI, double>::iterator, bool> IResult;
	if(factor > 1) factor = 1;
	
	if(factor > 0)
	{
		IResult = fValorisableIV.insert( pair<ZAI ,double>(zai, factor));
		if(IResult.second == false)
			IResult.first->second = factor;
	}
DBGL;
}



//________________________________________________________________________
//_______________________________ Evolution ______________________________
//________________________________________________________________________
void FabricationPlant::Evolution(double t)
{
	DBGL;
		// Check if the FabricationPlant has been created ...
	if(t == fInternalTime && t!=0) return;
		// Make the evolution for the FabricationPlant ...
	FabricationPlantEvolution(t);
		// ... And Finaly update the AbsoluteInternalTime
	fInternalTime = t;
	DBGL;
}

//________________________________________________________________________
void FabricationPlant::FabricationPlantEvolution(double t)
{
	DBGL;
	
	
	map<int ,double >::iterator it;
	for( it = fReactorNextStep.begin(); it!= fReactorNextStep.end(); it++ )
	{
		if( t + fFabricationTime >= fParc->GetReactor()[ (*it).first ]->GetCreationTime()
		   && t + fFabricationTime < fParc->GetReactor()[ (*it).first ]->GetCreationTime() + fParc->GetReactor()[ (*it).first ]->GetLifeTime())
		{
			if( abs((*it).second - t) < 3600 )
			{
				BuildFuelForReactor( (*it).first );
				(*it).second += fParc->GetReactor()[ (*it).first ]->GetCycleTime();
			}
			else if ( (*it).second - fParc->GetReactor()[ (*it).first ]->GetCycleTime() + fFabricationTime > t )
			{
				map<int ,IsotopicVector >::iterator it2 = fReactorFuturIV.find( (*it).first );
				if (it2 != fReactorFuturIV.end())
					(*it2).second = GetDecay((*it2).second, fFabricationTime - ((*it).second - t) );
			}
		}
	}
	
	
	DBGL;
}


	//________________________________________________________________________
void FabricationPlant::BuildFuelForReactor(int ReactorId)
{
	DBGL;
	
	EvolutionDataBase<IsotopicVector>* FuelType = fParc->GetReactor()[ReactorId]->GetFuelType();
	
	
	if(FuelType->GetFuelType() != "MOX")
	{
		cout << "!!Bad Trouble!! !!!FabricationPlant!!! Try to do MOX with a not MOXed DB "<< endl;
		fLog->fLog << "!!Bad Trouble!! !!!FabricationPlant!!! Try to do MOX with a not MOXed DB" << endl;
		exit (1);
	}
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
	
	double HMmass = fParc->GetReactor()[ReactorId]->GetHeavyMetalMass();
	double BU = fParc->GetReactor()[ReactorId]->GetBurnUp();
	IsotopicVector FullUsedStock;
	IsotopicVector stock;
	
	bool FuelBuild = false;
	while(FuelBuild == false)
	{
		double nPu_0 = 0;
		double MPu_0 = 0;
		{
			map<ZAI ,double>::iterator it;
			
			map<ZAI ,double> isotopicquantity = GetDecay( FullUsedStock , fFabricationTime).GetSpeciesComposition(94).GetIsotopicQuantity();
			for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
				nPu_0 += (*it).second;
			
			isotopicquantity = FullUsedStock.GetSpeciesComposition(94).GetIsotopicQuantity();
			for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
				MPu_0 += (*it).second*ZAImass.find( (*it).first )->second/Na*1e-6;
		}
		stock = GetStockToRecycle();
		if( stock.GetZAIIsotopicQuantity(ZAI(-1,-1,-1)) == 1 )
		{
			{
				EvolutiveProduct evolutiondb;
				pair<map<int, EvolutiveProduct>::iterator, bool> IResult;
				IResult = fReactorFuturDB.insert( pair<int, EvolutiveProduct>(ReactorId,evolutiondb) );
				if(IResult.second == false)
					IResult.first->second = evolutiondb;
			}
			{
				IsotopicVector EmptyIV;
				pair<map<int, IsotopicVector>::iterator, bool> IResult;
				IResult = fReactorFuturIV.insert( pair<int, IsotopicVector>(ReactorId,EmptyIV) );
				if(IResult.second == false)
					IResult.first->second = EmptyIV;
			}
			FuelBuild = true;
			fFractionToTake.clear();
		}
		else
		{
			double nPu_1 = 0;
			double MPu_1 = 0;
			double Sum_AlphaI_nPuI = 0;
			double Sum_AlphaI_nPuI0 = 0;
			{
				map<ZAI ,double>::iterator it;
				map<ZAI ,double> isotopicquantity = GetDecay( stock , fFabricationTime).GetSpeciesComposition(94).GetIsotopicQuantity();

				for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
				{
					if ((*it).first.A() >= 238 && (*it).first.A() <= 242)
					{ 
						nPu_1 += (*it).second;
						Sum_AlphaI_nPuI += FuelType->GetPFuelParameter()[(*it).first.A() -237]*(*it).second;
					}
				}
				
				isotopicquantity = stock.GetSpeciesComposition(94).GetIsotopicQuantity();
				for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
				if ((*it).first.A() >= 238 && (*it).first.A() <= 242)
				{ 
					MPu_1 += (*it).second * (ZAImass.find( (*it).first )->second)/Na*1e-6;
				}
				
				isotopicquantity = GetDecay( FullUsedStock , fFabricationTime).GetSpeciesComposition(94).GetIsotopicQuantity();
				for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
					if ((*it).first.A() >= 238 && (*it).first.A() <= 242)
					{ 
						Sum_AlphaI_nPuI0 += FuelType->GetPFuelParameter()[(*it).first.A() -237]*(*it).second;
					}				
			}
			
			double StockFactionToUse = 0;
			   
			double NT = HMmass*1e6 * Na / ZAImass.find( ZAI(92,238,0) )->second;
			
			double N1 = (BU - FuelType->GetPFuelParameter()[6]) * NT;
			double N2 = -Sum_AlphaI_nPuI0;
			double N3 = -FuelType->GetPFuelParameter()[0] * Na / ZAImass.find( ZAI(92,238,0) )->second * (HMmass*1e6 - MPu_0*1e6);
			
			double D1 = Sum_AlphaI_nPuI;
			double D2 = -FuelType->GetPFuelParameter()[0] * MPu_1*1e6 * Na / ZAImass.find( ZAI(92,238,0) )->second;
			
			StockFactionToUse = (N1 + N2 + N3) / (D1 + D2);
			
			if(StockFactionToUse < 0)
			{
				cout << "!!Bad Trouble!! !!!FabricationPlant!!! Oups Bug in calculating stock fraction to use "<< endl;
				fLog->fLog << "!!Bad Trouble!! !!!FabricationPlant!!! Oups Bug in calculating stock fraction to use" << endl;
				exit (1);
			}
			
			if( StockFactionToUse > 1 )
			{
				
				FullUsedStock += stock;
				RecycleStock(1);
				FuelBuild = false;
			}
			else
			{
				RecycleStock(StockFactionToUse);
				
				IsotopicVector IVBeginCycle;
				FuelBuild = true;
				
				ZAI U8 = ZAI(92,238,0);
				double U8_Quantity =  (HMmass - (MPu_0+StockFactionToUse*MPu_1 ))/ZAImass.find( ZAI(92,238,0))->second*Na/1e-6;
				fParc->AddGodIncome( U8, U8_Quantity );
				
				for(int i = (int)fFractionToTake.size()-1; i >= 0; i--)
				{
					IVBeginCycle += fStorage->GetStock()[fFractionToTake[i].first].GetSpeciesComposition(94)*( fFractionToTake[i].second );
					pair<IsotopicVector,IsotopicVector> reste = Separation( fStorage->GetStock()[fFractionToTake[i].first]*(fFractionToTake[i].second)
											       - fStorage->GetStock()[fFractionToTake[i].first].GetSpeciesComposition(94)*(fFractionToTake[i].second) );
					
					fStorage->TakeFractionFromStock(fFractionToTake[i].first,fFractionToTake[i].second);
					fParc->AddWaste(reste.second);
					fReUsable->AddToStock(reste.first);
				}
				fFractionToTake.clear();
				
				IVBeginCycle += U8_Quantity*U8;
				EvolutiveProduct evolutiondb = BuildEvolutiveDB(ReactorId, IVBeginCycle);
				
				{
					pair<map<int, EvolutiveProduct>::iterator, bool> IResult;
					IResult = fReactorFuturDB.insert( pair<int, EvolutiveProduct>(ReactorId,evolutiondb) );
					if(IResult.second == false)
						IResult.first->second = evolutiondb;
				}
				{
					pair<map<int, IsotopicVector>::iterator, bool> IResult;
					IResult = fReactorFuturIV.insert( pair<int, IsotopicVector>(ReactorId,IVBeginCycle) );
					if(IResult.second == false)
						IResult.first->second = IVBeginCycle;
				}
			}
		}
	}
	
	DBGL;
}



//________________________________________________________________________
//_________________________________ Decay ________________________________
//________________________________________________________________________
IsotopicVector FabricationPlant::GetDecay(IsotopicVector isotopicvector, double t)
{
	DBGL;
	IsotopicVector IV;
	
	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	
	map<ZAI ,double >::iterator it;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
	{
		if((*it).second > 0)
		{
 			IsotopicVector ivtmp = fDecayDataBase->Evolution(it->first, t) * (*it).second ;
			IV += ivtmp;
		}
	}
	DBGL;
	return IV;
}


//________________________________________________________________________
//_____________________________ Reactor & DB _____________________________
//________________________________________________________________________
EvolutiveProduct FabricationPlant::BuildEvolutiveDB(int ReactorId,IsotopicVector isotopicvector)
{
DBGL;
	EvolutionDataBase<IsotopicVector>* evolutiondb = fParc->GetReactor()[ReactorId]->GetFuelType();

	isotopicvector = GetDecay(isotopicvector, fFabricationTime);
	
	map<double, EvolutiveProduct> distances = evolutiondb->GetDistancesTo(isotopicvector);
	

	EvolutiveProduct EvolBuild = distances.begin()->second.GenerateDBFor(isotopicvector);
	
	

	return EvolBuild;
DBGL;
}

//________________________________________________________________________
void FabricationPlant::TakeReactorFuel(int Id)
{
	DBGL;
	
	IsotopicVector IV;
	map<int ,IsotopicVector >::iterator it2 = fReactorFuturIV.find( Id );
	if (it2 != fReactorFuturIV.end())
		(*it2).second = IV;
	DBGL;
}

//________________________________________________________________________
EvolutiveProduct FabricationPlant::GetReactorEvolutionDB(int ReactorId)
{
	DBGL;
	map< int,EvolutiveProduct >::iterator it = fReactorFuturDB.find(ReactorId);
	return (*it).second;
	DBGL;
}

//________________________________________________________________________
//_______________________________ Storage ________________________________
//________________________________________________________________________

IsotopicVector FabricationPlant::GetStockToRecycle()
{
	DBGL;
	IsotopicVector NextStock;
	int IdToTake = -1;
	
	if(fChronologicalTimePriority == true)
		IdToTake = (int)( fFractionToTake.size() );
	else
		IdToTake = (int)( fStorage->GetStock().size() -1 - fFractionToTake.size() );
	if(0 <= IdToTake && IdToTake <= (int)fStorage->GetStock().size()-1)
	{
		NextStock = fStorage->GetStock()[IdToTake];
		fFractionToTake.push_back( pair<int,double>(IdToTake,0.) );
	}
	else NextStock += ZAI(-1,-1,-1) *1;
	
	return NextStock;
	DBGL;
}

//________________________________________________________________________
void FabricationPlant::RecycleStock(double fraction)
{
DBGL;
	fFractionToTake.back().second = fraction;
DBGL;
}



//________________________________________________________________________
void FabricationPlant::DumpStock()
{
DBGL;
	
	
DBGL;
}

//________________________________________________________________________
pair<IsotopicVector, IsotopicVector> FabricationPlant::Separation(IsotopicVector isotopicvector)
{
DBGL;
	//[0] = re-use ; [1] = waste
	pair<IsotopicVector, IsotopicVector>	IVTmp;
	
	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for(it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	{
		map<ZAI ,double>::iterator it2;
		it2 = fValorisableIV.find((*it).first);
		if ( it2 != fValorisableIV.end() )
		{
			IVTmp.first.Add(	(*it).first, (*it).second * (*it2).second );		//re-use
			IVTmp.second.Add(	(*it).first, (*it).second * (1-(*it2).second) );	//waste
		}
		else IVTmp.second.Add(	(*it).first, (*it).second );	//waste
	}
DBGL;
	return IVTmp;
}









