#include "FabricationPlant.hxx"

#include "Storage.hxx"
#include "FabricationPlant.hxx"
#include "Reactor.hxx"
#include "EvolutionData.hxx"
#include "DataBank.hxx"
#include "IsotopicVector.hxx"
#include "CLASS.hxx"
#include "Defines.hxx"
#include "LogFile.hxx"




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
	

	
}

FabricationPlant::FabricationPlant(LogFile* log)
{
	
	fLog = log;
	fChronologicalTimePriority = false;
	fFabricationTime = -1;
	fUpdateReferenceDBatEachStep = false;
	fSubstitutionFuel = false;
	
		// Warning
	
	cout	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	cout	<< "\t Chronological Stock Priority set! "<< endl << endl;
	cout	<< "!!WARNING!! !!!FabricationPlant!!! You need to set the different stock manually as well as the Fabrication Time Manualy !! " << endl;
	
	fLog->fLog	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	fLog->fLog	<< "\t Chronological Stock Priority set! "<< endl << endl;
	fLog->fLog	<< "!!WARNING!! !!!FabricationPlant!!! You need to set the different stock manually as well as the Fabrication Time Manualy !! " << endl;
	
	
}

FabricationPlant::FabricationPlant(LogFile* log, Storage* storage, Storage* reusable, double fabircationtime)
{
	
	fLog = log;
	
	fChronologicalTimePriority = false;
	fUpdateReferenceDBatEachStep = false;
	fSubstitutionFuel = false;
	
	fFabricationTime = (cSecond)fabircationtime;
	fStorage = storage;
	fReUsable = reusable;
	
	
	cout	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	cout	<< "\t Chronological Stock Priority has been set! "<< endl;
	cout	<< "\t Fabrication time set to \t " << (double)(fFabricationTime/3600/24/365.25) << " year" << endl << endl;
	
	fLog->fLog	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	fLog->fLog	<< "\t Chronological Stock Priority has been set! "<< endl;
	fLog->fLog	<< "\t Fabrication time set to \t " << (double)(fFabricationTime/3600/24/365.25) << " year" << endl << endl;
	
	
}



	//________________________________________________________________________
FabricationPlant::~FabricationPlant()
{
	
	
}


	//________________________________________________________________________
void	FabricationPlant::AddValorisableIV(ZAI zai, double factor)
{
	
	pair<map<ZAI, double>::iterator, bool> IResult;
	if(factor > 1) factor = 1;
	
	if(factor > 0)
	{
		IResult = fValorisableIV.insert( pair<ZAI ,double>(zai, factor));
		if(IResult.second == false)
			IResult.first->second = factor;
	}
	
}



	//________________________________________________________________________
	//_______________________________ Evolution ______________________________
	//________________________________________________________________________
void FabricationPlant::Evolution(cSecond t)
{
	
		// Check if the FabricationPlant has been created ...
	if(t == fInternalTime && t != 0) return;
		// Make the evolution for the FabricationPlant ...
	FabricationPlantEvolution(t);
		// ... And Finaly update the AbsoluteInternalTime
	fInternalTime = t;
	
}

	//________________________________________________________________________
void FabricationPlant::FabricationPlantEvolution(cSecond t)
{
	
	
	
	map<int ,cSecond >::iterator it;
	for( it = fReactorNextStep.begin(); it!= fReactorNextStep.end(); it++ )
	{
		if( t + fFabricationTime >= fParc->GetReactor()[ (*it).first ]->GetCreationTime()
		   && t + fFabricationTime < fParc->GetReactor()[ (*it).first ]->GetCreationTime() + fParc->GetReactor()[ (*it).first ]->GetLifeTime())
		{
			if( (*it).second == t )
			{
				BuildFuelForReactor( (*it).first );
				(*it).second += fParc->GetReactor()[ (*it).first ]->GetCycleTime();
			}
			else if ( (*it).second - fParc->GetReactor()[ (*it).first ]->GetCycleTime() + fFabricationTime > t )
			{
				map<int ,IsotopicVector >::iterator it2 = fReactorFuturIV.find( (*it).first );
				if (it2 != fReactorFuturIV.end())
					(*it2).second = GetDecay((*it2).second, t - fInternalTime );
			}
		}
	}
	
	
	
}


	//________________________________________________________________________
void FabricationPlant::BuildFuelForReactor(int ReactorId)
{
	
	
	DataBank<IsotopicVector>* FuelType = fParc->GetReactor()[ReactorId]->GetFuelType();
	
	
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
	while(!FuelBuild)
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
		
		if( stock.GetZAIIsotopicQuantity(ZAI(-1,-1,-1)) == 1 ) // Not enought stock to build the needed fuel
		{
			if (!fSubstitutionFuel)
			{
				{
					EvolutionData evolutiondb;
					pair<map<int, EvolutionData>::iterator, bool> IResult;
					IResult = fReactorFuturDB.insert( pair<int, EvolutionData>(ReactorId,evolutiondb) );
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
			}
			else
			{
				{
					EvolutionData evolutiondb = fSubstitutionEvolutionData* HMmass;
					pair<map<int, EvolutionData>::iterator, bool> IResult;
					IResult = fReactorFuturDB.insert( pair<int, EvolutionData>(ReactorId,evolutiondb) );
					if(IResult.second == false)
						IResult.first->second = evolutiondb;
				}
				{
					IsotopicVector IV = fSubstitutionEvolutionData.GetIsotopicVectorAt(0)* HMmass;
					pair<map<int, IsotopicVector>::iterator, bool> IResult;
					IResult = fReactorFuturIV.insert( pair<int, IsotopicVector>(ReactorId, IV) );
					if(IResult.second == false)
						IResult.first->second = IV;
				}
				
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
						Sum_AlphaI_nPuI += FuelType->GetFuelParameter()[(*it).first.A() -237]*(*it).second;
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
						Sum_AlphaI_nPuI0 += FuelType->GetFuelParameter()[(*it).first.A() -237]*(*it).second;
					}
			}
			
			double StockFactionToUse = 0;
			
			double NT = HMmass*1e6 * Na / (ZAImass.find( ZAI(92,238,0) )->second*0.997 + ZAImass.find( ZAI(92,235,0) )->second*0.003 );
			
			double N1 = (BU - FuelType->GetFuelParameter()[6]) * NT;
			double N2 = -Sum_AlphaI_nPuI0;
			double N3 = -FuelType->GetFuelParameter()[0] * Na / (ZAImass.find( ZAI(92,238,0) )->second*0.997 + ZAImass.find( ZAI(92,235,0) )->second*0.003 ) * (HMmass*1e6 - MPu_0*1e6);
			
			double D1 = Sum_AlphaI_nPuI;
			double D2 = -FuelType->GetFuelParameter()[0] * MPu_1*1e6 * Na / (ZAImass.find( ZAI(92,238,0) )->second*0.997 + ZAImass.find( ZAI(92,235,0) )->second*0.003 ) ;
			
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
				ZAI U5 = ZAI(92,235,0);
				double U8_Quantity = (HMmass - (MPu_0+StockFactionToUse*MPu_1 ))/(ZAImass.find( ZAI(92,238,0) )->second*0.997 + ZAImass.find( ZAI(92,235,0) )->second*0.003 )*Na/1e-6;
				
				fParc->AddGodIncome( U8, U8_Quantity*0.997 );
				fParc->AddGodIncome( U5, U8_Quantity*0.003 );
				
				for(int i = (int)fFractionToTake.size()-1; i >= 0; i--)
				{
					IVBeginCycle += fStorage->GetStock()[fFractionToTake[i].first].GetSpeciesComposition(94)*( fFractionToTake[i].second );					
					fReUsable->AddToStock(fStorage->GetStock()[fFractionToTake[i].first]*(fFractionToTake[i].second)
							      - fStorage->GetStock()[fFractionToTake[i].first].GetSpeciesComposition(94)*(fFractionToTake[i].second));
					
					fStorage->TakeFractionFromStock(fFractionToTake[i].first,fFractionToTake[i].second);			
					
				}
				fFractionToTake.clear();
				
				IVBeginCycle += U8_Quantity*U8*0.997 + U8_Quantity*U5*0.003;
				EvolutionData evolutiondb = BuildEvolutiveDB(ReactorId, IVBeginCycle);
				
				{
					pair<map<int, EvolutionData>::iterator, bool> IResult;
					IResult = fReactorFuturDB.insert( pair<int, EvolutionData>(ReactorId,evolutiondb) );
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
	
	
}

void	FabricationPlant::SetSubstitutionFuel(EvolutionData fuel)
{
	
	fSubstitutionFuel = true;
	fSubstitutionEvolutionData = fuel / fuel.GetHMMass();
	
}


	//________________________________________________________________________
	//_________________________________ Decay ________________________________
	//________________________________________________________________________
IsotopicVector FabricationPlant::GetDecay(IsotopicVector isotopicvector, cSecond t)
{
	
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
	
	return IV;
}


	//________________________________________________________________________
	//_____________________________ Reactor & DB _____________________________
	//________________________________________________________________________
EvolutionData FabricationPlant::BuildEvolutiveDB(int ReactorId,IsotopicVector isotopicvector)
{
	
	DataBank<IsotopicVector>* evolutiondb = fParc->GetReactor()[ReactorId]->GetFuelType();
	
	isotopicvector = GetDecay(isotopicvector, fFabricationTime);
	
	EvolutionData EvolBuild;
	/*
	if( fUpdateReferenceDBatEachStep == true )
	{
		EvolutionData EvolBuild = evolutiondb->GenerateDB(isotopicvector,
								  fParc->GetReactor()[ReactorId]->GetCycleTime(),
								  fParc->GetReactor()[ReactorId]->GetPower());
	}
	else
	{
		map<double, EvolutionData> distances = evolutiondb->GetDistancesTo(isotopicvector);
		EvolBuild = distances.begin()->second.GenerateDBFor(isotopicvector);
		
	}
	 */
	EvolBuild = evolutiondb->GenerateEvolutionData(isotopicvector,
					    fParc->GetReactor()[ReactorId]->GetCycleTime(),
					    fParc->GetReactor()[ReactorId]->GetPower());
	 
	return EvolBuild;
	
}

	//________________________________________________________________________
void FabricationPlant::TakeReactorFuel(int Id)
{
	
	
	IsotopicVector IV;
	map<int ,IsotopicVector >::iterator it2 = fReactorFuturIV.find( Id );
	if (it2 != fReactorFuturIV.end())
		(*it2).second = IV;
	
}

	//________________________________________________________________________
EvolutionData FabricationPlant::GetReactorEvolutionDB(int ReactorId)
{
	
	map< int,EvolutionData >::iterator it = fReactorFuturDB.find(ReactorId);
	return (*it).second;
	
}

IsotopicVector FabricationPlant::GetFullFabrication()
{
	
	IsotopicVector tmp = 0*ZAI(0,0,0);
	
	map<int, IsotopicVector > reactorNextStep = fReactorFuturIV;
	map<int, IsotopicVector >::iterator it;
	for ( it = reactorNextStep.begin(); it != reactorNextStep.end(); it++)
		tmp += (*it).second;

	return tmp;
	
}

	//________________________________________________________________________
	//_______________________________ Storage ________________________________
	//________________________________________________________________________

IsotopicVector FabricationPlant::GetStockToRecycle()
{
	
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
	
}

	//________________________________________________________________________
void FabricationPlant::RecycleStock(double fraction)
{
	
	fFractionToTake.back().second = fraction;
	
}



	//________________________________________________________________________
void FabricationPlant::DumpStock()
{
	
	
	
	
}

	//________________________________________________________________________
pair<IsotopicVector, IsotopicVector> FabricationPlant::Separation(IsotopicVector isotopicvector)
{
	
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
	
	return IVTmp;
}









