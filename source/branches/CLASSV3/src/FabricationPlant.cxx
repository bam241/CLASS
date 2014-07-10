#include "FabricationPlant.hxx"

#include "Storage.hxx"
#include "Reactor.hxx"
#include "EvolutionData.hxx"
#include "DecayDataBank.hxx"
#include "PhysicModels.hxx"
#include "IsotopicVector.hxx"
#include "Scenario.hxx"
#include "CLASSLogger.hxx"
#include "CLASSConstante.hxx"




#include "TMatrixT.h"

#include <sstream>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

	//________________________________________________________________________
	//________________________________________________________________________
	//________________________________________________________________________
	//
	//		FabricationPlant
	//
	//
	//
	//
	//________________________________________________________________________
	//________________________________________________________________________
ClassImp(FabricationPlant)



FabricationPlant::FabricationPlant():CLASSFacility(16)
{
	SetName("F_FabricationPLant.");
	
	fReUsable = 0;
}


FabricationPlant::FabricationPlant(CLASSLogger* log, double fabricationtime):CLASSFacility(log, fabricationtime, 16)
{
	SetName("F_FabricationPLant.");

	fFiFo = false;
	fSubstitutionFuel = false;


	INFO	<< " A FabricationPlant has been define :" << endl;
	INFO	<< "\t Chronological Stock Priority has been set! "<< endl;
	INFO	<< "\t Fabrication time set to \t " << (double)(GetCycleTime()/3600/24/365.25) << " year" << endl << endl;
	
	


}



	//________________________________________________________________________
FabricationPlant::~FabricationPlant()
{
	
	
}


	//________________________________________________________________________
void	FabricationPlant::SetSeparartionEfficiencyIV(ZAI zai, double factor)
{

	pair<map<ZAI, double>::iterator, bool> IResult;
	if(factor > 1) factor = 1;
	
	if(factor > 0)
	{
		IResult =  fSeparationLostFraction.GetIsotopicQuantity().insert( pair<ZAI ,double>(zai, 1 - factor));
		if(!IResult.second)
			IResult.first->second = 1 - factor;
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
DBGL
	IsotopicVector fInsideIV;


	map<int ,cSecond >::iterator it;
	for( it = fReactorNextStep.begin(); it!= fReactorNextStep.end(); it++ )
	{
		if( t + GetCycleTime() >= GetParc()->GetReactor()[ (*it).first ]->GetCreationTime()
		   && t + GetCycleTime() < GetParc()->GetReactor()[ (*it).first ]->GetCreationTime() + GetParc()->GetReactor()[ (*it).first ]->GetLifeTime())
		{
			if( (*it).second == t )
			{
#pragma omp critical(FuelBuild)
				{BuildFuelForReactor( (*it).first );}
				(*it).second += GetParc()->GetReactor()[ (*it).first ]->GetCycleTime();
			}
			else if ( (*it).second - GetParc()->GetReactor()[ (*it).first ]->GetCycleTime() + GetCycleTime() > t )
			{
				map<int ,IsotopicVector >::iterator it2 = fReactorFuturIV.find( (*it).first );
				if (it2 != fReactorFuturIV.end())
					(*it2).second = GetDecay((*it2).second, t - fInternalTime );

				fInsideIV += (*it2).second;
			}
		}
	}
	
	
DBGL
}


	//________________________________________________________________________
void FabricationPlant::BuildFuelForReactor(int ReactorId)
{
	DBGV( " IN"<< "in ");
	if(fFissileStorage.size() == 0)
	{
		ERROR << " One need at least one Fissile storage to build fuel " << endl;
		ERROR << " Use AddFissileStorage to add a stock to provide fissil material... " << endl;
		ERROR << " One need at least one Fissile storage to build fuel " << endl;
		exit(1);
	}



	double R_HM_Mass	= GetParc()->GetReactor()[ ReactorId ]->GetHeavyMetalMass();
	double R_BU		= GetParc()->GetReactor()[ ReactorId ]->GetBurnUp();
	double R_CycleTime	= GetParc()->GetReactor()[ ReactorId ]->GetCycleTime();
	double R_Power		= GetParc()->GetReactor()[ ReactorId ]->GetPower();

	PhysicModels* FuelType = GetParc()->GetReactor()[ReactorId]->GetFuelType();

	fFissileList = FuelType->GetEquivalenceModel()->GetFissileList();
	BuildFissileArray();


	fFertileList = FuelType->GetEquivalenceModel()->GetFertileList();


	if(fFertileStorage.size() != 0)			// If the fertile need to be taken in stock
		BuildFertileArray();
	else						// if their is not stock and the fertile come from outside of the park
	{
		fFertileArray.push_back( fFertileList / fFertileList.GetTotalMass() * R_HM_Mass );
	}

	vector<double> LambdaArray =  FuelType->GetEquivalenceModel()->BuildFuel(R_BU, R_HM_Mass, fFissileArray, fFertileArray);


	if(LambdaArray[0] != -1)
	{
		IsotopicVector IV = BuildFuelFromEqModel(LambdaArray);
		EvolutionData EvolDB = FuelType->GenerateEvolutionData( GetDecay(IV,fCycleTime), R_CycleTime, R_Power);

		{
			pair<map<int, IsotopicVector>::iterator, bool> IResult;
			IResult = fReactorFuturIV.insert( pair<int, IsotopicVector>(ReactorId, IV) );
			if(!IResult.second)
			IResult.first->second = IV;
		}
		{
			pair<map<int, EvolutionData>::iterator, bool> IResult;
			IResult = fReactorFuturDB.insert( pair<int, EvolutionData>(ReactorId,EvolDB) );
			if(!IResult.second)
			IResult.first->second = EvolDB;
		}
		fInsideIV += IV;
		AddCumulativeIVIn(IV);

		return;
	}
	else
	{

		if (!fSubstitutionFuel)
		{
			{
				EvolutionData EmptyDB;
				pair<map<int, EvolutionData>::iterator, bool> IResult;
				IResult = fReactorFuturDB.insert( pair<int, EvolutionData>(ReactorId,EmptyDB) );
				if(!IResult.second)
					IResult.first->second = EmptyDB;
			}
			{
				IsotopicVector EmptyIV;
				pair<map<int, IsotopicVector>::iterator, bool> IResult;
				IResult = fReactorFuturIV.insert( pair<int, IsotopicVector>(ReactorId,EmptyIV) );
				if(!IResult.second)
					IResult.first->second = EmptyIV;
			}
		}
		else
		{

			IsotopicVector IV = fSubstitutionEvolutionData.GetIsotopicVectorAt(0);
			EvolutionData evolutiondb = fSubstitutionEvolutionData * R_HM_Mass / IV.GetTotalMass();

			IV = IV* R_HM_Mass / IV.GetTotalMass();
			{
				pair<map<int, IsotopicVector>::iterator, bool> IResult;
				IResult = fReactorFuturIV.insert( pair<int, IsotopicVector>(ReactorId, IV) );
				if(!IResult.second)
				IResult.first->second = IV;
			}
			{
				pair<map<int, EvolutionData>::iterator, bool> IResult;
				IResult = fReactorFuturDB.insert( pair<int, EvolutionData>(ReactorId,evolutiondb) );
				if(!IResult.second)
				IResult.first->second = evolutiondb;
			}
			GetParc()->AddGod( IV );
			fInsideIV += IV;
			AddCumulativeIVIn(IV);



		}
		return;
	}
DBGL
}



void FabricationPlant::BuildFissileArray()
{

	for(int i = 0; i < (int)fFissileStorage.size(); i++)
	{
		vector<IsotopicVector> IVArray = fFissileStorage[i]->GetIVArray();

		for(int j = 0; j < (int)IVArray.size(); j++)
		{

			IsotopicVector SeparatedIV = Separation(IVArray[j], fFissileList).first;
			IsotopicVector CooledSeparatedIV = GetDecay( SeparatedIV , GetCycleTime());

			fFissileArray.push_back( CooledSeparatedIV );
			fFissileArrayAdress.push_back( pair<int,int>(i,j) );
			fFissileArrayTime.push_back(fFissileStorage[i]->GetIVArrayArrivalTime()[j]);
		}

	}

	SortArray(0);
}


void FabricationPlant::BuildFertileArray()
{


	for(int i = 0; i < (int)fFertileStorage.size(); i++)
	{
		vector<IsotopicVector> IVArray = fFertileStorage[i]->GetIVArray();
		for(int j = 0; j < (int)IVArray.size(); j++)
		{

			IsotopicVector SeparatedIV = Separation(IVArray[j], fFertileList).first;
			IsotopicVector CooledSeparatedIV = GetDecay( SeparatedIV , GetCycleTime());


			fFertileArray.push_back( CooledSeparatedIV );
			fFertileArrayAdress.push_back( pair<int,int>(i,j) );
			fFertileArrayTime.push_back(fFertileStorage[i]->GetIVArrayArrivalTime()[j]);
		}

	}

	SortArray(1);

}

void FabricationPlant::SortArray(int i)
{


	vector<IsotopicVector>	IVArray;
	vector<cSecond>		TimeArray;
	vector< pair<int,int> >	AdressArray;

	if(i==0) //Fissile
	{
		IVArray		= fFissileArray;
		TimeArray	= fFissileArrayTime;
		AdressArray	= fFissileArrayAdress;
	}
	else if (i==1) //Fertile
	{
		IVArray		= fFertileArray;
		TimeArray	= fFertileArrayTime;
		AdressArray	= fFertileArrayAdress;

	}

	if(fFiFo)
	{
		for(int j = 0; j < (int)TimeArray.size(); j++)
		{
			for (int k = j+1; k < (int)TimeArray.size(); k++)
			{
				cSecond time_tmp = TimeArray[j];
				pair<int,int> Adress_tmp = AdressArray[j];
				IsotopicVector IV_tmp = IVArray[j];

				if(time_tmp > TimeArray[k])
				{
					TimeArray[j] = TimeArray[k];
					TimeArray[k] = time_tmp;

					AdressArray[j] = AdressArray[k];
					AdressArray[k] = Adress_tmp;

					IVArray[j] = IVArray[k];
					IVArray[k] = IV_tmp;
				}

			}
		}
	}
	else
	{
		for(int j = 0; j < (int)fFissileArrayTime.size(); j++)
		{
			for (int k = j+1; k < (int)TimeArray.size(); k++)
			{
				cSecond time_tmp = TimeArray[j];
				pair<int,int> Adress_tmp = AdressArray[j];
				IsotopicVector IV_tmp = IVArray[j];

				if(time_tmp < TimeArray[k])
				{
					TimeArray[j] = TimeArray[k];
					TimeArray[k] = time_tmp;

					AdressArray[j] = AdressArray[k];
					AdressArray[k] = Adress_tmp;

					IVArray[j] = IVArray[k];
					IVArray[k] = IV_tmp;
				}
				
			}
		}
	}


	if(i==0) //Fissile
	{
		fFissileArray		= IVArray;
		fFissileArrayTime	= TimeArray;
		fFissileArrayAdress	= AdressArray;
	}
	else if (i==1) //Fertile
	{
		fFertileArray = IVArray;
		fFertileArrayTime = TimeArray;
		fFertileArrayAdress = AdressArray;

	}


}




//________________________________________________________________________
void	FabricationPlant::SetSubstitutionFuel(EvolutionData fuel)
{
	
	fSubstitutionFuel = true;
	map<ZAI ,double>::iterator it;
	map<ZAI ,double> isotopicquantity = fuel.GetIsotopicVectorAt(0.).GetActinidesComposition().GetIsotopicQuantity();
	double M0 = 0;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
		M0 += (*it).second*cZAIMass.fZAIMass.find( (*it).first )->second/AVOGADRO*1e-6;
	fSubstitutionEvolutionData = fuel / M0;

}


	//________________________________________________________________________
	//_____________________________ Reactor & DB _____________________________
	//________________________________________________________________________
	//________________________________________________________________________
void FabricationPlant::TakeReactorFuel(int Id)
{
	
	
	IsotopicVector IV;
	map<int ,IsotopicVector >::iterator it2 = fReactorFuturIV.find( Id );
	AddCumulativeIVOut(it2->second);
	fInsideIV -= (*it2).second;


	if (it2 != fReactorFuturIV.end())
		(*it2).second = IV;


	map< int,EvolutionData >::iterator it = fReactorFuturDB.find(Id);
	(*it).second = EvolutionData();

}

//________________________________________________________________________
EvolutionData FabricationPlant::GetReactorEvolutionDB(int ReactorId)
{

	map< int,EvolutionData >::iterator it = fReactorFuturDB.find(ReactorId);
	return (*it).second;
}
	//________________________________________________________________________
	//_______________________________ Storage ________________________________
	//________________________________________________________________________
IsotopicVector FabricationPlant::BuildFuelFromEqModel(vector<double> LambdaArray)
{
DBGL
	IsotopicVector BuildedFuel;
	IsotopicVector Lost;

	for(int i = 0; i < (int)fFissileArray.size(); i++)
	{
		if(LambdaArray[i] != 0)
		{
			int Stor_N = fFissileArrayAdress[i].first;
			int IV_N = fFissileArrayAdress[i].second;

			pair<IsotopicVector, IsotopicVector> Separated_Lost;
			Separated_Lost = Separation( fFissileStorage[Stor_N]->GetIVArray()[IV_N]*LambdaArray[i], fFissileList);
			BuildedFuel += Separated_Lost.first;
			Lost += Separated_Lost.second;
		}
	}

	if(fFertileStorage.size() != 0)
	{
		for(int i = fFissileArray.size(); i < (int)(fFertileArray.size()+fFissileArray.size()); i++)
		{
			if(LambdaArray[i] != 0)
			{
				int Stor_N = fFertileArrayAdress[i].first;
				int IV_N = fFertileArrayAdress[i].second;

				pair<IsotopicVector, IsotopicVector> Separated_Lost;
				Separated_Lost = Separation( fFertileStorage[Stor_N]->GetIVArray()[IV_N]*LambdaArray[i], fFissileList);
				BuildedFuel += Separated_Lost.first;
				Lost += Separated_Lost.second;
			}
		}
	}
	else
		BuildedFuel += fFertileArray[0]*LambdaArray.back();

	DumpStock(LambdaArray);

DBGL
	return BuildedFuel;
}


	//________________________________________________________________________
void FabricationPlant::DumpStock(vector<double> LambdaArray)
{
DBGL

	for(int i = 0; i < (int)fFissileArray.size(); i++)
	{
		if(LambdaArray[i] != 0)
		{
			int Stor_N = fFissileArrayAdress[i].first;
			int IV_N = fFissileArrayAdress[i].second;
			fFissileStorage[Stor_N]->TakeFractionFromStock( IV_N, LambdaArray[i] );
		}
	}
	if(fFertileStorage.size() != 0)
	{
		for(int i = fFissileArray.size(); i < (int)(fFertileArray.size()+fFissileArray.size()); i++)
		{
			if(LambdaArray[i] != 0)
			{
				int Stor_N = fFertileArrayAdress[i].first;
				int IV_N = fFertileArrayAdress[i].second;

				fFertileStorage[Stor_N]->TakeFractionFromStock( IV_N, LambdaArray[i] );
			}
		}
	}
	else
		GetParc()->AddGod( fFertileArray[0]*LambdaArray.back() );




	//Clear the Building Array (Fissile and Fertile)
	fFissileArray.clear();
	fFissileArrayTime.clear();
	fFissileArrayAdress.clear();
	fFertileArray.clear();
	fFertileArrayTime.clear();
	fFertileArrayAdress.clear();

	fFertileList = fFissileList = IsotopicVector();

DBGL
}

	//________________________________________________________________________
pair<IsotopicVector, IsotopicVector> FabricationPlant::Separation(IsotopicVector isotopicvector, IsotopicVector ExtractedList)
{
	
		//[0] = re-use ; [1] = waste
	IsotopicVector LostPart  = isotopicvector.GetThisComposition(ExtractedList) * fSeparationLostFraction;

	IsotopicVector SeparatedPart  = isotopicvector.GetThisComposition(ExtractedList) - LostPart;
	LostPart = isotopicvector - SeparatedPart;


	return pair<IsotopicVector, IsotopicVector> (SeparatedPart, LostPart);
}



//________________________________________________________________________
//	Get Decay
//________________________________________________________________________
IsotopicVector FabricationPlant::GetDecay(IsotopicVector isotopicvector, cSecond t)
{

	IsotopicVector IV;

	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	{
		if((*it).second > 0)
		{
 			IsotopicVector ivtmp = fDecayDataBase->Evolution(it->first, t) * (*it).second ;
			IV += ivtmp;
		}
	}

	return IV;
	
}






