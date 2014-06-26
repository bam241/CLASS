#include "FabricationPlant.hxx"

#include "Storage.hxx"
#include "Reactor.hxx"
#include "EvolutionData.hxx"
#include "DecayDataBank.hxx"
#include "PhysicModels.hxx"
#include "IsotopicVector.hxx"
#include "CLASS.hxx"
#include "CLASSHeaders.hxx"
#include "LogFile.hxx"




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
template <class T>  T random(T a, T b) //peak random numebr between a and b
{
	double range = pow(2., 31);
	srand(time(NULL)); //initialize the srand
	return (T)a + (T)(b-a)*rand()/range;
}

ClassImp(FabricationPlant)



FabricationPlant::FabricationPlant():CLASSFacility()
{
	SetFacilityType(16);
	SetName("F_FabricationPLant.");
	
	fStorage = 0;
	fReUsable = 0;
}


FabricationPlant::FabricationPlant(LogFile* log, double fabricationtime):CLASSFacility(log, fabricationtime)
{
	SetFacilityType(16);
	SetName("F_FabricationPLant.");

	fChronologicalTimePriority = false;
	fUpdateReferenceDBatEachStep = false;
	fSubstitutionFuel = false;


	cout	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	cout	<< "\t Chronological Stock Priority has been set! "<< endl;
	cout	<< "\t Fabrication time set to \t " << (double)(GetCycleTime()/3600/24/365.25) << " year" << endl << endl;
	
	GetLog()->fLog	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	GetLog()->fLog	<< "\t Chronological Stock Priority has been set! "<< endl;
	GetLog()->fLog	<< "\t Fabrication time set to \t " << (double)(GetCycleTime()/3600/24/365.25) << " year" << endl << endl;
	
	


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
		if(!IResult.second)
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
	
	IsotopicVector EmptyIV;
	fInsideIV = EmptyIV;


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
	
	
	
}


	//________________________________________________________________________

void FabricationPlant::BuildFuelForReactor(int ReactorId)
{

	if(fFissileStorage.size() == 0)
	{
		cout << "!!Error!! !!!FabricationPlant!!! One need at least one Fissile storage to build fuel " << endl;
		cout << "!!Error!! !!!FabricationPlant!!! use AddFissileStorage to add a stock to provide fissil material... " << endl;
		GetLog()->fLog << "!!Error!! !!!FabricationPlant!!! One need at least one Fissile storage to build fuel " << endl;
		exit(1);
	}



//	PhysicModels* FuelType = GetParc()->GetReactor()[ReactorId]->GetFuelType();
}

//________________________________________________________________________
void	FabricationPlant::SetSubstitutionFuel(EvolutionData fuel)
{
	
	fSubstitutionFuel = true;
	double Na = 6.02214129e23;	//N Avogadro
	map<ZAI ,double>::iterator it;
	map<ZAI ,double> isotopicquantity = fuel.GetIsotopicVectorAt(0.).GetActinidesComposition().GetIsotopicQuantity();
	double M0 = 0;
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
		M0 += (*it).second*cZAIMass.fZAIMass.find( (*it).first )->second/Na*1e-6;
	fSubstitutionEvolutionData = fuel / M0;

}


	//________________________________________________________________________
	//_____________________________ Reactor & DB _____________________________
	//________________________________________________________________________
EvolutionData FabricationPlant::BuildEvolutiveDB(int ReactorId,IsotopicVector isotopicvector)
{
	
	FuelDataBank* evolutiondb = GetParc()->GetReactor()[ReactorId]->GetFuelType();
	
	isotopicvector = GetDecay(isotopicvector, GetCycleTime());
	
	EvolutionData EvolBuild;

#pragma omp single
	{
		EvolBuild = evolutiondb->GenerateEvolutionData(isotopicvector, GetParc()->GetReactor()[ReactorId]->GetCycleTime(), GetParc()->GetReactor()[ReactorId]->GetPower());
	}
	return EvolBuild;
	
}

	//________________________________________________________________________
void FabricationPlant::TakeReactorFuel(int Id)
{
	
	
	IsotopicVector IV;
	map<int ,IsotopicVector >::iterator it2 = fReactorFuturIV.find( Id );
	AddCumulativeIVOut(it2->second);
	fInsideIV -= (*it2).second;


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
	
	if(fChronologicalTimePriority )
		IdToTake = (int)( fFractionToTake.size() );
	else
		IdToTake = (int)( fStorage->GetIVArray().size() -1 - fFractionToTake.size() );

	if(0 <= IdToTake && IdToTake < (int)fStorage->GetIVArray().size())
	{
		NextStock = fStorage->GetIVArray()[IdToTake];
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

	for(int i = (int)fFractionToTake.size()-1; i >= 0; i--)
	{

		pair<IsotopicVector, IsotopicVector> SeparatedIV;
		IsotopicVector IV_in_Stock_i = fStorage->GetIVArray()[fFractionToTake[i].first];
		double Fraction_Taken_from_Stock_i = fFractionToTake[i].second ;

		SeparatedIV = Separation( Fraction_Taken_from_Stock_i * ( IV_in_Stock_i - IV_in_Stock_i.GetSpeciesComposition(94) ) );

		fReUsable->AddIV(SeparatedIV.first);
		GetParc()->AddWaste(SeparatedIV.second);




		fStorage->TakeFractionFromStock(fFractionToTake[i].first,fFractionToTake[i].second);


	}
	fFractionToTake.clear();




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
		else
		{
			IVTmp.first.Add(	(*it).first, (*it).second );	//re-use
//			IVTmp.second.Add(	(*it).first, (*it).second * (1-(*it2).second) );	//waste -> Empty
		}
	}
	
	return IVTmp;
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






