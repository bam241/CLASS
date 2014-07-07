#include "PWR_THU_FabricationPlant.hxx"

#include "Storage.hxx"
#include "Reactor.hxx"
#include "EvolutionData.hxx"
#include "FuelDataBank.hxx"
#include "DecayDataBank.hxx"
#include "IsotopicVector.hxx"
#include "Scenario.hxx"
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
	//		MyFabricationPlant
	//	a FP specific to one type of Reactor and one type of Fuel
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

PWR_THU_FabricationPlant::PWR_THU_FabricationPlant():FabricationPlant()
{

}

PWR_THU_FabricationPlant::PWR_THU_FabricationPlant(LogFile* log)
{
	SetFacilityType(16);

	SetLog(log);
	fFiFo = false;
	SetCycleTime(-1);
	fUpdateReferenceDBatEachStep = false;
	fSubstitutionFuel = false;
	fStorage = 0;
	fReUsable = 0;

	cout	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	cout	<< "\t Chronological Stock Priority set! "<< endl << endl;
	cout	<< "!!WARNING!! !!!FabricationPlant!!! You need to set the different stock manually as well as the Fabrication Time Manualy !! " << endl;
	GetLog()->fLog	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	GetLog()->fLog	<< "\t Chronological Stock Priority set! "<< endl << endl;
	GetLog()->fLog	<< "!!WARNING!! !!!FabricationPlant!!! You need to set the different stock manually as well as the Fabrication Time Manualy !! " << endl;
	
	

}

PWR_THU_FabricationPlant::PWR_THU_FabricationPlant(LogFile* log, Storage* storage, Storage* reusable, double fabircationtime)
{
	SetFacilityType(16);

	SetLog(log);
	
	fFiFo = false;
	fUpdateReferenceDBatEachStep = false;
	fSubstitutionFuel = false;


	SetCycleTime((cSecond)fabircationtime );
	fStorage = storage;
	fReUsable = reusable;
	
	
	cout	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	cout	<< "\t Chronological Stock Priority has been set! "<< endl;
	cout	<< "\t Fabrication time set to \t " << (double)(GetCycleTime()/3600/24/365.25) << " year" << endl << endl;
	
	GetLog()->fLog	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	GetLog()->fLog	<< "\t Chronological Stock Priority has been set! "<< endl;
	GetLog()->fLog	<< "\t Fabrication time set to \t " << (double)(GetCycleTime()/3600/24/365.25) << " year" << endl << endl;
	
	


}



	//________________________________________________________________________
PWR_THU_FabricationPlant::~PWR_THU_FabricationPlant()
{
	
	
}





void PWR_THU_FabricationPlant::BuildFuelForReactor(int ReactorId)
{
	//cout<<"INFO : This is a specific FabricationPlant"<<endl<<"GOOD JOB!!!"<<endl;
	FuelDataBank* FuelType = GetParc()->GetReactor()[ReactorId]->GetFuelType();
	string ReactorType ="PWR";	
	if(FuelType->GetFuelType() != "THU" || ReactorType !="PWR")//Check if the reactor is the right type and use the right type of fuel
	{
		cout << "!!Bad Trouble!! !!!FabricationPlant!!! Try to do MOX with a not MOXed DB "<< endl;
		GetLog()->fLog << "!!Bad Trouble!! !!!FabricationPlant!!! Try to do MOX with a not MOXed DB" << endl;
		exit (1);
	}	

	double HMmass = GetParc()->GetReactor()[ReactorId]->GetHeavyMetalMass();
	double BU = GetParc()->GetReactor()[ReactorId]->GetBurnUp();
	IsotopicVector FullUsedStock;
	IsotopicVector stock;
	
	bool FuelBuild = false;
	while(!FuelBuild)
	{
		double nU_0 = 0;
		double MU_0 = 0;
		{
			map<ZAI ,double>::iterator it;
			
			map<ZAI ,double> isotopicquantity = GetDecay( FullUsedStock , GetCycleTime()).GetSpeciesComposition(92).GetIsotopicQuantity();
			for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
				nU_0 += (*it).second;
			
			isotopicquantity = FullUsedStock.GetSpeciesComposition(92).GetIsotopicQuantity();
			for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
				MU_0 += (*it).second*cZAIMass.fZAIMass.find( (*it).first )->second/AVOGADRO*1e-6;
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
					if(!IResult.second)
						IResult.first->second = evolutiondb;
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
				{
					EvolutionData evolutiondb = fSubstitutionEvolutionData* HMmass;
					pair<map<int, EvolutionData>::iterator, bool> IResult;
					IResult = fReactorFuturDB.insert( pair<int, EvolutionData>(ReactorId,evolutiondb) );
					if(!IResult.second)
						IResult.first->second = evolutiondb;
				}
				{
					IsotopicVector IV = fSubstitutionEvolutionData.GetIsotopicVectorAt(0)* HMmass;
					pair<map<int, IsotopicVector>::iterator, bool> IResult;
					IResult = fReactorFuturIV.insert( pair<int, IsotopicVector>(ReactorId, IV) );
					if(!IResult.second)
						IResult.first->second = IV;
				}
				
			}
			FuelBuild = true;
			fFractionToTake.clear();
		}
		else
		{
			double nU_1 = 0;
			double MU_1 = 0;
			double Sum_AlphaI_nUI = 0;
			double Sum_AlphaI_nUI0 = 0;
			{
				map<ZAI ,double>::iterator it;
				map<ZAI ,double> isotopicquantity = GetDecay( stock , GetCycleTime()).GetSpeciesComposition(92).GetIsotopicQuantity();
				
				for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
				{
					if ((*it).first.A() >= 232 && (*it).first.A() <= 236)
					{
						nU_1 += (*it).second;
						Sum_AlphaI_nUI += FuelType->GetFuelParameter()[(*it).first.A() -231]*(*it).second;
					}
				}
				
				isotopicquantity = stock.GetSpeciesComposition(92).GetIsotopicQuantity();
				for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
					if ((*it).first.A() >= 232 && (*it).first.A() <= 236)
					{
						MU_1 += (*it).second * (cZAIMass.fZAIMass.find( (*it).first )->second)/AVOGADRO*1e-6;
					}
				
				isotopicquantity = GetDecay( FullUsedStock , GetCycleTime()).GetSpeciesComposition(92).GetIsotopicQuantity();
				for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
					if ((*it).first.A() >= 232 && (*it).first.A() <= 236)
					{
						Sum_AlphaI_nUI0 += FuelType->GetFuelParameter()[(*it).first.A() -231]*(*it).second;
					}
			}
			
			double StockFactionToUse = 0;
			
			double NT = HMmass*1e6 * AVOGADRO / (cZAIMass.fZAIMass.find( ZAI(90,232,0) )->second);
			
			double N1 = (BU - FuelType->GetFuelParameter()[6]) * NT;
			double N2 = -Sum_AlphaI_nUI0;
			double N3 = -FuelType->GetFuelParameter()[0] * AVOGADRO / (cZAIMass.fZAIMass.find( ZAI(90,232,0) )->second) * (HMmass*1e6 - MU_0*1e6);
			
			double D1 = Sum_AlphaI_nUI;
			double D2 = -FuelType->GetFuelParameter()[0] * MU_1*1e6 * AVOGADRO / (cZAIMass.fZAIMass.find( ZAI(90,232,0) )->second) ;
			
			StockFactionToUse = (N1 + N2 + N3) / (D1 + D2);
			
			if(StockFactionToUse < 0)
			{
				stock.GetActinidesComposition().Print();

				cout << "!!Bad Trouble!! !!!FabricationPlant!!! Oups Bug in calculating stock fraction to use "<< endl;
				GetLog()->fLog << "!!Bad Trouble!! !!!FabricationPlant!!! Oups Bug in calculating stock fraction to use" << endl;
				RecycleStock(0.);
				FuelBuild = false;
				
			}
			else if( StockFactionToUse > 1 )
			{
				//cout<<"TEST TEST TEST TEST TEST Not enough Pu"<<endl;

				FullUsedStock += stock;
				RecycleStock(1);
				FuelBuild = false;
			}
			else
			{
				RecycleStock(StockFactionToUse);
				
				IsotopicVector IVBeginCycle;
				FuelBuild = true;
				
				ZAI Th = ZAI(90,232,0);
				//cout<<"TEST TEST TEST TEST TEST (HMmass - (MU_0+StockFactionToUse*MU_1 ) : "<<(HMmass - (MU_0+StockFactionToUse*MU_1 ))<<endl;
				//cout<<"TEST TEST TEST TEST TEST HMmass : "<<HMmass <<endl;
				//cout<<"TEST TEST TEST TEST TEST MU_0 : "<<MU_0<<endl;
				//cout<<"TEST TEST TEST TEST TEST StockFactionToUse*MU_1 : "<<StockFactionToUse<<"*"<<MU_1<<endl;
				double Th_Quantity = (HMmass - (MU_0+StockFactionToUse*MU_1 ))/(cZAIMass.fZAIMass.find( ZAI(90,232,0) )->second)*AVOGADRO/1e-6;
				
				GetParc()->AddGodIncome( Th, Th_Quantity);
				
				for(int i = (int)fFractionToTake.size()-1; i >= 0; i--)
				{
					IVBeginCycle += fStorage->GetIVArray()[fFractionToTake[i].first].GetSpeciesComposition(92)*( fFractionToTake[i].second );					
					fReUsable->AddIV(fStorage->GetIVArray()[fFractionToTake[i].first]*(fFractionToTake[i].second)
							      - fStorage->GetIVArray()[fFractionToTake[i].first].GetSpeciesComposition(92)*(fFractionToTake[i].second));
					
					fStorage->TakeFractionFromStock(fFractionToTake[i].first,fFractionToTake[i].second);			
					
				}
				fFractionToTake.clear();
				
				IVBeginCycle += Th_Quantity*Th;
				EvolutionData evolutiondb = BuildEvolutiveDB(ReactorId, IVBeginCycle);
				
				{
					pair<map<int, EvolutionData>::iterator, bool> IResult;
					IResult = fReactorFuturDB.insert( pair<int, EvolutionData>(ReactorId,evolutiondb) );
					if(!IResult.second)
						IResult.first->second = evolutiondb;
				}
				{
					pair<map<int, IsotopicVector>::iterator, bool> IResult;
					IResult = fReactorFuturIV.insert( pair<int, IsotopicVector>(ReactorId,IVBeginCycle) );
					if(!IResult.second)
						IResult.first->second = IVBeginCycle;

					AddCumulativeIVIn(IVBeginCycle);
					fInsideIV += IVBeginCycle;


				}
			}
		}
	}

}
