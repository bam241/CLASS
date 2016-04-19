#include "FabricationPlant.hxx"

#include "Storage.hxx"
#include "Reactor.hxx"
#include "EvolutionData.hxx"
#include "PhysicsModels.hxx"
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
#include <ctime>        
#include <cstdlib>  

ClassImp(FabricationPlant)

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

//________________________________________________________________________
FabricationPlant::FabricationPlant():CLASSFacility(16)
{
	SetName("F_FabricationPLant.");
	
	fReUsable = 0;
	fIsReusable = false;
}

//________________________________________________________________________
FabricationPlant::FabricationPlant(CLASSLogger* log, double fabricationtime):CLASSFacility(log, fabricationtime, 16)
{
DBGL
	SetName("F_FabricationPLant.");

	fStorageManagement = kpLiFo;
	fIsSeparationManagement = true;
	fSubstitutionFuel = false;
	fSubstitutionFissile = false;
	fIsReplaceFissileStock = false;

	fReUsable = 0;
	fIsReusable = false;

	INFO	 <<  " A FabricationPlant has been define :" << endl;
	INFO	 <<  "\t Chronological Stock Priority has been set! " <<  endl;
	INFO	 <<  "\t Fabrication time set to \t " << (double)(GetCycleTime()/cYear) << " year" << endl << endl;
DBGL
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
		IResult = fSeparationLostFraction.GetIsotopicQuantity().insert( pair<ZAI ,double>(zai, 1 - factor));
		if(!IResult.second)
			IResult.first->second = 1 - factor;
	}
	
}

//________________________________________________________________________
//_______________________________ Evolution ______________________________
//________________________________________________________________________

//________________________________________________________________________
void FabricationPlant::Evolution(cSecond t)
{
	
		// Check if the FabricationPlant has been created ...
	if(t ==  fInternalTime && t != 0) return;
		// Make the evolution for the FabricationPlant ...
	FabricationPlantEvolution(t);
		//Update Inside IsotopicVector
	UpdateInsideIV();
		// ... And Finaly update the AbsoluteInternalTime
	fInternalTime = t;
	
}

//________________________________________________________________________
void FabricationPlant::FabricationPlantEvolution(cSecond t)
{
DBGL
	map<int ,cSecond >::iterator it;
	for( it = fReactorNextStep.begin(); it !=  fReactorNextStep.end(); it++ )
	{
		double R_CreactionTime = GetParc()->GetReactor()[ (*it).first ]->GetCreationTime();
		double R_LifeTime = GetParc()->GetReactor()[ (*it).first ]->GetLifeTime();

		int ReactorId = (*it).first;
		pair<CLASSFuel, double> R_Fuel = GetParc()->GetReactor()[ReactorId]->GetFuelPlan()->GetFuelAt( t + GetCycleTime() );
		double R_BU = R_Fuel.second;
		double R_Power = GetParc()->GetReactor()[ReactorId]->GetPower();
		double R_HMMass = GetParc()->GetReactor()[ReactorId]->GetHeavyMetalMass();
		cSecond R_CycleTime = (cSecond) (R_BU*1e9 / (R_Power) * R_HMMass * 3600 * 24);


		if( R_CycleTime < GetCycleTime())
		{
			ERROR << "Reactor Cycle Time is shorter than Fabrication Time of the fuel, we cannot deal it upto now!!!" <<  endl;
			exit(1);
		}

		if( t + GetCycleTime() >= R_CreactionTime
		   && t + GetCycleTime() < R_CreactionTime + R_LifeTime)
		{
			if( (*it).second ==  t )
			{
#pragma omp critical(FuelBuild)
				{
					if( R_Fuel.first.GetPhysicsModels() )
					{
						BuildFuelForReactor( (*it).first, t );
					}
					(*it).second += R_CycleTime;
				}

			}
			else if ( (*it).second - R_CycleTime + GetCycleTime() >= t && (*it).second - R_CycleTime  < t )
			{
				map<int ,IsotopicVector >::iterator it2 = fReactorFuturIV.find( (*it).first );
				if (it2 != fReactorFuturIV.end())
					(*it2).second = GetDecay((*it2).second, t - fInternalTime );		
			}
		}
	}
	
	
DBGL
}

//________________________________________________________________________
void FabricationPlant::UpdateInsideIV()
{
	DBGL
	fInsideIV = IsotopicVector();

	map< int,IsotopicVector >::iterator it;
	for( it = fReactorFuturIV.begin(); it != fReactorFuturIV.end(); it++ )
		fInsideIV += (*it).second;

	DBGL
}

//________________________________________________________________________
void FabricationPlant::BuildFuelForReactor(int ReactorId, cSecond t)
{
	DBGL
	if(fFissileStorage.size() ==  0)
	{
		ERROR << " One need at least one Fissile storage to build fuel " << endl;
		ERROR << " Use AddFissileStorage to add a stock to provide fissil material... " << endl;
		ERROR << " One need at least one Fissile storage to build fuel " << endl;
		exit(1);
	}



	double R_HM_Mass = GetParc()->GetReactor()[ ReactorId ]->GetHeavyMetalMass();
	double R_CycleTime = GetParc()->GetReactor()[ ReactorId ]->GetCycleTime();
	double R_Power	 = GetParc()->GetReactor()[ ReactorId ]->GetPower();

	pair<CLASSFuel, double > FuelBU = GetParc()->GetReactor()[ReactorId]->GetFuelPlan()->GetFuelAt(t+GetCycleTime()) ;
	PhysicsModels FuelType = *FuelBU.first.GetPhysicsModels();
	double R_BU	      = FuelBU.second;
	
	fFissileList = FuelType.GetEquivalenceModel()->GetFissileList();
	BuildFissileArray();

	// If there is not enough Fissile its possible to take from an infinite fissile stock
	if( !fIsReplaceFissileStock && fSubstitutionFissile )//if it is defined and wanted
	{	IsotopicVector CooledSeparatedIV = GetDecay( fSubstitutionFissileIV , GetCycleTime());
		fFissileArray.push_back( CooledSeparatedIV/ CooledSeparatedIV.GetTotalMass() * R_HM_Mass );
	}


	fFertileList = FuelType.GetEquivalenceModel()->GetFertileList();


	if(fFertileStorage.size() != 0)			// If the fertile need to be taken in stock
		BuildFertileArray();
	else						// if their is not stock and the fertile come from outside of the park
	{
		fFertileArray.push_back( fFertileList / fFertileList.GetTotalMass() * R_HM_Mass );
		DBGV("Fertile Array size : " << fFertileArray.size())
	}

	
	vector<double> LambdaArray = FuelType.GetEquivalenceModel()->BuildFuel(R_BU, R_HM_Mass, fFissileArray, fFertileArray);

	double  LambdaSum = 0;
	for(int i = 0; i < (int)fFissileArray.size();i++)
		LambdaSum += LambdaArray[i];

	if(LambdaArray[0] != -1 && LambdaSum > 0 )
	{		
		DBGV("Building process suceed: ")

		IsotopicVector IV = BuildFuelFromEqModel(LambdaArray);
		EvolutionData EvolDB = FuelType.GenerateEvolutionData( GetDecay(IV,fCycleTime), R_CycleTime, R_Power);

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
		DBGL
		return;
	}
	else
	{
		DBGV("Building process failed: ")
		if (!fSubstitutionFuel)
		{
			DBGV("Reactor not loaded ")
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
			DBGV("Using substitute : ")
			if(fSubstitutionFissile)//if the build fail possibility to take a fissile material from an infinite fissile stock (if infinite stock defined)
			{
				DBGV("->From an infinite stock")
				//make the user specified fissil composition to decay the fabrication time
				IsotopicVector CooledSeparatedIV = GetDecay(fSubstitutionFissileIV, GetCycleTime());
				//Building the fuel :
				double MolarFissileContent = FuelType.GetEquivalenceModel()->GetFissileMolarFraction(CooledSeparatedIV,fFertileList, R_BU);
				IsotopicVector BuiltFuel = MolarFissileContent*fSubstitutionFissileIV/fSubstitutionFissileIV.GetSumOfAll() + (1-MolarFissileContent)*fFertileList/fFertileList.GetSumOfAll();
				IsotopicVector IV = BuiltFuel/ BuiltFuel.GetTotalMass() * R_HM_Mass;

				//Generating the EvolutionData
				EvolutionData EvolDB = FuelType.GenerateEvolutionData( GetDecay(IV,fCycleTime), R_CycleTime, R_Power);
		
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

				GetParc()->AddOutIncome(IV);
				fInsideIV += IV;
				AddCumulativeIVIn(IV);

				DBGL
			}
			else
			{
				DBGV("->From a fixed data base")
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
				GetParc()->AddOutIncome( IV );
				fInsideIV += IV;
				AddCumulativeIVIn(IV);
			}

		}
		DBGL
		ResetArrays();

		return;
	}
DBGL
}

//________________________________________________________________________
void FabricationPlant::SortArray()
{
    map < string , vector <IsotopicVector> >::iterator it_s_vIV;
    
    for( it_s_vIV = fStreamArray.begin();  it_s_vIV != fStreamArray.end(); it_s_vIV++)
    {
        
        vector<IsotopicVector> IVArray;
        vector<cSecond>	 TimeArray;
        vector< pair<int,int> >	 AdressArray;
        
        IVArray		= fStreamArray[(*it_s_vIV).first] ;
        TimeArray	= fStreamArrayTime[(*it_s_vIV).first] ;
        AdressArray	= fStreamArrayAdress[(*it_s_vIV).first] ;
        
        if(fFiFo)
        {
            for(int j = 0; j < (int)TimeArray.size(); j++)
            {
                for (int k = j+1; k < (int)TimeArray.size(); k++)
                {
                    cSecond time_tmp 		= TimeArray[j];
                    pair<int,int> Adress_tmp 	= AdressArray[j];
                    IsotopicVector IV_tmp 		= IVArray[j];
                    
                    if(time_tmp > TimeArray[k])
                    {
                        TimeArray[j] 	= TimeArray[k];
                        TimeArray[k] 	= time_tmp;
                        
                        AdressArray[j] 	= AdressArray[k];
                        AdressArray[k] 	= Adress_tmp;
                        
                        IVArray[j] 	= IVArray[k];
                        IVArray[k] 	= IV_tmp;
                    }
                }
            }
        }
        else
        {
            for(int j = 0; j < (int)TimeArray.size(); j++)
            {
                
                for (int k = j+1; k < (int)TimeArray.size(); k++)
                {
                    cSecond time_tmp 		= TimeArray[j];
                    pair<int,int> Adress_tmp 	= AdressArray[j];
                    IsotopicVector IV_tmp 		= IVArray[j];
                    
                    if(time_tmp < TimeArray[k])
                    {
                        TimeArray[j] 	= TimeArray[k];
                        TimeArray[k] 	= time_tmp;
                        
                        AdressArray[j] 	= AdressArray[k];
                        AdressArray[k] 	= Adress_tmp;
                        
                        IVArray[j] 	= IVArray[k];
                        IVArray[k] 	= IV_tmp;
                    }	
                }
            }
        }
        fStreamArray[(*it_s_vIV).first]		= IVArray;
        fStreamArrayTime[(*it_s_vIV).first]	= TimeArray;
        fStreamArrayAdress[(*it_s_vIV).first]	= AdressArray;
    }
}

//________________________________________________________________________
void FabricationPlant::SetSubstitutionFuel(EvolutionData fuel)
{
	
	fSubstitutionFuel = true;

	double M0 = cZAIMass.GetMass( fuel.GetIsotopicVectorAt(0.).GetActinidesComposition() );
	fSubstitutionEvolutionData = fuel / M0;

}

//________________________________________________________________________
void FabricationPlant::SortFiFo(vector<IsotopicVector>	&IVArray, vector<cSecond> &TimeArray, vector< pair<int,int> > &AdressArray)
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

//________________________________________________________________________
void FabricationPlant::SortLiFo(vector<IsotopicVector>	&IVArray, vector<cSecond> &TimeArray, vector< pair<int,int> > &AdressArray)
{
    for(int j = 0; j < (int)TimeArray.size(); j++)
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


//________________________________________________________________________
void FabricationPlant::SortMix(vector<IsotopicVector>	&IVArray, vector<cSecond> &TimeArray, vector< pair<int,int> > &AdressArray)
{
    
    //Sort by anti-chronoligical order (youger first)
    SortLiFo(IVArray, TimeArray, AdressArray);
    /*******Store it ******/
    vector<IsotopicVector>		Saved_IVArray		= IVArray	 ;
    vector<cSecond> 			Saved_TimeArray		= TimeArray	 ;
    vector< pair<int,int> > 	Saved_AdressArray	= AdressArray;
    
    int IVsize = (int)IVArray.size();
    /*******Then reset the vectors ******/
    IVArray.clear();
    TimeArray.clear();
    AdressArray.clear();
    
    
    int HalfSize = floor( (double)IVsize/2. );
    
    int old = IVsize;
    
    bool isYoung=true;//change to false to begin with an old isotopicvector
    
    int RemainingIV = -1;
    
    
    for(int young = 0 ; young < HalfSize ; young++)
    {
        
        if(isYoung)
        {
            IVArray.push_back(Saved_IVArray[young]);
            TimeArray.push_back(Saved_TimeArray[young]);
            AdressArray.push_back(Saved_AdressArray[young]);
            
            isYoung=!isYoung;
            RemainingIV = young+1; // +1 ? -> The next young will be the +1
        }
        if(!isYoung)
        {
            old--;
            
            IVArray.push_back(Saved_IVArray[old]);
            TimeArray.push_back(Saved_TimeArray[old]);
            AdressArray.push_back(Saved_AdressArray[old]);
            
            isYoung=!isYoung;
            RemainingIV = old-1; // -1 ? -> The next old will be the -1
            
        }
        
    }
    
    if(  (double)IVsize/2. - (double)HalfSize  > 0.0 ) //if odd number of isotopic vector : one isotopic vector is still missing add it @ the end
    {
        
        IVArray.push_back(Saved_IVArray[RemainingIV]);
        TimeArray.push_back(Saved_TimeArray[RemainingIV]);
        AdressArray.push_back(Saved_AdressArray[RemainingIV]);
        
    }
    
    
}

//________________________________________________________________________
void FabricationPlant::SortRandom(vector<IsotopicVector>	&IVArray, vector<cSecond> &TimeArray, vector< pair<int,int> > &AdressArray)
{
    int SizeOfIVArray = IVArray.size();
    
    /*********Create a Random list of vector position**********/
    srand ( unsigned ( std::time(0) ) );
    vector<int> RandomPosition;
    for (int i=0 ; i < (int) SizeOfIVArray ; ++i) 
        RandomPosition.push_back(i); 
    
    random_shuffle(RandomPosition.begin(), RandomPosition.end());
    
    /*******Store old vectors ******/
    vector<IsotopicVector>		Saved_IVArray		= IVArray	 ;
    vector<cSecond> 			Saved_TimeArray		= TimeArray	 ; 
    vector< pair<int,int> > 	Saved_AdressArray	= AdressArray;
    
    /*******Asign values ******/
    
    for (int i=0 ; i < (int) SizeOfIVArray ; ++i) 
    {
        IVArray[i]		=	Saved_IVArray[RandomPosition[i]];
        TimeArray[i]	=	Saved_TimeArray[RandomPosition[i]];
        AdressArray[i]	=	Saved_AdressArray[RandomPosition[i]];
    }
    
}	



//________________________________________________________________________
//_____________________________ Reactor & DB _____________________________
//________________________________________________________________________

//________________________________________________________________________
void FabricationPlant::TakeReactorFuel(int Id)
{
DBGL
	IsotopicVector IV;
	map<int ,IsotopicVector >::iterator it2 = fReactorFuturIV.find( Id );

	AddCumulativeIVOut(it2->second);

	if (it2 != fReactorFuturIV.end())
		(*it2).second = IV;

	map< int,EvolutionData >::iterator it = fReactorFuturDB.find(Id);
	(*it).second.DeleteEvolutionDataCopy();

	UpdateInsideIV();
DBGL
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

	int StockCorrection = 0;
	if( !fIsReplaceFissileStock && fSubstitutionFissile)
	{	
		StockCorrection = 1;
		BuildedFuel += fFissileArray.back()*LambdaArray[fFissileArray.size()-1];
	}

	for(int i = 0; i < (int)fFissileArray.size() - StockCorrection ; i++)
	{
		if(LambdaArray[i] != 0)
		{
			int Stor_N = fFissileArrayAdress[i].first;
			int IV_N = fFissileArrayAdress[i].second;

			pair<IsotopicVector, IsotopicVector> Separated_Lost;
			Separated_Lost = Separation( fFissileStorage[Stor_N]->GetIVArray()[IV_N]*LambdaArray[i], fFissileList );
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
				Separated_Lost = Separation( fFertileStorage[Stor_N]->GetIVArray()[IV_N]*LambdaArray[i], fFertileList);
				BuildedFuel += Separated_Lost.first;
				Lost += Separated_Lost.second;
			}
		}
	}
	else
		BuildedFuel += fFertileArray[0]*LambdaArray.back();

	if(fIsReusable)
		fReUsable->AddIV(Lost);
	else
		GetParc()->AddWaste(Lost);

	DumpStock(LambdaArray);

DBGL
	return BuildedFuel;
}

//________________________________________________________________________
void FabricationPlant::DumpStock(vector<double> LambdaArray)
{
DBGL
	int StockCorrection = 0;
	if( !fIsReplaceFissileStock && fSubstitutionFissile)
	{	StockCorrection = 1;
		GetParc()->AddOutIncome( fFissileArray.back()*LambdaArray[fFissileArray.size()-1] );
	}
	for(int i = 0; i < (int)fFissileArray.size() - StockCorrection; i++)
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
		GetParc()->AddOutIncome( fFertileArray[0]*LambdaArray.back() );

	ResetArrays();


DBGL
}

//________________________________________________________________________
void FabricationPlant::ResetArrays()
{
		//Clear the Building Array (Fissile and Fertile)
	fFissileArray.clear();
	fFissileArrayTime.clear();
	fFissileArrayAdress.clear();
	fFertileArray.clear();
	fFertileArrayTime.clear();
	fFertileArrayAdress.clear();

	fFertileList = fFissileList = IsotopicVector();
}

//________________________________________________________________________
pair<IsotopicVector, IsotopicVector> FabricationPlant::Separation(IsotopicVector isotopicvector, IsotopicVector ExtractedList)
{
DBGL
	IsotopicVector SeparatedPart;
	IsotopicVector LostPart;

	if(fIsSeparationManagement)
	{
		//[0] = re-use ; [1] = waste
		IsotopicVector LostInReprocessing  = isotopicvector.GetThisComposition(ExtractedList) * fSeparationLostFraction;
		SeparatedPart  = isotopicvector.GetThisComposition(ExtractedList) - LostInReprocessing;
		LostPart = isotopicvector - SeparatedPart;
	}
	else
	{
		//[0] = re-use ; [1] = waste
		//IsotopicVector LostInReprocessing  = isotopicvector.GetThisComposition(ExtractedList) * fSeparationLostFraction;
		SeparatedPart  = isotopicvector;
		LostPart = isotopicvector - SeparatedPart;
	}
DBGL
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

