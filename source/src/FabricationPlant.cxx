#include "FabricationPlant.hxx"

#include "Storage.hxx"
#include "Reactor.hxx"
#include "EvolutionData.hxx"
#include "DecayDataBank.hxx"
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
#include <stdarg.h>

#include <ctime>
#include <cstdlib>

ClassImp(FabricationPlant)

//________________________________________________________________________
//________________________________________________________________________
//
//		FabricationPlant
//________________________________________________________________________
//________________________________________________________________________


FabricationPlant::FabricationPlant():CLASSFacility(16)
{
    SetName("F_FabricationPLant.");
    
    fReUsable = 0;
    fIsReusable = false;

    IsotopicVector MaxEfficiency;
    SetSeparationEfficiency(MaxEfficiency,  -1e300);
}

//________________________________________________________________________
FabricationPlant::FabricationPlant(CLASSLogger* log, double fabricationtime):CLASSFacility(log, fabricationtime, 16)
{
    DBGL
    SetName("F_FabricationPLant.");
    
    fStorageManagement = kpLiFo;
    fIsSeparationManagement = true;
    fSubstitutionFuel = false;
    
    fReUsable = 0;
    fIsReusable = false;

    IsotopicVector MaxEfficiency;
    SetSeparationEfficiency(MaxEfficiency,  -1e300);

    INFO	<< " A FabricationPlant has been define :" << endl;
    INFO	<< "\t Chronological Stock Priority has been set! "<< endl;
    INFO	<< "\t Fabrication time set to \t " << (double)(GetCycleTime()/cYear) << " year" << endl << endl;
    DBGL
}

//________________________________________________________________________
FabricationPlant::~FabricationPlant()
{
    
    
}

//________________________________________________________________________
void    FabricationPlant::SetSeparationEfficiency(IsotopicVector IV,  cSecond TimeOfSep)
{
    DBGL
    IsotopicVector SeparationEfficiency;
    IsotopicVector UnityVector;

    map<ZAI, double>::iterator it;
    map<ZAI, double> IVb = IV.GetIsotopicQuantity();

    for(it = IVb.begin() ; it != IVb.end() ; it++ )
    {   
        if(it->second > 1 )
        {
            ERROR << " Efficiency must be below one";
            exit(1);
        }
        else if (it->second < 0)
        {
            ERROR << " Efficiency must be positive ";
            exit(1);
        }

        UnityVector.Add(it->first,1);
    } 



    SeparationEfficiency = UnityVector - IV ;

    fSeparationStrategy.insert( pair <cSecond, IsotopicVector>(TimeOfSep, SeparationEfficiency )  );

}

//________________________________________________________________________
IsotopicVector  FabricationPlant::GetSeparationEfficiencyAt(cSecond time)
{

    DBGV("Getting Separation Efficiency at time : "<< (time/3600/24./365.25) <<" y")
    map<cSecond , IsotopicVector>::iterator itlow;

    itlow = fSeparationStrategy.lower_bound(time);
    DBGV("Wanted time : " << itlow->first /3600/24./365.25 <<" y" )

    if(itlow->first == time){

        IsotopicVector IV =  itlow->second;
        DBGV(IV.sPrint())
        return IV;
    }
    else if(itlow != fSeparationStrategy.begin()){
        
        IsotopicVector IV =   (--itlow)->second;
         DBGV(IV.sPrint());
        return IV ;
    }
    else {
        
        IsotopicVector IV =  fSeparationStrategy.begin()->second;
         DBGV(IV.sPrint())

        return IV;

    }

}




//________________________________________________________________________
//_______________________________ Evolution ______________________________
//________________________________________________________________________

//________________________________________________________________________
void FabricationPlant::Evolution(cSecond t)
{
    
    // Check if the FabricationPlant has been created ...
    if(t == fInternalTime && t != 0) return;
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
    for( it = fReactorNextStep.begin(); it!= fReactorNextStep.end(); it++ )
    {
        double R_CreactionTime = GetParc()->GetReactor()[ (*it).first ]->GetCreationTime();
        double R_LifeTime 	   = GetParc()->GetReactor()[ (*it).first ]->GetLifeTime();
        
        int ReactorId = (*it).first;
        ScheduleEntry* R_Entry = GetParc()->GetReactor()[ReactorId]->GetScheduler()->GetEntryAt(t + GetCycleTime());
        double R_BU 		   = R_Entry->GetBurnUp();
        double R_Power 	       = R_Entry->GetPower();
        double R_HMMass 	   = R_Entry->GetHeavyMetalMass();
        cSecond R_CycleTime    = (cSecond) (R_BU*1e9 / (R_Power) * R_HMMass * 3600 * 24);

        DBGL
        if( R_CycleTime < GetCycleTime())
        {
            ERROR << "Reactor Cycle Time is shorter than Fabrication Time of the fuel, we cannot deal it upto now!!!"<< endl;
            exit(1);
        }
        
        if( t + GetCycleTime() >= R_CreactionTime
           && t + GetCycleTime() < R_CreactionTime + R_LifeTime)
        {
            if( (*it).second == t )
            {
#pragma omp critical(FuelBuild)
                {
                    if( R_Entry->GetReactorModel()->GetPhysicsModels() )
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

    ScheduleEntry* R_Entry  = GetParc()->GetReactor()[ ReactorId ]->GetScheduler()->GetEntryAt(t + GetCycleTime());
    double R_BU             = R_Entry->GetBurnUp();
    double R_Power          = R_Entry->GetPower();
    double R_HM_Mass        = R_Entry->GetHeavyMetalMass();
    PhysicsModels* FuelType = R_Entry->GetReactorModel()->GetPhysicsModels();
    cSecond R_CycleTime     = (cSecond) (R_BU*1e9 / (R_Power) * R_HM_Mass * 3600 * 24);

    if( !fSeparationStrategy.empty() )
      fSeparationLostFraction = GetSeparationEfficiencyAt(t);

    map < string , vector <IsotopicVector> >::iterator it_s_vIV;
    map < string , vector <double> >::iterator it_s_vD;
    map < string , bool >::iterator it_s_B;
    
    fStreamList = FuelType->GetEquivalenceModel()->GetAllStreamList();
    
    BuildArray(ReactorId, t +  GetCycleTime()); // Checker chez les stocks si les StreamList sont présentes
    // Grosse map qui contient tous les IV
    // séparation + refroidissement virtuel
    // Construit les stocks de matière infini (=taille du réacteur)
    
    // string="MA, .." LambdaArray = tableau sur les IV
    map < string , vector <double> > LambdaArray =  FuelType->GetEquivalenceModel()->BuildFuel(R_BU, R_HM_Mass, fStreamArray, fStreamListFPMassFractionMin, fStreamListFPMassFractionMax, fStreamListFPPriority, fStreamListFPIsBuffer);
    
     fFuelCanBeBuilt 	  = true;
    double  LambdaSum  = 0;
    
     map < string, IsotopicVector>::iterator it_s_IV;

     //initialize error map
    for( it_s_IV = fStreamList.begin();  it_s_IV != fStreamList.end(); it_s_IV++)
        fErrorOnLambda[(*it_s_IV).first] = false;
    

    for( it_s_IV = fStreamList.begin();  it_s_IV != fStreamList.end(); it_s_IV++)
    {
        string StreamName =  it_s_IV->first ;

        if(find(LambdaArray[StreamName].begin(), LambdaArray[StreamName].end(), -1 )!= LambdaArray[StreamName].end()) //There is IV but with  enought material
            fErrorOnLambda[StreamName] = true;

        for(int i = 0; i < (int)LambdaArray[StreamName].size();i++)
            LambdaSum += LambdaArray[StreamName][i];
        
    }
    
    for( it_s_B = fErrorOnLambda.begin();  it_s_B != fErrorOnLambda.end(); it_s_B++)
    {
        if(fErrorOnLambda[(*it_s_B).first]){fFuelCanBeBuilt = false;}
    }
    
    if(fFuelCanBeBuilt && LambdaSum > 0 )
    {
        DBGV("Building process from initial stocks has succeeded : ")
        IsotopicVector IV 		= BuildFuelFromEqModel(LambdaArray);
        IsotopicVector LoadedIV 	           = GetDecay(IV,fCycleTime);
        
        double Pu8   = LoadedIV.GetZAIIsotopicQuantity(94,238,0);
        double Pu9   = LoadedIV.GetZAIIsotopicQuantity(94,239,0);
        double Pu10 = LoadedIV.GetZAIIsotopicQuantity(94,240,0);
        double Pu11 = LoadedIV.GetZAIIsotopicQuantity(94,241,0);
        double Pu12 = LoadedIV.GetZAIIsotopicQuantity(94,242,0);
        double Am1 = LoadedIV.GetZAIIsotopicQuantity(95,241,0);
        double U5    = LoadedIV.GetZAIIsotopicQuantity(92,235,0);
        double U8    = LoadedIV.GetZAIIsotopicQuantity(92,238,0);
        double Ntot  = Pu8 + Pu9 + Pu10 + Pu11 + Pu12 +Am1 + U5 +U8;
        
        Pu8         = Pu8/Ntot;
        Pu9         = Pu9/Ntot;
        Pu10       = Pu10/Ntot;
        Pu11       = Pu11/Ntot;
        Pu12       = Pu12/Ntot;
        Am1        = Am1/Ntot;
        U5          = U5/Ntot;
        U8          = U8/Ntot;

        double eU5      = U5/(1-(Pu8 + Pu9 + Pu10 + Pu11 + Pu12 +Am1));
        double wPu      = Pu8 + Pu9 + Pu10 + Pu11 + Pu12 +Am1;

        Pu8         = Pu8/wPu;
        Pu9         = Pu9/wPu;
        Pu10       = Pu10/wPu;
        Pu11       = Pu11/wPu;
        Pu12       = Pu12/wPu;
        Am1        = Am1/wPu;
        
        double m_U5         = LoadedIV.GetZAIIsotopicQuantity(92,235,0)*235e-06/6.023e23/LoadedIV.GetTotalMass();
        double m_U8         = LoadedIV.GetZAIIsotopicQuantity(92,238,0)*238e-06/6.023e23/LoadedIV.GetTotalMass();
        double m_Pu8        = LoadedIV.GetZAIIsotopicQuantity(94,238,0)*238e-06/6.023e23/LoadedIV.GetTotalMass();
        double m_Pu9        = LoadedIV.GetZAIIsotopicQuantity(94,239,0)*239e-06/6.023e23/LoadedIV.GetTotalMass();
        double m_Pu10       = LoadedIV.GetZAIIsotopicQuantity(94,240,0)*240e-06/6.023e23/LoadedIV.GetTotalMass();
        double m_Pu11       = LoadedIV.GetZAIIsotopicQuantity(94,241,0)*241e-06/6.023e23/LoadedIV.GetTotalMass();
        double m_Pu12       = LoadedIV.GetZAIIsotopicQuantity(94,242,0)*242e-06/6.023e23/LoadedIV.GetTotalMass();
        double m_Am1        = LoadedIV.GetZAIIsotopicQuantity(95,241,0)*241e-06/6.023e23/LoadedIV.GetTotalMass();

        double M_Pu8    = LoadedIV.GetZAIIsotopicQuantity(94,238,0)*238e-06/6.023e23;
        double M_Pu9    = LoadedIV.GetZAIIsotopicQuantity(94,239,0)*239e-06/6.023e23;
        double M_Pu10  = LoadedIV.GetZAIIsotopicQuantity(94,240,0)*240e-06/6.023e23;
        double M_Pu11  = LoadedIV.GetZAIIsotopicQuantity(94,241,0)*241e-06/6.023e23;
        double M_Pu12  = LoadedIV.GetZAIIsotopicQuantity(94,242,0)*242e-06/6.023e23;
        double M_Am1   = LoadedIV.GetZAIIsotopicQuantity(95,241,0)*241e-06/6.023e23;

        double m_Pu = M_Pu8 + M_Pu9 + M_Pu10 + M_Pu11 + M_Pu12 + M_Am1;
        double w_mPu = (m_Pu8 + m_Pu9 + m_Pu10 + m_Pu11 + m_Pu12 + m_Am1); 

        //cout<<"Frac "<<IV.GetTotalMass()<<" "<<t/3600./24./365.25<<" "<<w_mPu<<"   "<<wPu<<"   "<<U5<<"   "<<"    "<<U8<<"    "<<Pu8<<"    "<<Pu9<<"    "<<Pu10<<"    "<<Pu11<<"    "<<Pu12<<"    "<<Am1<<endl; 
        //cout<<"Masse "<<IV.GetTotalMass()<<" "<<t/3600./24./365.25<<" "<<m_Pu<<" "<<w_mPu<<" "<<m_U5<<" "<<m_U8<<"   "<<m_Pu8<<"   "<<m_Pu9<<"   "<<m_Pu10<<"    "<<m_Pu11<<"    "<<m_Pu12<<"    "<<m_Am1<<endl; 
        //cout<<endl;
   
        EvolutionData EvolDB = FuelType->GenerateEvolutionData(GetDecay(IV,fCycleTime), R_CycleTime, R_Power);

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
        DBGV("Building process from initial stocks failed: ")
        //If the building fails : possibility to take a materials from an infinite stock composed of a substitution IV defined by user
        bool IsSubstitutionMaterials = false;
        
        for( it_s_B = fErrorOnLambda.begin();  it_s_B != fErrorOnLambda.end(); it_s_B++)
        {
            if(fErrorOnLambda[(*it_s_B).first] && fSubstitutionMaterialFromIV[(*it_s_B).first])
                IsSubstitutionMaterials = true;
            else if(fErrorOnLambda[(*it_s_B).first] && !fSubstitutionMaterialFromIV[(*it_s_B).first])
                IsSubstitutionMaterials = false;
        }
        
        if(IsSubstitutionMaterials)
        {
            DBGV("Using substitute : -> From infinite substitutions IV ")
            
            //Make the user specified composition to decay the fabrication time
            
            map < string , IsotopicVector> CooledSeparatedIV;
            
            for( it_s_B = fSubstitutionMaterialFromIV.begin();  it_s_B != fSubstitutionMaterialFromIV.end(); it_s_B++)
                CooledSeparatedIV[(*it_s_B).first] = GetDecay(fSubstitutionIV[(*it_s_B).first], GetCycleTime());
            
            for( it_s_vIV = fStreamArray.begin();  it_s_vIV != fStreamArray.end(); it_s_vIV++)
            {
                if(fErrorOnLambda[it_s_vIV->first])
                {
                    fStreamArray[it_s_vIV->first].clear();
                    fStreamArray[it_s_vIV->first].push_back(CooledSeparatedIV[(*it_s_vIV).first]);

                }
            }
            
            //Building the fuel :
            for( it_s_vD = LambdaArray.begin();  it_s_vD != LambdaArray.end(); it_s_vD++)
                LambdaArray[(*it_s_vD).first].clear();
            
            LambdaArray 			= FuelType->GetEquivalenceModel()->BuildFuel(R_BU, R_HM_Mass, fStreamArray, fStreamListFPMassFractionMin, fStreamListFPMassFractionMax, fStreamListFPPriority, fStreamListFPIsBuffer);
            IsotopicVector IV 		= BuildFuelFromEqModel(LambdaArray);
            
            //Generating the EvolutionData
            EvolutionData EvolDB = FuelType->GenerateEvolutionData(GetDecay(IV,fCycleTime), R_CycleTime, R_Power);
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
        else if(fSubstitutionFuel)
        {
            DBGV("Using substitute : -> From a fixed data base")
            IsotopicVector IV 		= fSubstitutionEvolutionData.GetIsotopicVectorAt(0);
            EvolutionData evolutiondb 	= fSubstitutionEvolutionData * R_HM_Mass / IV.GetTotalMass();
            
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
        else
        {
            DBGV("No Alternative Solution. Reactor not loaded. ")
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
        
        DBGL
        ResetArrays();
        
        return;
    }
    DBGL
}
//________________________________________________________________________
void FabricationPlant::BuildArray(int ReactorId, cSecond ReactorLoadingTime)
{
    DBGL
    ScheduleEntry* R_Entry  = GetParc()->GetReactor()[ ReactorId ]->GetScheduler()->GetEntryAt(ReactorLoadingTime);
    double R_HM_Mass         = R_Entry->GetHeavyMetalMass();

   
    vector <IsotopicVector>  StreamArray;
    vector <cSecond> 	   StreamArrayTime;
    vector < pair<int,int> >  StreamArrayAdress;
    
    map < string , IsotopicVector>::iterator it;
    for( it = fStreamList.begin();  it != fStreamList.end(); it++)
    {
        if(fInfiniteMaterialFromList[(*it).first])
        {
            IsotopicVector IV = (*it).second / (*it).second.GetTotalMass() * R_HM_Mass;
            StreamArray.push_back(IV);
            StreamArrayAdress.push_back(pair<int,int>(0,0));
            StreamArrayTime.push_back(0);
        }
        else
        {   
            for(int j = 0; j < (int)fStorage[(*it).first].size(); j++)
            {
                vector<IsotopicVector> IVArray = fStorage[(*it).first][j]->GetIVArray();
                for(int k = 0; k < (int)IVArray.size(); k++)
                {
                    IsotopicVector SeparatedIV = Separation(IVArray[k], (*it).second).first;
                    if(Norme(SeparatedIV) != 0)
                    {
                        IsotopicVector CooledSeparatedIV = GetDecay( SeparatedIV , GetCycleTime());
                        StreamArray.push_back( CooledSeparatedIV );
                        StreamArrayAdress.push_back( pair<int,int>(j,k) );
                        StreamArrayTime.push_back(fStorage[(*it).first][j]->GetIVArrayArrivalTime()[k]);
                    }
                }
            }
        }

        fStreamArray[(*it).first] 	         = StreamArray;			StreamArray.clear();
        fStreamArrayAdress[(*it).first]	= StreamArrayAdress;		StreamArrayAdress.clear();
        fStreamArrayTime[(*it).first]	= StreamArrayTime; 		StreamArrayTime.clear();
    }
    
    SortArray();
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
        
        if( (int)IVArray.size() > 1 ) // no need to sort if there is only one IV in stock !!
        {
        	
     	  switch(fStorageManagement)
     	  {
     	      case kpFiFo: SortFiFo(IVArray, TimeArray, AdressArray);
     	          break;
     	          
     	      case kpLiFo: SortLiFo(IVArray, TimeArray, AdressArray);
     	          break;
     	          
     	      case kpMix: SortMix(IVArray, TimeArray, AdressArray);
     	          break;
     	          
     	      case kpRand: SortRandom(IVArray, TimeArray, AdressArray);
     	          break;
     	          
     	      default:
     	          ERROR<<" Posibble Storage Management are"<<endl;
     	          ERROR<<" YourFabPlant->SetStorageManagement(key); //with key can be either"<<endl;
     	          ERROR<<"\tkFiFo : First In First Out (i.e the older storage first)"<<endl;
     	          ERROR<<"\tkLiFo : Last In First Out  (i.e the youger storage first)"<<endl;
     	          ERROR<<"\tkMix : IVs are sorted that way : "<<"\n"<<"First: The younger , Second: The older, Third: The second younger ,4th : the second older ...."<<endl;
     	          ERROR<<"\tkRand : IVs order in storage is random"<<endl;
     	          
     	          exit(1);
     	  }
        
        fStreamArray[(*it_s_vIV).first]		= IVArray;
        fStreamArrayTime[(*it_s_vIV).first]	= TimeArray;
        fStreamArrayAdress[(*it_s_vIV).first]	= AdressArray;

       } 
    }
    
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
    vector<IsotopicVector>	Saved_IVArray		= IVArray	 ;
    vector<cSecond> 		Saved_TimeArray	= TimeArray	 ;
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
    vector<IsotopicVector>	Saved_IVArray		= IVArray	 ;
    vector<cSecond> 		Saved_TimeArray	= TimeArray	 ;
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
//	Substitution Fuel
//________________________________________________________________________

void FabricationPlant::SetSubstitutionFuel(EvolutionData fuel)
{
    
    fSubstitutionFuel = true;
    
    double M0 = cZAIMass.GetMass( fuel.GetIsotopicVectorAt(0.).GetActinidesComposition() );
    fSubstitutionEvolutionData = fuel / M0;
    
}

//________________________________________________________________________
//_____________________________ Reactor & DB _____________________________
void FabricationPlant::TakeReactorFuel(int Id)
{
    DBGL
    IsotopicVector IV;
    map<int ,IsotopicVector >::iterator it2 = fReactorFuturIV.find( Id );
    
    AddCumulativeIVOut(it2->second);
    
    if (it2 != fReactorFuturIV.end())
        (*it2).second = IV;
    
    map< int,EvolutionData >::iterator it = fReactorFuturDB.find(Id);
    (*it).second = EvolutionData();
    
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

//________________________________________________________________________
IsotopicVector FabricationPlant::BuildFuelFromEqModel(map <string , vector<double> > LambdaArray)
{
    DBGL
    IsotopicVector BuildedFuel;
    IsotopicVector Lost;
    
    map < string , IsotopicVector>::iterator it;
    
    for( it = fStreamList.begin();  it != fStreamList.end(); it++)
    {
        for(int i = 0; i < (int)fStreamArray[(*it).first].size(); i++)
        {
            if(fInfiniteMaterialFromList[(*it).first])
                BuildedFuel += fStreamArray[(*it).first][i]*LambdaArray[(*it).first][i];
            
            else if (fSubstitutionMaterialFromIV[(*it).first] && fErrorOnLambda[(*it).first])
                BuildedFuel += fStreamArray[(*it).first][i]*LambdaArray[(*it).first][i];
            else
            {
                if(LambdaArray[(*it).first][i] != 0)
                {
                    int Stor_N 	= fStreamArrayAdress[(*it).first][i].first;
                    int IV_N 	= fStreamArrayAdress[(*it).first][i].second;
                    
                    pair<IsotopicVector, IsotopicVector> Separated_Lost;
                    Separated_Lost = Separation( fStorage[(*it).first][Stor_N]->GetIVArray()[IV_N]*LambdaArray[(*it).first][i], (*it).second );
                    BuildedFuel += Separated_Lost.first;
                    Lost += Separated_Lost.second;
                    
                }
            }
        }
    }
    
    if(fIsReusable)
        fReUsable->AddIV(Lost);
    else
        GetParc()->AddWaste(Lost);
    
    DumpStock(LambdaArray);
    
    DBGL
    return BuildedFuel;
}

//________________________________________________________________________
void FabricationPlant::DumpStock(map <string , vector<double> > LambdaArray)
{
    DBGL
    
    map < string , IsotopicVector>::iterator it;
    
    for( it = fStreamList.begin();  it != fStreamList.end(); it++)
    {
        for(int i = 0; i < (int)fStreamArray[(*it).first].size(); i++)
        {
            if(fInfiniteMaterialFromList[(*it).first])
            {DBGL
                GetParc()->AddOutIncome( fStreamArray[(*it).first][0]*LambdaArray[(*it).first][i] );
                DBGL
            }	
            
            else
            {		
                if(LambdaArray[(*it).first][i] != 0 && fFuelCanBeBuilt)
                {
                    DBGL
                    int Stor_N = fStreamArrayAdress[(*it).first][i].first;
                    DBGL
                    int IV_N = fStreamArrayAdress[(*it).first][i].second;
                    DBGL
                    fStorage[(*it).first][Stor_N]->TakeFractionFromStock( IV_N, LambdaArray[(*it).first][i] );
                }
            }
        }
    }

    ResetArrays();
        
    DBGL
}

//________________________________________________________________________
void FabricationPlant::ResetArrays()
{
    //Clear the Building Array 
    map < string , IsotopicVector>::iterator it;
    
    for( it = fStreamList.begin();  it != fStreamList.end(); it++)
    {		
        fStreamList[(*it).first]= IsotopicVector();
        fStreamArray[(*it).first].clear();
        fStreamArrayTime[(*it).first].clear();
        fStreamArrayAdress[(*it).first].clear();
        
    }
    
    
}

//________________________________________________________________________
pair<IsotopicVector, IsotopicVector> FabricationPlant::Separation(IsotopicVector isotopicvector, IsotopicVector ExtractedList)
{
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
        SeparatedPart  = isotopicvector;
        LostPart = isotopicvector - SeparatedPart;
    }
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

//________________________________________________________________________
void FabricationPlant::AddStorage(string keyword, Storage* Stock, double MassFractionMin, double MassFractionMax, int Priority)
{
    fStorage[keyword].push_back(Stock);

    if (MassFractionMin>MassFractionMax)
    {
            ERROR << " Mass fraction min of material : "<<keyword<<"has to be lower than mass fraction max." <<endl;
            exit(1);
    }
    if (MassFractionMin<0.0)
    {
            ERROR << " Mass fraction min of material : "<<keyword <<"has to be higher than zero." <<endl;
            exit(1);
    }
   if (MassFractionMax>1.0)
    {
            ERROR << " Mass fraction max of material : "<<keyword <<"has to be lower than one." <<endl;
            exit(1);
    }
    
    map < int , string >::iterator it_i_s;
    map < string , bool >::iterator it_s_B;
    for( it_i_s = fStreamListFPPriority.begin();  it_i_s != fStreamListFPPriority.end(); it_i_s++)
    { 
            if ((*it_i_s).first == Priority && (*it_i_s).second != keyword)
            {
                    ERROR << "2 materials have the same position in priority order : "<<(*it_i_s).second<<" and "<<keyword<< endl;
                    exit(1);

            }
    }

    for( it_s_B = fStreamListFPIsBuffer.begin();  it_s_B != fStreamListFPIsBuffer.end(); it_s_B++)
    { 
            if ((*it_s_B).first == keyword && fStreamListFPIsBuffer[(*it_s_B).first]==1)
            {
                    ERROR << "The material : "<<keyword<<" is already defined as buffer. "<< endl;
                    exit(1);
            }
    }
    
    fStreamListFPMassFractionMin[keyword] = MassFractionMin;
    fStreamListFPMassFractionMax[keyword] = MassFractionMax;


    fStreamListFPPriority[Priority]           = keyword;
    fStreamListFPIsBuffer[keyword]        = 0;

} 
//________________________________________________________________________
void FabricationPlant::AddInfiniteStorage(string keyword, double MassFractionMin, double MassFractionMax, int Priority)
{
    Storage* Stock;
    fStorage[keyword].push_back(Stock);

    if (MassFractionMin>MassFractionMax)
    {
            ERROR << " Mass fraction min of material : "<<keyword<<"has to be lower than mass fraction max." <<endl;
            exit(1);
    }
    if (MassFractionMin<0.0)
    {
            ERROR << " Mass fraction min of material : "<<keyword <<"has to be higher than zero." <<endl;
            exit(1);
    }
   if (MassFractionMax>1.0)
    {
            ERROR << " Mass fraction max of material : "<<keyword <<"has to be lower than one." <<endl;
            exit(1);
    }

    map < int , string >::iterator it_i_s;
    map < string , bool >::iterator it_s_B;
    for( it_i_s = fStreamListFPPriority.begin();  it_i_s != fStreamListFPPriority.end(); it_i_s++)
    { 
            if ((*it_i_s).first == Priority && (*it_i_s).second != keyword)
            {
                    ERROR << "2 materials have the same position in priority order : "<<(*it_i_s).second<<" and "<<keyword<< endl;
                    exit(1);

            }
    }
    for( it_s_B = fStreamListFPIsBuffer.begin();  it_s_B != fStreamListFPIsBuffer.end(); it_s_B++)
    { 
            if ((*it_s_B).first == keyword && fStreamListFPIsBuffer[(*it_s_B).first]==1)
            {
                    ERROR << "The material : "<<keyword<<" is already defined as buffer. "<< endl;
                    exit(1);
            }
    }
    fStreamListFPMassFractionMin[keyword] = MassFractionMin;
    fStreamListFPMassFractionMax[keyword] = MassFractionMax;

    fStreamListFPPriority[Priority]         = keyword;
    fStreamListFPIsBuffer[keyword]      = 0;

    fInfiniteMaterialFromList[keyword] = true;
} 
//________________________________________________________________________
void FabricationPlant::AddFuelBuffer(string keyword)
{
    Storage* Stock;
    fStorage[keyword].push_back(Stock);

     //Test if there is no buffer already defined//
    map < string , bool >::iterator it_s_B;
    for( it_s_B = fStreamListFPIsBuffer.begin();  it_s_B != fStreamListFPIsBuffer.end(); it_s_B++)
    { 
            if ((*it_s_B).second == 1)
            {
                    ERROR << (*it_s_B).first<<" is already defined as a buffer. "<< endl;
                    ERROR << " Fuel can't have more than one buffer in current algorithm. "<< endl;      
                    exit(1);
            }
    }
    map < int , string >::iterator it_i_s;
    for( it_i_s = fStreamListFPPriority.begin();  it_i_s != fStreamListFPPriority.end(); it_i_s++)
    { 
            if ((*it_i_s).second == keyword)
            {
                    ERROR << "The material "<<keyword<<" can't be defined as a material with a priority order and a buffer. "<< endl;
                    exit(1);

            }
    }

    fStreamListFPIsBuffer[keyword]        = 1;
    fInfiniteMaterialFromList[keyword]    = true;
} 

//________________________________________________________________________
void FabricationPlant::AddFuelBuffer(string keyword, Storage* Stock)
{
    fStorage[keyword].push_back(Stock);
    
     //Test if there is no buffer already defined//
    map < string , bool >::iterator it_s_B;
    for( it_s_B = fStreamListFPIsBuffer.begin();  it_s_B != fStreamListFPIsBuffer.end(); it_s_B++)
    { 
            if ((*it_s_B).second == 1)
            {
                    ERROR << (*it_s_B).first<<" is already defined as a buffer. "<< endl;
                    ERROR << " Fuel can't have more than one buffer in current algorithm. "<< endl;      
                    exit(1);
            }
    }
    
    map < int , string >::iterator it_i_s;
    for( it_i_s = fStreamListFPPriority.begin();  it_i_s != fStreamListFPPriority.end(); it_i_s++)
    { 
            if ((*it_i_s).second == keyword)
            {
                    ERROR << "The material "<<keyword<<" can't be defined as a material with a priority order and a buffer. "<< endl;
                    exit(1);

            }
    }

    fStreamListFPIsBuffer[keyword]          = 1;
    fInfiniteMaterialFromList[keyword]      = true;
} 

