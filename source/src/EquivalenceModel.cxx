#include "EquivalenceModel.hxx"
#include "StringLine.hxx"
#include "CLASSMethod.hxx"

#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "TSystem.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodCuts.h"

#include "CLASSReader.hxx"

//________________________________________________________________________
//________________________________________________________________________
EquivalenceModel::EquivalenceModel(): CLASSObject()
{
    EquivalenceModel::LoadKeyword();
    fread = false;
}
//________________________________________________________________________
EquivalenceModel::EquivalenceModel(CLASSLogger* log): CLASSObject(log)
{
    EquivalenceModel::LoadKeyword();
    fread = false;
}
//________________________________________________________________________
EquivalenceModel::~EquivalenceModel()
{

}
//________________________________________________________________________
map <string , vector<double> > EquivalenceModel::BuildFuel(double BurnUp, double HMMass, map < string , vector <IsotopicVector> > StreamArray,  map < string , double> StreamListMassFractionMin, map < string , double> StreamListMassFractionMax, map < int , string> StreamListPriority, map < string , bool> StreamListIsBuffer)
{

}

//________________________________________________________________________
void EquivalenceModel::LoadKeyword() 
{
    DBGL
    fKeyword.insert( pair<string, EQM_MthPtr>( "k_zail", & EquivalenceModel::ReadZAIlimits) );
    fKeyword.insert( pair<string, EQM_MthPtr>( "k_reactor", & EquivalenceModel::ReadType) );
    fKeyword.insert( pair<string, EQM_MthPtr>( "k_fuel", & EquivalenceModel::ReadType) );
    fKeyword.insert( pair<string, EQM_MthPtr>( "k_massfractionmin", & EquivalenceModel::ReadEqMinFraction) );
    fKeyword.insert( pair<string, EQM_MthPtr>( "k_massfractionmax", & EquivalenceModel::ReadEqMaxFraction) );
    fKeyword.insert( pair<string, EQM_MthPtr>( "k_list", & EquivalenceModel::ReadList) );
    fKeyword.insert( pair<string, EQM_MthPtr>( "k_specpower", & EquivalenceModel::ReadSpecificPower) );
    fKeyword.insert( pair<string, EQM_MthPtr>( "k_zainame", & EquivalenceModel::ReadZAIName) );
    fKeyword.insert( pair<string, EQM_MthPtr>( "k_maxburnup", & EquivalenceModel::ReadMaxBurnUp) ); 
    fKeyword.insert( pair<string, EQM_MthPtr>( "k_predictortype", & EquivalenceModel::ReadPredictorType) );
    fKeyword.insert( pair<string, EQM_MthPtr>( "k_output", & EquivalenceModel::ReadOutput) );
    fKeyword.insert( pair<string, EQM_MthPtr>( "k_buffer", & EquivalenceModel::ReadBuffer) ); 
    DBGL
}

//________________________________________________________________________
void EquivalenceModel::ReadNFO()
{
    DBGL
    std::cout << fTMVANFOFilePath << endl;
      ifstream NFO(fTMVANFOFilePath.c_str());
    
    if(!NFO)
    {
        ERROR << "Can't find/open file " << fTMVANFOFilePath << endl;
        exit(0);
    }
    
    do
    {
        string line;
        getline(NFO,line);

        EquivalenceModel::ReadLine(line);
        
    } while(!NFO.eof());
    
    DBGL
}

//________________________________________________________________________
void EquivalenceModel::ReadLine(string line)
{
    DBGL
    
    if (!fread)
    {
        int pos = 0;
        string keyword = tlc(StringLine::NextWord(line, pos, ' '));

        map<string, EQM_MthPtr>::iterator it = fKeyword.find(keyword);
        
        if(it != fKeyword.end())
            (this->*(it->second))( line );
        
        fread = true;
        ReadLine(line);
        
    }
    
    fread = false;
    
    DBGL
}
//________________________________________________________________________
void EquivalenceModel::ReadType(const string &line)
{
    DBGL
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_fuel" && keyword != "k_reactor" ) // Check the keyword
    {
        ERROR << " Bad keyword : " << keyword << " Not found !" << endl;
        exit(1);
    }
    if( keyword ==  "k_fuel" )
        fDBFType = StringLine::NextWord(line, pos, ' ');
    else if( keyword ==  "k_reactor" )
        fDBRType = StringLine::NextWord(line, pos, ' ');
    
    DBGL
}
//________________________________________________________________________
void EquivalenceModel::ReadZAIlimits(const string &line)
{
    DBGL
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_zail" )   // Check the keyword
    {
        ERROR << " Bad keyword : \"k_zail\" not found !" << endl;
        exit(1);
    }
    
    int Z   = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    int A   = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    int I   = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    
    double downLimit    = atof(StringLine::NextWord(line, pos, ' ').c_str());
    double upLimit  = atof(StringLine::NextWord(line, pos, ' ').c_str());
    
    if (upLimit < downLimit)
    {
        double tmp  = upLimit;
        upLimit     = downLimit;
        downLimit   = tmp;
    }
    fZAILimits.insert(pair<ZAI, pair<double, double> >(ZAI(Z,A,I), pair<double,double>(downLimit, upLimit)));
    DBGL
}
//________________________________________________________________________
void EquivalenceModel::ReadList(const string &line)
{
    DBGL
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_list" )   // Check the keyword
    {
        ERROR << " Bad keyword : \"k_list\" not found !" << endl;
        exit(1);
    }
    string ListName= StringLine::NextWord(line, pos, ' ');
    int Z       = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    int A       = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    int I       = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    double Q    = atof(StringLine::NextWord(line, pos, ' ').c_str());
    fStreamList[ListName].Add(Z, A, I, Q);
    
    DBGL
}
//________________________________________________________________________
void EquivalenceModel::ReadEqMinFraction(const string &line)
{
    DBGL
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_massfractionmin" )    // Check the keyword
    {
        ERROR << " Bad keyword : \"k_massfractionmin\" not found !" << endl;
        exit(1);
    }
    string ListName= StringLine::NextWord(line, pos, ' ');
    double Q     = atof(StringLine::NextWord(line, pos, ' ').c_str());
    fStreamListEqMMassFractionMin[ListName] = Q;

    DBGL
}

//________________________________________________________________________
void EquivalenceModel::ReadEqMaxFraction(const string &line)
{
    DBGL
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_massfractionmax" )    // Check the keyword
    {
        ERROR << " Bad keyword : \"k_massfractionmax\" not found !" << endl;
        exit(1);
    }
    string ListName= StringLine::NextWord(line, pos, ' ');
    double Q     = atof(StringLine::NextWord(line, pos, ' ').c_str());
    fStreamListEqMMassFractionMax[ListName] = Q;

    DBGL
}

//________________________________________________________________________
void EquivalenceModel::ReadSpecificPower(const string &line)
{
    DBGL
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_specpower")   // Check the keyword
    {
        ERROR << " Bad keyword : \"k_specpower\" Not found !" << endl;
        exit(1);
    }
    
    fSpecificPower = atof(StringLine::NextWord(line, pos, ' ').c_str());
    
    DBGL
}
//________________________________________________________________________
void EquivalenceModel::ReadZAIName(const string &line)
{
    DBGL
    
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_zainame" )    // Check the keyword
    {
        ERROR << " Bad keyword : \"k_zainame\" not found !" << endl;
        exit(1);
    }
    
    int Z = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    int A = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    int I = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    
    string name = StringLine::NextWord(line, pos, ' ');
    
    fMapOfTMVAVariableNames.insert( pair<ZAI,string>( ZAI(Z, A, I), name ) );

    DBGL    
}
//________________________________________________________________________
void EquivalenceModel::ReadMaxBurnUp(const string &line)
{
    DBGL
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_maxburnup" )  // Check the keyword
    {
        ERROR << " Bad keyword : \"k_maxburnup\" not found !" << endl;
        exit(1);
    }
    
    fMaximalBU = atof(StringLine::NextWord(line, pos, ' ').c_str());

    DBGL
}

//________________________________________________________________________
void EquivalenceModel::SetLambdaToErrorCode(vector<double>& lambda)
{
    DBGL
    if (lambda.size() == 0) //then we have to add an element to send the error code to the fab (case for no storage in stream)
    {
        lambda.push_back(-1);
    }

    else // other errors (no enough material or too many steps)
    {
        for ( int i = 0; i < (int)lambda.size(); i++)
        {
            lambda[i] = -1;
        }
    }
    DBGL
}
//________________________________________________________________________
void EquivalenceModel::ReadPredictorType(const string &line)
{
    DBGL
    
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_predictortype" )  // Check the keyword
    {
        ERROR << " Bad keyword : \"k_predictortype\" not found !" << endl;
        exit(1);
    }
        
    fPredictorType = StringLine::NextWord(line, pos, ' ');
    
    DBGL    
}
//________________________________________________________________________
void EquivalenceModel::ReadOutput(const string &line)
{
    DBGL
    
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_output" ) // Check the keyword
    {
        ERROR << " Bad keyword : \"k_output\" not found !" << endl;
        exit(1);
    }
        
    fOutput = StringLine::NextWord(line, pos, ' ');
    
    DBGL    
}
//________________________________________________________________________
void EquivalenceModel::ReadBuffer(const string &line)
{
    DBGL
    
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_buffer" ) // Check the keyword
    {
        ERROR << " Bad keyword : \"k_buffer\" not found !" << endl;
        exit(1);
    }
        
    fBuffer = StringLine::NextWord(line, pos, ' ');
    
    DBGL    
}


//________________________________________________________________________
void EquivalenceModel::StocksTotalMassCalculation(map < string , vector <IsotopicVector> > const& Stocks)
{
    DBGL
    // Calculating total mass of stock once and for all
    double TotalMassInStocks = 0;
    map < string , vector <IsotopicVector> >::const_iterator it_s_vIV;

    for ( it_s_vIV = Stocks.begin();  it_s_vIV != Stocks.end(); it_s_vIV++)
    {
        fTotalMassInStocks[ it_s_vIV->first ] = 0;
        fLambdaMax[ it_s_vIV->first ] = 0;
    }
    for (  it_s_vIV = Stocks.begin();   it_s_vIV != Stocks.end();  it_s_vIV++)
    {
        TotalMassInStocks = 0;
        for (int i = 0; i < (int)Stocks.at((* it_s_vIV).first).size(); i++)
        {
            TotalMassInStocks  +=  Stocks.at( it_s_vIV->first )[i].GetTotalMass();
        }
        fLambdaMax[(*it_s_vIV).first] = Stocks.at( it_s_vIV->first ).size();
        fTotalMassInStocks[ it_s_vIV->first ] = TotalMassInStocks * 1e6; // in grams
    }
    DBGL
}

//________________________________________________________________________
void EquivalenceModel::ConvertMassToLambdaVector(string MaterialDenomination, vector<double>& lambda, double MaterialMassNeeded, vector <IsotopicVector>  Stocks)
{
    DBGL

    if (Stocks.size() == 0)
    {
        return;
    }
    if (Stocks.size() != lambda.size())
    {
        ERROR << "Stocks vector size MUST be the same as lamba size!!" << endl;
        exit(1);
    }
    double Lambda_tot = 0;

    // Calculation of Lambda tot associated to the required mass MaterialMassNeeded
    for ( int i = 0; i < (int)Stocks.size(); i++)
    {
        if ( MaterialMassNeeded >= (Stocks[i].GetTotalMass() * 1e6))
        {
            Lambda_tot +=  1;
            MaterialMassNeeded -=  (Stocks[i].GetTotalMass() * 1e6);
        }
        else
        {
            Lambda_tot +=  MaterialMassNeeded / (Stocks[i].GetTotalMass() * 1e6);
            break;
        }
    }
    // Calculate lambda vector associated to the lambda tot
    if (Lambda_tot > (int)lambda.size() )
    {
        cout << Lambda_tot << "  " << lambda.size() << endl;
        ERROR << " FATAL ERROR " << endl;
        exit(0);
    }

    for (int i = 0 ; i < (int)lambda.size() ; i++) //set to 0 all non touched value (to be sure)
        lambda[i] = 0  ;

    int IntegerPart         = floor( Lambda_tot );
    double DecimalPart  = Lambda_tot - IntegerPart;

    for (int i = 0  ; i < IntegerPart; i++ )
        lambda[i] = 1;

    lambda[IntegerPart] = DecimalPart;
    DBGL
}
//________________________________________________________________________
bool EquivalenceModel::isIVInDomain(IsotopicVector IV)
{
    DBGL
    bool IsInDomain = true;

    if (fZAILimits.empty())
    {
        WARNING << "Fresh Fuel variation domain is not set" << endl;
        WARNING << "CLASS has no clue if the computed evolution for this fresh fuel is correct" << endl;
        WARNING << "Proceed finger crossed !!" << endl;
        return true;
    }

    else
    {
        IsotopicVector IVNorm = IV / IV.GetSumOfAll();
        for (map< ZAI, pair<double, double> >::iterator Domain_it = fZAILimits.begin(); Domain_it != fZAILimits.end(); Domain_it++)
        {
            double ThatZAIProp = IVNorm.GetIsotopicQuantity()[Domain_it->first];
            double ThatZAIMin = Domain_it->second.first;
            double ThatZAIMax = Domain_it->second.second;
            if ( (ThatZAIProp > ThatZAIMax) || (ThatZAIProp < ThatZAIMin) )
            {
                IsInDomain = false;

                WARNING << "Fresh fuel out of model range" << endl;
                WARNING << "\t AT LEAST this ZAI is accused to be outrange :" << endl;
                WARNING << "\t\t" << Domain_it->first.Z() << " " << Domain_it->first.A() << " " << Domain_it->first.I() << endl;
                WARNING << "\t\t min = " << ThatZAIMin << " value = " << ThatZAIProp << " max = " << ThatZAIMax << endl;
                WARNING << "\t IV accused :" << endl << endl;
                WARNING << IVNorm.sPrint() << endl;
                break;
            }
        }
    }
    DBGL
    return IsInDomain;

}
