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
    freaded = false;
}
//________________________________________________________________________
EquivalenceModel::EquivalenceModel(CLASSLogger* log): CLASSObject(log)
{
    freaded = false;
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

    if(Stocks.size == 0)
    {
        return;
    }
    if(Stocks.size != lambda.size())
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
            double ThatZAIMin   = Domain_it->second.first;
            double ThatZAIMax   = Domain_it->second.second;
            if ( (ThatZAIProp > ThatZAIMax) || (ThatZAIProp <   ThatZAIMin) )
            {
                IsInDomain = false;

                WARNING << "Fresh fuel out of model range" << endl;
                WARNING << "\t AT LEAST this ZAI is accused to be outrange :" << endl;
                WARNING << "\t\t" << Domain_it->first.Z() << " " << Domain_it->first.A() << " " << Domain_it->first.I() << endl;
                WARNING << "\t\t min = " << ThatZAIMin   << " value = " << ThatZAIProp << " max = " << ThatZAIMax << endl;
                WARNING << "\t IV accused :" << endl << endl;
                WARNING << IVNorm.sPrint() << endl;
                break;
            }
        }
    }
    DBGL
    return IsInDomain;

}
