//
//  XSModel.cxx
//  CLASSSource
//
//  Created by BaM on 04/05/2014.
//  Copyright (c) 2014 BLG. All rights reserved.
//

#include "XSModel.hxx"
#include "StringLine.hxx"
#include "CLASSMethod.hxx"
#include "ZAI.hxx"


using namespace std;



XSModel::XSModel():CLASSObject()
{
    DBGL
    XSModel::LoadKeyword();
    DBGL
    fread = false;
    DBGL

}
//________________________________________________________________________

XSModel::XSModel(CLASSLogger* log):CLASSObject(log)
{
    
    DBGL
    XSModel::LoadKeyword();
    DBGL
    fread = false;
}
//________________________________________________________________________
XSModel::~XSModel()
{

}
//________________________________________________________________________
void XSModel::ReadNFO()
{
    DBGL
    ifstream NFO(fInformationFile.c_str());
    
    if(!NFO)
    {
        ERROR << "Can't find/open file " << fInformationFile << endl;
        exit(0);
    }
    
    do
    {
        string line;
        getline(NFO,line);

        XSModel::ReadLine(line);
        
    } while(!NFO.eof());
    
    DBGL
}

//________________________________________________________________________
void XSModel::ReadLine(string line)
{
    DBGL
    
    if (!fread)
    {
        int pos = 0;
        string keyword = tlc(StringLine::NextWord(line, pos, ' '));

        map<string, XSM_MthPtr>::iterator it = fKeyword.find(keyword);
        
        if(it != fKeyword.end())
            (this->*(it->second))( line );
        
        fread = true;
        ReadLine(line);
        
    }
    
    fread = false;
    
    DBGL
}

//________________________________________________________________________
void XSModel::LoadKeyword()
{
    DBGL
    fKeyword.insert( pair<string, XSM_MthPtr>( "k_zail",    & XSModel::ReadZAIlimits));
    fKeyword.insert( pair<string, XSM_MthPtr>( "k_reactor",    & XSModel::ReadType)     );
    fKeyword.insert( pair<string, XSM_MthPtr>( "k_fuel",    & XSModel::ReadType)     );
    fKeyword.insert( pair<string, XSM_MthPtr>( "k_mass",    & XSModel::ReadRParam)     );
    fKeyword.insert( pair<string, XSM_MthPtr>( "k_power",    & XSModel::ReadRParam)     );
    DBGL
}

//________________________________________________________________________
void XSModel::ReadRParam(const string &line)
{
    DBGL
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_power" && keyword != "k_mass" )    // Check the keyword
    {
        ERROR << " Bad keyword : " << keyword << " Not found !" << endl;
        exit(1);
    }
    if( keyword ==  "k_mass" )
        fDBHMMass = atof(StringLine::NextWord(line, pos, ' ').c_str());
    else if( keyword ==  "k_power" )
        fDBPower = atof(StringLine::NextWord(line, pos, ' ').c_str());
    
    DBGL
}

//________________________________________________________________________

void XSModel::ReadType(const string &line)
{
    DBGL
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_fuel" && keyword != "k_reactor" )    // Check the keyword
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
void XSModel::ReadZAIlimits(const string &line)
{
    DBGL
    int pos = 0;
    string keyword = tlc(StringLine::NextWord(line, pos, ' '));
    if( keyword != "k_zail" )    // Check the keyword
    {
        ERROR << " Bad keyword : \"k_zail\" not found !" << endl;
        exit(1);
    }
    
    int Z = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    int A = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    int I = atoi(StringLine::NextWord(line, pos, ' ').c_str());
    
    double downLimit = atof(StringLine::NextWord(line, pos, ' ').c_str());
    double upLimit = atof(StringLine::NextWord(line, pos, ' ').c_str());

    if (upLimit < downLimit)
    {
        double tmp = upLimit;
        upLimit = downLimit;
        downLimit = tmp;
    }
    fZAILimits.insert(pair<ZAI, pair<double, double> >(ZAI(Z,A,I), pair<double,double>(downLimit, upLimit)));
    DBGL
}

//________________________________________________________________________
bool XSModel::isIVInDomain(IsotopicVector IV)
{
    DBGL
    bool IsInDomain = true;
    
    if(fZAILimits.empty())
    {
     WARNING << "Fresh Fuel variation domain is not set" << endl;
     WARNING << "CLASS has no clue if the computed evolution for this fresh fuel is correct" << endl;
     WARNING << "Proceed finger crossed !!" << endl;
     return true;
    }
    
    else
    {
        IsotopicVector IVNorm = IV /IV.GetSumOfAll();
        for (map< ZAI,pair<double,double> >::iterator Domain_it = fZAILimits.begin(); Domain_it != fZAILimits.end(); Domain_it++)
        {
            double ThatZAIProp = IVNorm.GetIsotopicQuantity()[Domain_it->first]    ;
            double ThatZAIMin  = Domain_it->second.first;
            double ThatZAIMax  = Domain_it->second.second;
            if( (ThatZAIProp > ThatZAIMax) || (ThatZAIProp <  ThatZAIMin) )
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
