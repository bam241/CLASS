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


using namespace std;



XSModel::XSModel():CLASSObject()
{
	LoadKeyword();
}


XSModel::XSModel(CLASSLogger* log):CLASSObject(log)
{
	LoadKeyword();
}	

void XSModel::LoadKeyword()
{
	DBGL
	fKeyword.insert( pair<string, MthPtr>( "k_zail",	& XSModel::ReadZAIlimits));
	fKeyword.insert( pair<string, MthPtr>( "k_reactor",	& XSModel::ReadType)	);
	fKeyword.insert( pair<string, MthPtr>( "k_fuel",	& XSModel::ReadType)	);
	fKeyword.insert( pair<string, MthPtr>( "k_mass",	& XSModel::ReadRParam)	);
	fKeyword.insert( pair<string, MthPtr>( "k_power",	& XSModel::ReadRParam)	);
	DBGL
}

void XSModel::ReadRParam(const string &line)
{
	DBGL
	int start = 0;
	string type = tlc(StringLine::NextWord(line, start, ' '));
	if( type != "k_power" || type != "k_mass" )	// Check the keyword
	{
		ERROR << " Bad keyword : " << type << " Not found !" << endl;
		exit(1);
	}
	if( type == "k_mass" )
		fDBHMMass = atof(StringLine::NextWord(line, start, ' ').c_str());
	else if( type == "k_mass" )
		fDBPower = atof(StringLine::NextWord(line, start, ' ').c_str());
	
	DBGL
}


void XSModel::ReadType(const string &line)
{
	DBGL
	int start = 0;
	string type = tlc(StringLine::NextWord(line, start, ' '));
	if( type != "k_fuel" || type != "k_reactor" )	// Check the keyword
	{
		ERROR << " Bad keyword : " << type << " Not found !" << endl;
		exit(1);
	}
	if( type == "k_fuel" )
		fDBFType = StringLine::NextWord(line, start, ' ');
	else if( type == "k_reactor" )
		fDBRType = StringLine::NextWord(line, start, ' ');

	DBGL
}


void XSModel::ReadZAIlimits(const string &line)
{
	DBGL
	int start = 0;
	if( tlc(StringLine::NextWord(line, start, ' ')) != "k_zail" )	// Check the keyword
	{
		ERROR << " Bad keyword : \"k_zail\" not found !" << endl;
		exit(1);
	}
	
	int Z = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int A = atoi(StringLine::NextWord(line, start, ' ').c_str());
	int I = atoi(StringLine::NextWord(line, start, ' ').c_str());
	
	double upLimit = atof(StringLine::NextWord(line, start, ' ').c_str());
	double downLimit = atof(StringLine::NextWord(line, start, ' ').c_str());
	
	DBGL
}


bool XSModel::isIVInDomain(IsotopicVector IV)
{
DBGL
	bool IsInDomain=true;

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
		for (map< ZAI,pair<double,double> >::iterator Domain_it=fZAILimits.begin(); Domain_it!=fZAILimits.end(); Domain_it++)
		{		
			double ThatZAIProp = IVNorm.GetIsotopicQuantity()[Domain_it->first]	;
			double ThatZAIMin  = Domain_it->second.first;
			double ThatZAIMax  = Domain_it->second.second;
			if( (ThatZAIProp > ThatZAIMax) || (ThatZAIProp <  ThatZAIMin) )
			{ 	
				IsInDomain = false;	

				WARNING<<"Fresh fuel out of model range"<<endl;
				WARNING<<"\t AT LEAST this ZAI is accused to be outrange :"<<endl;
				WARNING<<"\t\t"<<Domain_it->first.Z()<<" "<<Domain_it->first.A()<<" "<<Domain_it->first.I()<<endl;
				WARNING<<"\t\t min="<<ThatZAIMin<<" value="<<ThatZAIProp<<" max="<<ThatZAIMax<<endl;
				WARNING<<"\t IV accused :"<<endl<<endl;
				WARNING<<IVNorm.sPrint()<<endl;
				break;
			}
		}	
	}
DBGL
return IsInDomain;

}