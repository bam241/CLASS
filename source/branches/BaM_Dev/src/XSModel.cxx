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
	XSModel::LoadKeyword();
}


XSModel::XSModel(CLASSLogger* log):CLASSObject(log)
{
	
	XSModel::LoadKeyword();
}


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
		
	} while(NFO.eof());
	
	DBGL
}

//________________________________________________________________________
void XSModel::ReadLine(string line)
{
	DBGL
	
	if (!freaded)
	{
		int start = 0;
		string keyword = tlc(StringLine::NextWord(line, start, ' '));
		
		(this->*fKeyword[ keyword ])(line);
		freaded = true;
		
		ReadLine(line);
		
	}
	
	freaded = false;
	
	DBGL
}


void XSModel::LoadKeyword()
{
	DBGL
	fKeyword.insert( pair<string, MthPtr>( "k_zail",	& XSModel::ReadZAIlimits));
	fKeyword.insert( pair<string, MthPtr>( "k_reactor",	& XSModel::ReadType)	 );
	fKeyword.insert( pair<string, MthPtr>( "k_fuel",	& XSModel::ReadType)	 );
	fKeyword.insert( pair<string, MthPtr>( "k_mass",	& XSModel::ReadRParam)	 );
	fKeyword.insert( pair<string, MthPtr>( "k_power",	& XSModel::ReadRParam)	 );
	DBGL
}

void XSModel::ReadRParam(const string &line)
{
	DBGL
	int start = 0;
	string keyword = tlc(StringLine::NextWord(line, start, ' '));
	if( keyword != "k_power" || keyword != "k_mass" )	// Check the keyword
	{
		ERROR << " Bad keyword : " << keyword << " Not found !" << endl;
		exit(1);
	}
	if( keyword == "k_mass" )
		fDBHMMass = atof(StringLine::NextWord(line, start, ' ').c_str());
	else if( keyword == "k_mass" )
		fDBPower = atof(StringLine::NextWord(line, start, ' ').c_str());
	
	DBGL
}


void XSModel::ReadType(const string &line)
{
	DBGL
	int start = 0;
	string keyword = tlc(StringLine::NextWord(line, start, ' '));
	if( keyword != "k_fuel" || keyword != "k_reactor" )	// Check the keyword
	{
		ERROR << " Bad keyword : " << keyword << " Not found !" << endl;
		exit(1);
	}
	if( keyword == "k_fuel" )
		fDBFType = StringLine::NextWord(line, start, ' ');
	else if( keyword == "k_reactor" )
		fDBRType = StringLine::NextWord(line, start, ' ');
	
	DBGL
}


void XSModel::ReadZAIlimits(const string &line)
{
	DBGL
	int start = 0;
	string keyword = tlc(StringLine::NextWord(line, start, ' '));
	if( keyword != "k_zail" )	// Check the keyword
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