//
//  XSModel.cxx
//  CLASSSource
//
//  Created by BaM on 04/05/2014.
//  Copyright (c) 2014 BLG. All rights reserved.
//

#include "XSModel.hxx"

using namespace std;

XSModel::XSModel():CLASSObject()
{

}


XSModel::XSModel(CLASSLogger* log):CLASSObject(log)
{

}	

bool XSModel::isIVInDomain(IsotopicVector IV)
{DBGL
	bool IsInDomain=true;

	if(fFreshFuelDomain.empty())
	{
	 WARNING << "Fresh Fuel variation domain is not set" << endl;
	 WARNING << "CLASS has no clue if the computed evolution for this fresh fuel is correct" << endl;
	 WARNING << "Proceed finger crossed !!" << endl;
	 return true;
	}

	else
	{ 
		IsotopicVector IVNorm = IV /IV.GetSumOfAll();
		for (map< ZAI,pair<double,double> >::iterator Domain_it=fFreshFuelDomain.begin(); Domain_it!=fFreshFuelDomain.end(); Domain_it++)
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