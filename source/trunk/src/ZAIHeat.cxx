#include "ZAIHeat.hxx"

//#include "CLASSConstante.hxx"

#include "stdlib.h"
#include <fstream>
#include <string>

#include "IsotopicVector.hxx"
#include "StringLine.hxx"

#include "IrradiationModel.hxx"

	//________________________________________________________________________
	//
	//		ZAI
	//
	//
	//
	//
	//________________________________________________________________________
	//____________________________InClass Operator____________________________
	//________________________________________________________________________

//________________________________________________________________________
ZAIHeat::ZAIHeat()
{
	string  CLASSPATH = getenv("CLASS_PATH");
	string	HeatToxFile = CLASSPATH + "/data/HeatTox.dat";

	ifstream infile(HeatToxFile.c_str());	

	if(!infile.good())
	{	
		cout << "Error in ZAIHeat : can't find/open file " << HeatToxFile << endl;
		exit(1);
	}

	IrradiationModel* IM= new IrradiationModel();
	IM->LoadDecay();

	int Z,A,I;
	string Name;
	double WperBq,SvperBq;
	while (!infile.eof())
	{	
		string line;
		stringstream ossline;
		getline(infile, line);
		ossline<<line;
		
		if(StringLine::IsDouble(line.substr(0,1))) //else is a comment
		{	ossline>>Z>>A>>I>>WperBq>>SvperBq;
			fZAIHeat.insert( pair< ZAI,double >( ZAI(Z,A,I),IM->GetDecayConstant(ZAI(Z,A,I))*WperBq) );
		}	
	}

	delete IM;
	infile.close();
}	
//________________________________________________________________________
ZAIHeat::~ZAIHeat()
{
	fZAIHeat.clear();
}
//________________________________________________________________________
double ZAIHeat::GetHeat(ZAI zai ) const
{

	map<ZAI ,double> ZAIHeat = fZAIHeat;
	map<ZAI ,double>::iterator it;
	it = ZAIHeat.find(zai);

	if ( it != ZAIHeat.end() )
	{
		return it->second;
	}
	else
	{
		return 0;
	}

}
//________________________________________________________________________
double ZAIHeat::GetHeat(const IsotopicVector IV) const
{
	double TotalHeat = 0;
	
	map<ZAI ,double >::iterator it;
	map<ZAI ,double > isotopicquantity = IV.GetIsotopicQuantity();
	
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		TotalHeat += (*it).second * GetHeat( (*it).first ) ;
	
	
	return TotalHeat;
	
}
