#include "ZAITox.hxx"

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
ZAITox::ZAITox()
{
	string  CLASSPATH = getenv("CLASS_PATH");
	string	HeatToxFile = CLASSPATH + "/data/HeatTox.dat";

	ifstream infile(HeatToxFile.c_str());	

	if(!infile.good())
	{	
		cout << "Error in ZAITox : can't find/open file " << HeatToxFile << endl;
		exit(1);
	}

	IrradiationModel* IM = new IrradiationModel();
	IM->SetZAIThreshold(1);
	IM->LoadDecay();

	int Z,A,I;
	string Name;
	double WperBq,SvperBq;
	while (!infile.eof())
	{	
		string line;
		stringstream ossline;
		getline(infile, line);
		ossline << line;
		
		if(StringLine::IsDouble(line.substr(0,1))) //else is a comment
		{	ossline>>Z>>A>>I>>WperBq>>SvperBq;
			fZAITox.insert( pair< ZAI,double >( ZAI(Z,A,I),IM->GetDecayConstant(ZAI(Z,A,I))*SvperBq) );
		}	
	}

	delete IM;
	infile.close();
}	
//________________________________________________________________________
ZAITox::~ZAITox()
{
	fZAITox.clear();
}
//________________________________________________________________________
double ZAITox::GetRadioTox(ZAI zai ) const
{

	map<ZAI ,double> ZAITox = fZAITox;
	map<ZAI ,double>::iterator it;
	it = ZAITox.find(zai);

	if ( it != ZAITox.end() )
	{
		return it->second;
	}
	else
	{
		return 0;
	}

}
//________________________________________________________________________
double ZAITox::GetRadioTox(const IsotopicVector IV) const
{
	double TotalTox = 0;
	
	map<ZAI ,double >::iterator it;
	map<ZAI ,double > isotopicquantity = IV.GetIsotopicQuantity();
	
	for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		TotalTox += (*it).second * GetRadioTox( (*it).first ) ;
	
	
	return TotalTox;
	
}
