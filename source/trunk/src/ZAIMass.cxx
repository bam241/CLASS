#include "ZAIMass.hxx"
#include "CLASSHeaders.hxx"
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


ZAIMass::ZAIMass()
{
	string  CLASSPATH = getenv("CLASS_PATH");
	string	MassDataFile = CLASSPATH +"/source/data/Mass.dat";
	ifstream infile(MassDataFile.c_str());	
	if(!infile.good())
	{	
		cout<<"ZAIMass Error.\n can't find/open file "<<MassDataFile<<endl;
		exit(1);
	}

	int Z,A;
	string Name;
	double MassUnity,MassDec,error;
	while (infile>>Z>>A>>Name>>MassUnity>>MassDec>>error)
	{
		double Masse=MassUnity+MassDec*1e-6;
		fZAIMass.insert( pair< ZAI,double >( ZAI(Z,A,0), Masse ) );
	}

	infile.close();

}


ZAIMass::~ZAIMass()
{
	fZAIMass.clear();
}

