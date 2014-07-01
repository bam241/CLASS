#include "EQM_LIN_PWR_MOX.hxx"

#include <vector>

#include "StringLine.hxx"
#include "LogFile.hxx"





EQM_LIN_PWR_MOX::EQM_LIN_PWR_MOX(string WeightPath):EquivalenceModel()
{
	fWeightPath =  WeightPath;

	ifstream DataDB(fWeightPath.c_str());							// Open the File
	if(!DataDB)
		cout << "!!Warning!! !!!EQM QUAD PWR MOX!!! \n Can't open \"" << fWeightPath << "\"\n" << endl;

	string line;
	int start = 0;	// First Get Fuel Parameter
	getline(DataDB, line);

	if( StringLine::NextWord(line, start, ' ') != "PARAM")
	{
		cout << "!!Bad Trouble!! !!!EQM QUAD PWR MOX!!! Bad Database file : " <<  fWeightPath << " Can't find the Parameter of the DataBase"<< endl;
		exit (1);
	}
	while(start < (int)line.size())
		fFuelParameter.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));

	cout << "!!INFO!! !!!EQM QUAD PWR MOX!!! " <<  fFuelParameter.size() << " have been read"<< endl;



	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	// ADD ENrichment of the U reading !!!!!!!!!!!!!!!!!!!!!!!!!!!!		       //
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//

	ZAI U8(92,238,0);
	ZAI U5(92,235,0);
	double U5_enrich= 0.0025;
	fFertileList = U5*U5_enrich + U8*(1-U5_enrich);


	ZAI Pu8(94,238,0);
	ZAI Pu9(94,239,0);
	ZAI Pu0(94,240,0);
	ZAI Pu1(94,241,0);
	ZAI Pu2(94,242,0);
	fFissileList = Pu8*1+Pu9*1+Pu0*1+Pu1*1+Pu2*1;
	
	
}

EQM_LIN_PWR_MOX::~EQM_LIN_PWR_MOX()
{

}

//________________________________________________________________________
vector<double> EQM_LIN_PWR_MOX::BuildFuel(double BurnUp, double HMMass,vector<IsotopicVector> FissilArray, vector<IsotopicVector> FertilArray)
{

	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	// ADD ENrichment of the U check !!!!!!!!!!!!!!!!!!!!!!!!!!!!		       //
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//


	vector<double> lambda;

	return lambda;
}
