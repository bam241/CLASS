#include "EQM_QUAD_PWR_MOX.hxx"

#include <vector>

#include "StringLine.hxx"
#include "LogFile.hxx"




EQM_QUAD_PWR_MOX::EQM_QUAD_PWR_MOX(string WeightPath):EquivalenceModel()
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


EQM_QUAD_PWR_MOX::~EQM_QUAD_PWR_MOX()
{

}




double EQM_QUAD_PWR_MOX::GetFissileMolarFraction(IsotopicVector Fissile,IsotopicVector Fertile,double BurnUp)
{


	ZAI ZAIList[6] = {ZAI(94,238,0), ZAI(94,239,0), ZAI(94,240,0), ZAI(94,241,0), ZAI(94,242,0), ZAI(95,241,0)  };

	vector<double> PuCompo;
	double Sum = Fissile.GetSumOfAll();


	for(int i = 0; i< 5; i++)
		PuCompo.push_back( Fissile.GetZAIIsotopicQuantity(ZAIList[i])/Sum *100 );

	PuCompo[2] += Fissile.GetZAIIsotopicQuantity(ZAIList[5])/Sum *100;
	double A = 0;

	if(PuCompo[0] <= PuCompo[2] && PuCompo[0] <= PuCompo[4] && PuCompo[1] + PuCompo[3] >= 40 && PuCompo[1] >0 )
	{
		int par = 0;
		for(int j = 0 ; j < 5 ; j++)
		{
			A += fFuelParameter[par]   * PuCompo[j]/100 ;
			par++;
			for(int i = j ; i < 5 ; i++)
			{
				A += fFuelParameter[par] *PuCompo[i]/100 *PuCompo[j]/100;
				par++;
			}
		}
		A += fFuelParameter[par];
	}

	return A;
}

