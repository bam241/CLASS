#include "EQM_LIN_PWR_MOX.hxx"

#include "CLASSConstante.hxx"

#include <vector>

#include "StringLine.hxx"
#include "CLASSLogger.hxx"
#include "IsotopicVector.hxx"



EQM_LIN_PWR_MOX::EQM_LIN_PWR_MOX(string WeightPath):EquivalenceModel(new CLASSLogger("EQM_LIN_PWR_MOX.log"))
{
	fWeightPath =  WeightPath;

	ifstream DataDB(fWeightPath.c_str());							// Open the File
	if(!DataDB)
		WARNING << "Can't open \"" << fWeightPath << "\"\n" << endl;

	string line;
	int start = 0;	// First Get Fuel Parameter
	getline(DataDB, line);

	if( StringLine::NextWord(line, start, ' ') != "PARAM")
	{
		ERROR << " Bad Database file : " <<  fWeightPath << " Can't find the Parameter of the DataBase " << endl;
		exit (1);
	}
	while(start < (int)line.size())
		fFuelParameter.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));

	INFO << " " << (int)fFuelParameter.size() << " parameters have been read " << endl;



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


EQM_LIN_PWR_MOX::EQM_LIN_PWR_MOX(CLASSLogger* log, string WeightPath):EquivalenceModel(log)
{
	fWeightPath =  WeightPath;

	ifstream DataDB(fWeightPath.c_str());							// Open the File
	if(!DataDB)
		WARNING << " Can't open \"" << fWeightPath << "\"\n" << endl;

	string line;
	int start = 0;	// First Get Fuel Parameter
	getline(DataDB, line);

	if( StringLine::NextWord(line, start, ' ') != "PARAM")
	{
		ERROR << " Bad Database file : " <<  fWeightPath << " Can't find the Parameter of the DataBase"<< endl;
		exit (1);
	}
	while(start < (int)line.size())
		fFuelParameter.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));

	INFO << fFuelParameter.size() << " have been read"<< endl;



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
	// ADD ENrichment of the U check ++ Check Un seul fertile !!!!		       //
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//
	//-----------------------------------------------------------------------------//


	vector<double> lambda;
	for(int i = 0; i < (int) (FissilArray.size() + FertilArray.size()); i++)
		lambda.push_back(0);


	double Na = 6.02214129e23;	//N Avogadro

	IsotopicVector FullUsedStock;
	IsotopicVector stock;

	bool FuelBuild = false;
	if(FissilArray.size() == 0)
	{
		for(int i = 0; i < (int)lambda.size(); i++)
			lambda[i] = -1;

		FuelBuild = true;
	}
	int N_FissilStock_OnCheck = 0;

	while(!FuelBuild)
	{

		double nPu_0 = 0;
		double MPu_0 = 0;
		{
			map<ZAI ,double>::iterator it;

			map<ZAI ,double> isotopicquantity = FullUsedStock.GetSpeciesComposition(94).GetIsotopicQuantity();
			for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
				nPu_0 += (*it).second;

			isotopicquantity = (FullUsedStock.GetSpeciesComposition(94) + ZAI(94,241,0)*FullUsedStock.GetZAIIsotopicQuantity(95,241,0)).GetIsotopicQuantity(); //Add the 241Am as 241Pu... the Pu is not old in the Eq Model but is in the FissileArray....;
			for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
				MPu_0 += (*it).second*cZAIMass.fZAIMass.find( (*it).first )->second/Na*1e-6;
		}
		stock = FissilArray[N_FissilStock_OnCheck];
		double nPu_1 = 0;
		double MPu_1 = 0;
		double Sum_AlphaI_nPuI = 0;
		double Sum_AlphaI_nPuI0 = 0;
		{
			map<ZAI ,double>::iterator it;
			map<ZAI ,double> isotopicquantity = stock.GetSpeciesComposition(94).GetIsotopicQuantity();

			for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
			{
				if ((*it).first.A() >= 238 && (*it).first.A() <= 242)
				{
					nPu_1 += (*it).second;
					Sum_AlphaI_nPuI += fFuelParameter[(*it).first.A() -237]*(*it).second;
				}
			}

			isotopicquantity = (stock.GetSpeciesComposition(94) + ZAI(94,241,0)*stock.GetZAIIsotopicQuantity(95,241,0)).GetIsotopicQuantity(); //Add the 241Am as 241Pu... the Pu is not old in the Eq Model but is in the FissileArray....
			for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
				if ((*it).first.A() >= 238 && (*it).first.A() <= 242)
				{
					MPu_1 += (*it).second * (cZAIMass.fZAIMass.find( (*it).first )->second)/Na*1e-6;
				}

			isotopicquantity = FullUsedStock.GetSpeciesComposition(94).GetIsotopicQuantity();
			for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
				if ((*it).first.A() >= 238 && (*it).first.A() <= 242)
				{
					Sum_AlphaI_nPuI0 += fFuelParameter[(*it).first.A() -237]*(*it).second;
				}
		}

		double StockFactionToUse = 0;

		double NT = HMMass*1e6 * Na / (cZAIMass.GetMass( ZAI(92,238,0) ) * 0.997
					       + cZAIMass.GetMass( ZAI(92,235,0) ) * 0.003 );

		double N1 = (BurnUp - fFuelParameter[6]) * NT;
		double N2 = -Sum_AlphaI_nPuI0;
		double N3 = -fFuelParameter[0] * Na / (cZAIMass.fZAIMass.find( ZAI(92,238,0) )->second*0.997
						       + cZAIMass.fZAIMass.find( ZAI(92,235,0) )->second*0.003 )
		* (HMMass*1e6 - MPu_0*1e6);

		double D1 = Sum_AlphaI_nPuI;
		double D2 = -fFuelParameter[0] * MPu_1*1e6 * Na / (cZAIMass.fZAIMass.find( ZAI(92,238,0) )->second*0.997
								   + cZAIMass.fZAIMass.find( ZAI(92,235,0) )->second*0.003 ) ;

		StockFactionToUse = (N1 + N2 + N3) / (D1 + D2);

		if(StockFactionToUse < 0)
		{
			WARNING << "!!!FabricationPlant!!! Oups Bug in calculating stock fraction to use "<< endl;
			lambda[N_FissilStock_OnCheck] = 0.;
			N_FissilStock_OnCheck++;
			FuelBuild = false;
		}
		else if( StockFactionToUse > 1 )
		{

			FullUsedStock += stock;
			lambda[N_FissilStock_OnCheck] = 1;
			N_FissilStock_OnCheck++;
			FuelBuild = false;
		}
		else
		{
			lambda[N_FissilStock_OnCheck] = StockFactionToUse;

			FuelBuild = true;

			double U8_Quantity = (HMMass - (MPu_0+StockFactionToUse*MPu_1 ))/(cZAIMass.fZAIMass.find( ZAI(92,238,0) )->second*0.997 + cZAIMass.fZAIMass.find( ZAI(92,235,0) )->second*0.003 )*Na/1e-6;

			lambda.back() = U8_Quantity / FertilArray[0].GetSumOfAll();
		}


		if( N_FissilStock_OnCheck == (int) FissilArray.size() )	// Check if the last Fissil stock has been tested... quit if so...
		{
			for(int i = 0; i < (int)lambda.size(); i++)
				lambda[i] = -1;
			
			FuelBuild = true;
		}
	}
	
	
	
	return lambda;
}
