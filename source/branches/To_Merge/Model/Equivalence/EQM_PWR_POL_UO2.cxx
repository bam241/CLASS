#include "EquivalenceModel.hxx"
#include "EQM_PWR_POL_UO2.hxx"
#include "CLASSLogger.hxx"
#include "StringLine.hxx"


// ________________________________________________________________________
// EQM_PWR_POL_UO2
//
// ________________________________________________________________________




//Constructor(s)
EQM_PWR_POL_UO2::EQM_PWR_POL_UO2(string PathToWeightFile):EquivalenceModel(new CLASSLogger("EQM_PWR_POL_UO2.log"))
{
	
	// Fertile
	ZAI U8(92 ,238 ,0) ;
	// Fissile
	ZAI U5(92 ,235 ,0) ;
	// ...

	fStreamList["Fissile"] = U5*1;
	fStreamList["Fertile"] = U8*1;

	ReadWeightFile(PathToWeightFile);

}
// _______________________________________________________________________
EQM_PWR_POL_UO2::EQM_PWR_POL_UO2(CLASSLogger* log,string PathToWeightFile):EquivalenceModel(log)
{
	
	// Fertile
	ZAI U8(92 ,238 ,0) ;
	// Fissile
	ZAI U5(92 ,235 ,0) ;
	// ...

	fStreamList["Fissile"] = U5*1;
	fStreamList["Fertile"] = U8*1;

	ReadWeightFile(PathToWeightFile);

}
// _______________________________________________________________________
void EQM_PWR_POL_UO2::ReadWeightFile(string PathToWeightFile)
{
	ifstream DataDB(PathToWeightFile.c_str());							// Open the File
	if(!DataDB)
		WARNING << " Can't open \"" << PathToWeightFile << "\"" << endl;
	
	string line;
	int start = 0;	// First Get Fuel Parameter
	getline(DataDB, line);
	
	if( StringLine::NextWord(line, start, ' ') != "PARAM")
	{
		ERROR << "Bad Database file : " <<  PathToWeightFile << " Can't find the Parameter of the DataBase " << endl;
		exit (1);
	}
	fParam_Bu_0 = atof(StringLine::NextWord(line, start, ' ').c_str()) ;
	fParam_Bu = atof(StringLine::NextWord(line, start, ' ').c_str());
	fParam_BuSquare = atof(StringLine::NextWord(line, start, ' ').c_str());
	
	INFO << "Weight parameters has been read " << endl;
	INFO  << "\t U enrichment = " << fParam_Bu_0 << " + " << fParam_Bu << "*Burnup + " <<  fParam_BuSquare << "*Burnup*Burnup" << endl;
}
// _______________________________________________________________________
map < string , double>  EQM_PWR_POL_UO2::GetMolarFraction( map < string , IsotopicVector> IVStream , double BurnUp )
{
	map < string , double> MolarFraction;

	double U5 = fParam_Bu_0 + fParam_Bu*BurnUp + fParam_BuSquare*BurnUp*BurnUp;

	MolarFraction["Fissile"] = U5;
	MolarFraction["Fertile"] = 1. - U5;

	return MolarFraction ;
	
}