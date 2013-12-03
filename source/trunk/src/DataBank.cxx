#include "DataBank.hxx"

#include "IsotopicVector.hxx"
#include "CLASSHeaders.hxx"
#include "EvolutionData.hxx"
#include "LogFile.hxx"
#include "StringLine.hxx"

#include <TGraph.h>
#include <TString.h>


#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>


//________________________________________________________________________
//
//		DataBank
//
//
//
//
//________________________________________________________________________

double ReactionRateWeightedDistance(IsotopicVector IV1, EvolutionData DB )
{
	
	double d2 = 0;
	double XS_total = 0;
	IsotopicVector IV2 = DB.GetIsotopicVectorAt(0.).GetActinidesComposition();
	IsotopicVector IVtmp = IV1 + IV2;
	map<ZAI ,double> IVtmpIsotopicQuantity = IVtmp.GetIsotopicQuantity();
	map<ZAI ,double >::iterator it;
	
	for( it = IVtmpIsotopicQuantity.begin(); it != IVtmpIsotopicQuantity.end(); it++)
	{
		double XS = 0;
		
		for(int i = 1; i < 4 ; i++)
			XS += DB.GetGetXSForAt(0., (*it).first, i);
		
		double Z1 = IV1.GetZAIIsotopicQuantity( (*it).first );
		double Z2 = IV2.GetZAIIsotopicQuantity( (*it).first );
		d2 += pow( (Z1-Z2)*XS , 2 );
		XS_total += (Z1+Z2)*XS/2;
	}
	
	
	return sqrt(d2)/XS_total;
}

double ReactionRateWeightedDistance(EvolutionData DB, IsotopicVector IV1  )
{
	return ReactionRateWeightedDistance( IV1, DB );
}


template<class T>
DataBank<T>::DataBank()
{
}

//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________

template<>
DataBank<ZAI>::DataBank(LogFile* Log, string DB_index_file, bool setlog, bool olfreadmethod)
{
	
	SetLog(Log);
	IsLog(setlog);
	fDataBaseIndex = DB_index_file;
	
	fOldReadMethod = olfreadmethod;
	fUseOldGeneration = false;
	
	// Warning
	if(PrintLog())
	{
		cout	<< "!!INFO!! !!!DataBank<ZAI>!!! A EvolutionData<ZAI> has been define :" << endl;
		cout	<< "\t His index is : \"" << DB_index_file << "\"" << endl << endl;
		
		GetLog()->fLog 	<< "!!INFO!! !!!DataBank<ZAI>!!! A EvolutionData<ZAI> has been define :" << endl;
		GetLog()->fLog	<< "\t His index is : \"" << DB_index_file << "\"" << endl << endl;
	}
	
}

//________________________________________________________________________
template<>
DataBank<ZAI>::~DataBank()
{
	
}

//________________________________________________________________________

template<>
IsotopicVector	DataBank<ZAI>::Evolution(const ZAI& zai, double dt)
{
	
	IsotopicVector	returnIV;
	
	map<ZAI ,EvolutionData >::iterator it = fDataBank.find(zai);
	
	if (it == fDataBank.end() )
	{
		ifstream DB_index(fDataBaseIndex.c_str());
		if( !DB_index)
		{
			cout << "!!!EVOLUTIVE DB!!!! Can't open \"" << fDataBaseIndex << "\"" << endl;
			GetLog()->fLog << "!!!EVOLUTIVE DB!!!! Can't open \"" << fDataBaseIndex << "\"" << endl;
			exit (1);
		}
		bool zaifind = false;
		string tmp;
		getline(DB_index,tmp);							// Read first line
		while (!DB_index.eof())
		{
			string line;
			int start=0;
			getline(DB_index,line);							// Read first line
			string first=StringLine::NextWord(line,start);				// Read first word
			
			if(first.size()==0) break;						// If First word is null.... quit
			
			int rZ=atoi(first.c_str());						// Get Z
			int rA=atoi(StringLine::NextWord(line,start).c_str());			// Get A
			int rI=atoi(StringLine::NextWord(line,start).c_str());			// Get Isomeric State
			
			if(rZ == zai.Z() && rA == zai.A() && rI == zai.I() )
			{
				string file_name = StringLine::NextWord(line,start);
				EvolutionData evolutionproduct = EvolutionData(GetLog(), file_name);
#pragma omp critical(DBupdate)
				{fDataBank.insert( pair<ZAI ,EvolutionData >(zai, evolutionproduct) );}
				returnIV = evolutionproduct.GetIsotopicVectorAt(dt);
				zaifind = true;
			}
		}
		
		if(zaifind == false)
		{
			GetLog()->fLog << "!!Warning!! !!!EVOLUTIVE DB!!! Oups... Can't Find the ZAI : " ;
			GetLog()->fLog << zai.Z() << " " << zai.A() << " "	<< zai.I() << "!!! It will be considered as stable !!" << endl;
			
			EvolutionData evolutionproduct = EvolutionData(GetLog()," " , false, zai);
			{fDataBank.insert( pair<ZAI, EvolutionData >(zai, evolutionproduct) );}
			returnIV = evolutionproduct.GetIsotopicVectorAt(dt);
			
			
		}
		
		
	}
	else	returnIV = (*it).second.GetIsotopicVectorAt(dt);
	
	return returnIV;
}

template<>
bool DataBank<ZAI>::IsDefine(const ZAI& zai) const
{
	
	map<ZAI ,EvolutionData > evolutiondb = (*this).GetDataBank();
	if (evolutiondb.find(zai) != evolutiondb.end())
		return true;
	else
		return false;
	
}
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________

template<>
void DataBank<ZAI>:: BuildEqns(double t, double *N, double *dNdt)
{}





//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
template<>
void DataBank<IsotopicVector>::ReadDataBase();

template<>
TMatrixT<double> DataBank<IsotopicVector>::GetFissionXsMatrix(EvolutionData EvolutionDataStep,double BUStep);

template<>
TMatrixT<double> DataBank<IsotopicVector>::GetCaptureXsMatrix(EvolutionData EvolutionDataStep,double BUStep);

template<>
TMatrixT<double> DataBank<IsotopicVector>::Getn2nXsMatrix(EvolutionData EvolutionDataStep,double BUStep);

template<>
TMatrixT<double> DataBank<IsotopicVector>::ExtractXS(EvolutionData EvolutionDataStep,double BUStep);

template<>
void DataBank<IsotopicVector>::BuildDecayMatrix();

template<>
EvolutionData DataBank<IsotopicVector>::OldGenerateEvolutionData(IsotopicVector isotopicvector, double cycletime, double Power);

template<>
string DataBank<IsotopicVector>::GetDecay(string DecayModes, double &BR,int &Iso, int &StartPos);

template<>
void DataBank<IsotopicVector>::SetTheMatrixToZero();

template<>
void DataBank<IsotopicVector>::ResetTheMatrix();

template<>
void DataBank<IsotopicVector>::SetTheNucleiVectorToZero();

template<>
void DataBank<IsotopicVector>::ResetTheNucleiVector();

template<>
void DataBank<IsotopicVector>::BuildEqns(double t, double *N, double *dNdt);

template<>
void DataBank<IsotopicVector>::SetTheMatrix(TMatrixT<double> BatemanMatrix);

template<>
void DataBank<IsotopicVector>::SetTheNucleiVector(TMatrixT<double> NEvolutionMatrix);

template<>
TMatrixT<double> DataBank<IsotopicVector>::GetTheNucleiVector();

//________________________________________________________________________
template<>
DataBank<IsotopicVector>::DataBank():DynamicalSystem()
{
	fTheNucleiVector = 0;
	fTheMatrix = 0;

	fOldReadMethod = true;
	fUseOldGeneration = false;
	fUseRK4EvolutionMethod = true;
	fDistanceType = 0;
	fShorstestHalflife = 3600.;
	
	fDataDirectoryName = getenv("CLASS_PATH");
	fDataDirectoryName += "/source/data/";
	fDataFileName = "chart.JEF3T";
	
	BuildDecayMatrix();
	
	fNVar = findex_inver.size();
	SetForbidNegativeValue();
	
}


template<>
DataBank<IsotopicVector>::DataBank(LogFile* Log, string DB_index_file, bool setlog, bool olfreadmethod):DynamicalSystem()
{
	
	SetLog(Log);
	IsLog(setlog);
	
	fTheNucleiVector = 0;
	fTheMatrix = 0;
	
	fDataDirectoryName = getenv("CLASS_PATH");
	fDataDirectoryName += "/source/data/";
	fDataFileName = "chart.JEF3T";
	
	fDataBaseIndex = DB_index_file;
	fOldReadMethod = olfreadmethod;
	fUseRK4EvolutionMethod = true;

	fUseOldGeneration = false;
	fDistanceType = 0;
	
	fShorstestHalflife = 3600.;
	
	BuildDecayMatrix();
	ReadDataBase();
	
	fNVar = findex_inver.size();
	SetForbidNegativeValue();
	
	if(PrintLog())
	{
		// Warning
		cout	<< "!!INFO!! !!!DataBank<IsotopicVector>!!! A EvolutionData<ZAI> has been define :" << endl;
		cout	<< "\t His index is : \"" << DB_index_file << "\"" << endl;
		cout	<< "\t " << fDataBank.size() << " EvolutionData have been read."<< endl << endl;
		
		GetLog()->fLog 	<< "!!INFO!! !!!DataBank<IsotopicVector>!!! A EvolutionData<ZAI> has been define :" << endl;
		GetLog()->fLog	<< "\t His index is : \"" << DB_index_file << "\"" << endl;
		GetLog()->fLog	<< "\t " << fDataBank.size() << " EvolutionData have been read."<< endl << endl;
	}
	
}

template<>
DataBank<IsotopicVector>::~DataBank()
{
}

//________________________________________________________________________

template<>
void DataBank<IsotopicVector>::SetFissionEnergy(ZAI zai, double E)
{
	pair<map<ZAI, double>::iterator, bool> IResult;
	IResult = fFissionEnergy.insert( pair<ZAI ,double>(zai, E));
	if(IResult.second == false)
		IResult.first->second = E;
	
}

template<>
void DataBank<IsotopicVector>::SetFissionEnergy(string FissionEnergyFile)
{
	ifstream FissionFile(FissionEnergyFile.c_str());	// Open the File
	if(!FissionFile)				//check if file is correctly open
	{
		cout << "!!Warning!! !!!DataBank!!! \n Can't open \"" << FissionFile << "\"\n" << endl;
		GetLog()->fLog << "!!Warning!! !!!DataBank!!! \n Can't open \"" << FissionFile << "\"\n" << endl;
	}
	
	do {
		int Z = 0;
		int A = 0;
		int I = 0;
		double E = 0;
		FissionFile >> Z >> A >> I >> E;
		SetFissionEnergy(Z, A, I, E);
	} while (!FissionFile.eof());
}

//________________________________________________________________________



template<>
void DataBank<IsotopicVector>::ReadDataBase()
{
	
	ifstream DataDB(fDataBaseIndex.c_str());							// Open the File
	if(!DataDB)
	{
		cout << "!!Warning!! !!!DataBank!!! \n Can't open \"" << fDataBaseIndex << "\"\n" << endl;
		GetLog()->fLog << "!!Warning!! !!!DataBank!!! \n Can't open \"" << fDataBaseIndex << "\"\n" << endl;
	}
	vector<double> vTime;
	vector<double> vTimeErr;
	
	string line;
	int start = 0;
	
	
	
	// First Get Fuel Type
	getline(DataDB, line);
	if( StringLine::NextWord(line, start, ' ') != "TYPE")
	{
		cout << "!!Bad Trouble!! !!!DataBank!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the type of the DataBase"<< endl;
		GetLog()->fLog << "!!Bad Trouble!! !!!DataBank!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the type of the DataBase"<< endl;
		exit (1);
	}
	fFuelType = StringLine::NextWord(line, start, ' ');
	// First Get Fuel Parameter
	getline(DataDB, line);
	start = 0;
	if( StringLine::NextWord(line, start, ' ') != "PARAM")
	{
		cout << "!!Bad Trouble!! !!!DataBank!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the Parameter of the DataBase"<< endl;
		GetLog()->fLog << "!!Bad Trouble!! !!!DataBank!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the Parameter of the DataBase"<< endl;
		exit (1);
	}
	while(start < (int)line.size())
		fFuelParameter.push_back(atof(StringLine::NextWord(line, start, ' ').c_str()));
	
	//Then Get All the Database
	
	while (!DataDB.eof())
	{
		getline(DataDB, line);
		if(line != "")
		{
			EvolutionData* evolutionproduct = new EvolutionData(GetLog(), line, fOldReadMethod);
			IsotopicVector ivtmp  = evolutionproduct->GetIsotopicVectorAt(0.).GetActinidesComposition();
			fDataBank.insert( pair<IsotopicVector, EvolutionData >(ivtmp , (*evolutionproduct) ));
		}
	}
	
}
//________________________________________________________________________



template<>
map<double, EvolutionData> DataBank<IsotopicVector>::GetDistancesTo(IsotopicVector isotopicvector, double t) const
{
	
	map<double, EvolutionData> distances;
	
	map<IsotopicVector, EvolutionData > evolutiondb = fDataBank;
	
	map<IsotopicVector, EvolutionData >::iterator it;
	for( it = evolutiondb.begin(); it != evolutiondb.end(); it++ )
	{
		pair<map<double, EvolutionData>::iterator, bool> IResult;
		
		double D = Distance(isotopicvector.GetActinidesComposition(), (*it).second.GetIsotopicVectorAt(t).GetActinidesComposition()/ Norme( (*it).second.GetIsotopicVectorAt(t).GetActinidesComposition() )*Norme(isotopicvector.GetActinidesComposition())
							,fDistanceType, fDistanceParameter);
		
		IResult = distances.insert( pair<double, EvolutionData>( D , (*it).second ) );
	}
	
	return distances;
	
}

template<>
EvolutionData DataBank<IsotopicVector>::GetClosest(IsotopicVector isotopicvector, double t) const
{
	
	map<IsotopicVector, EvolutionData > evolutiondb = fDataBank;
	
	double distance = Distance(isotopicvector.GetActinidesComposition(),
							   evolutiondb.begin()->second.GetIsotopicVectorAt(t).GetActinidesComposition()
							   / evolutiondb.begin()->second.GetIsotopicVectorAt(t).GetActinidesComposition().GetSumOfAll()
							   * isotopicvector.GetActinidesComposition().GetSumOfAll(),
							   fDistanceType, fDistanceParameter);
	
	EvolutionData CloseEvolData = evolutiondb.begin()->second ;
	
	map<IsotopicVector, EvolutionData >::iterator it;
	for( it = evolutiondb.begin(); it != evolutiondb.end(); it++ )
	{
		pair<map<double, EvolutionData>::iterator, bool> IResult;
		double D = Distance(isotopicvector.GetActinidesComposition(),
							(*it).second.GetIsotopicVectorAt(t).GetActinidesComposition()
							/  (*it).second.GetIsotopicVectorAt(t).GetActinidesComposition().GetSumOfAll()
							* isotopicvector.GetActinidesComposition().GetSumOfAll(),
							fDistanceType, fDistanceParameter);
		if (D< distance)
		{
			distance = D;
			CloseEvolData = (*it).second;
		}
	}
	return CloseEvolData;
	
}





//________________________________________________________________________
//________________________________________________________________________

template<>
void DataBank<IsotopicVector>::CalculateDistanceParameter()
{
	
	if(fDistanceType!=1){
		cout << "!!Warning!! !!!CalculateDistanceParameter!!!"
		<< " Distance Parameter will be calculate even if the distance type is not the good one. Any Distance Parameters given by the user will be overwriten"<<endl;
		
		GetLog()->fLog << "!!Warning!! !!!CalculateDistanceParameter!!!"
		<< " Distance Parameter will be calculate even if the distance type is not the good one. Any Distance Parameters given by the user will be overwriten"<<endl;
	}
	
	fDistanceParameter.Clear();
	
	//We calculate the weight for the distance calculation.
	map<IsotopicVector ,EvolutionData >::iterator it;
	map<IsotopicVector ,EvolutionData > databank = (*this).GetDataBank();
	int NevolutionDatainDataBank=0;
	
	for( it = databank.begin(); it != databank.end(); it++ ){
		NevolutionDatainDataBank++;
		map<ZAI ,double>::iterator itit;
		map<ZAI ,double> isovector=(*it).first.GetIsotopicQuantity();
		for(itit=isovector.begin(); itit != isovector.end(); itit++){//Boucle sur ZAI
			ZAI TmpZAI=(*itit).first;
			double TmpXS=0;
			for(int i=1;i<4;i++){		//Loop on Reactions 1==fission, 2==capture, 3==n2n
				TmpXS+=	(*it).second.GetGetXSForAt(0,TmpZAI,i);
			}
			fDistanceParameter.Add(TmpZAI,TmpXS);
		}
		
		
	}
	fDistanceParameter.Multiply((double)1.0/NevolutionDatainDataBank);
	
	
	GetLog()->fLog <<"!!INFO!! Distance Parameters "<<endl;
	map<ZAI ,double >::iterator it2;
	for(it2 = fDistanceParameter.GetIsotopicQuantity().begin();it2 != fDistanceParameter.GetIsotopicQuantity().end(); it2++)
	{
		GetLog()->fLog << (*it2).first.Z() << " ";
		GetLog()->fLog << (*it2).first.A() << " ";
		GetLog()->fLog << (*it2).first.I() << " ";
		GetLog()->fLog << ": " << (*it2).second;
		GetLog()->fLog << endl;
	}
	GetLog()->fLog << endl;
	
	
	
}

//________________________________________________________________________
template<>
void DataBank<IsotopicVector>::SetDistanceParameter(IsotopicVector DistanceParameter){
	
	fDistanceParameter=DistanceParameter;
	
	GetLog()->fLog <<"!!INFO!! Distance Parameters "<<endl;
	map<ZAI ,double >::iterator it2;
	for(it2 = fDistanceParameter.GetIsotopicQuantity().begin();it2 != fDistanceParameter.GetIsotopicQuantity().end(); it2++)
	{
		GetLog()->fLog << (*it2).first.Z() << " ";
		GetLog()->fLog << (*it2).first.A() << " ";
		GetLog()->fLog << (*it2).first.I() << " ";
		GetLog()->fLog << ": " << (*it2).second;
		GetLog()->fLog << endl;
	}
	GetLog()->fLog << endl;
	
}

//________________________________________________________________________
template<>
void DataBank<IsotopicVector>::SetDistanceType(int DistanceType)
{
	
	fDistanceType=DistanceType;
	if(fDistanceType==1){
		CalculateDistanceParameter();
	}
	else if(fDistanceType == 2 && Norme(fDistanceParameter)==0){
		// This is so bad!! You will probably unsynchronize all the reactor....
		cout << "!!Warning!! !!!DistanceType!!!"
		<< " Distance use weight defined by user for each isotope, but no weight have been given" << endl<<"Use SetDistanceParameter()"<<endl;
		
		GetLog()->fLog << "!!Warning!! !!!DistanceType!!!"
		<< " Distance use weight defined by user for each isotope, but no weight have been given" << endl<<"Use SetDistanceParameter()"<<endl;
		exit(1);
	}
	else if (fDistanceType != 0 && fDistanceType != 1 && fDistanceType != 2 ){
		cout << "!!ERROR!! !!!DistanceType!!!"
		<< " Distancetype defined by the user isn't recognized by the code"<<endl;
		
		GetLog()->fLog << "!!ERROR!! !!!DistanceType!!!"
		<< " Distancetype defined by the user isn't recognized by the code"<<endl;
		exit(1);
	}
	
}



template<>
EvolutionData DataBank<IsotopicVector>::GenerateEvolutionData(IsotopicVector isotopicvector, double cycletime, double Power)
{
	SetTheMatrixToZero();
	SetTheNucleiVectorToZero();
	
	if(fUseOldGeneration)
	{
		EvolutionData GeneratedDB = OldGenerateEvolutionData( isotopicvector,  cycletime, Power);
		fDataBankCalculated.insert( pair< IsotopicVector, EvolutionData > ( GeneratedDB.GetIsotopicVectorAt(0.), GeneratedDB) );
		return GeneratedDB;
		
	}
	string ReactorType;
	
	
	vector< TMatrixT<double> > NMatrix ;//  TMatrixT<double>(decayindex.size(),1))
	{	// Filling the t=0 State;
		map<ZAI, double > isotopicquantity = isotopicvector.GetIsotopicQuantity();
		TMatrixT<double>  N_0Matrix =  TMatrixT<double>( findex.size(),1) ;
		
		map<ZAI, double >::iterator it ;
		for(int i = 0; i < (int)findex.size(); i++)
			N_0Matrix[i] = 0;
		
		for(it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		{
			
			map<ZAI, int >::iterator it2;
			
			if( (*it).first.Z() < 90 )
				it2 = findex_inver.find( ZAI(-2,-2,-2) );
			else it2 = findex_inver.find( (*it).first );
			
			if(it2 == findex_inver.end() )		//If not in index should be TMP, can't be fast decay for new Fuel !!!
				it2 = findex_inver.find( ZAI(-3,-3,-3) );
			N_0Matrix[ (*it2).second ][0] = (*it).second ;


		}
		NMatrix.push_back(N_0Matrix);
		
	}
	

	//-------------------------//
	//--- Perform Evolution ---//
	//-------------------------//
	EvolutionData EvolutionDataStep = GetClosest(isotopicvector.GetActinidesComposition(), 0.);	//GetCLosest at the begining of evolution
	ReactorType = EvolutionDataStep.GetReactorType();
	
	double Na = 6.02214129e23;	//N Avogadro
	double M_ref = 0;
	double M = 0;
	double Power_ref =  EvolutionDataStep.GetPower();
	{
		map<ZAI, double >::iterator it ;
		
		
		IsotopicVector IVtmp = isotopicvector.GetActinidesComposition() + EvolutionDataStep.GetIsotopicVectorAt(0.).GetActinidesComposition();
		map<ZAI, double >isotopicquantity = IVtmp.GetIsotopicQuantity();
		
		for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
		{
			M_ref += EvolutionDataStep.GetIsotopicVectorAt(0.).GetActinidesComposition().GetZAIIsotopicQuantity( (*it).first )*cZAIMass.fZAIMass.find( (*it).first )->second/Na*1e-6;
			M += isotopicvector.GetActinidesComposition().GetZAIIsotopicQuantity( (*it).first )*cZAIMass.fZAIMass.find( (*it).first )->second/Na*1e-6;
		}
	}

	int DBTimeStepN = EvolutionDataStep.GetFissionXS().begin()->second->GetN();
	double* DBTimeStep = EvolutionDataStep.GetFissionXS().begin()->second->GetX();
	
	int InsideStep = 10;
	
	int NStep = (DBTimeStepN);
	double timevector[NStep];
	timevector[0] = 0;

	double  Flux[NStep];
	
	TMatrixT<double> SigmaPhi = TMatrixT<double>(findex.size()*3+1,NStep); // Store the XS and the flux trought the evolution calculation.
	
	TMatrixT<double> FissionEnergy = TMatrixT<double>(findex.size(),1);
	{
		map< ZAI, int >::iterator it;
		for(it = findex_inver.begin(); it != findex_inver.end(); it++)
		{
			map< ZAI, double >::iterator it2 = fFissionEnergy.find(it->first);
			if(it2 == fFissionEnergy.end())
			{
				if(it->first.Z() > 90)
					FissionEnergy[it->second][0] = 1.9679e6*it->first.A()-2.601e8; // //simple linear fit to known values ;extrapolation to unknown isotopes
				else FissionEnergy[it->second][0] = 0;
			}
			else
					FissionEnergy[it->second][0] = it2->second;

		}
	}
	
	vector< TMatrixT<double> > FissionXSMatrix; //The Fisison XS Matrix
	vector< TMatrixT<double> > CaptureXSMatrix; //The Capture XS Matrix
	vector< TMatrixT<double> > n2nXSMatrix;	 //The n2N XS Matrix

	for(int i = 0; i < NStep-1; i++)
	{
		double TStepMax = ( (DBTimeStep[i+1]-DBTimeStep[i] ) ) * Power_ref/M_ref / Power*M ;

		
		TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
		TMatrixT<double> BatemanReactionMatrix = TMatrixT<double>(findex.size(),findex.size());

		TMatrixT<double> NEvolutionMatrix = TMatrixT<double>(findex.size(),1);
		NEvolutionMatrix = NMatrix.back();
		
		
		
		FissionXSMatrix.push_back(GetFissionXsMatrix(EvolutionDataStep, DBTimeStep[i])); //Feel the reaction Matrix
		CaptureXSMatrix.push_back(GetCaptureXsMatrix(EvolutionDataStep, DBTimeStep[i])); //Feel the reaction Matrix
		n2nXSMatrix.push_back(Getn2nXsMatrix(EvolutionDataStep, DBTimeStep[i])); //Feel the reaction Matrix

		// ----------------   Evolution

		BatemanReactionMatrix = FissionXSMatrix[i];
		BatemanReactionMatrix += CaptureXSMatrix[i];
		BatemanReactionMatrix += n2nXSMatrix[i];

		if(fUseRK4EvolutionMethod)
		{
			for(int k=0; k < InsideStep; k++)
			{
				double ESigmaN = 0;
				for (int j = 0; j < (int)findex.size() ; j++)
					ESigmaN -= FissionXSMatrix[i][j][j]*NEvolutionMatrix[j][0]*1.6e-19*FissionEnergy[j][0];
				// Update Flux
				double Flux_k = Power/ESigmaN;

				if(k==0)
					Flux[i]=Flux_k;

				BatemanMatrix = BatemanReactionMatrix;
				BatemanMatrix *= Flux_k;
				BatemanMatrix += fDecayMatrix ;

				SetTheMatrixToZero();
				SetTheNucleiVectorToZero();

				SetTheMatrix(BatemanMatrix);
				SetTheNucleiVector(NEvolutionMatrix);


				RungeKutta(fTheNucleiVector, timevector[i]+TStepMax/InsideStep*k, timevector[i]+TStepMax/InsideStep*(k+1),  fNVar);
				NEvolutionMatrix = GetTheNucleiVector();

			}
			NEvolutionMatrix = GetTheNucleiVector();
			NMatrix.push_back(NEvolutionMatrix);
		}
		else
		{

			for(int k=0; k < InsideStep; k++)
			{
				double ESigmaN = 0;
				for (int j = 0; j < (int)findex.size() ; j++)
					ESigmaN -= FissionXSMatrix[i][j][j]*NEvolutionMatrix[j][0]*1.6e-19*FissionEnergy[j][0];
				// Update Flux
				double Flux_k = Power/ESigmaN;

				if(k==0)
					Flux[i]=Flux_k;

				BatemanMatrix = BatemanReactionMatrix;
				BatemanMatrix *= Flux_k;
				BatemanMatrix += fDecayMatrix ;
				BatemanMatrix *= TStepMax/InsideStep ;


				TMatrixT<double> IdMatrix = TMatrixT<double>(findex.size(),findex.size());
				for(int j = 0; j < (int)findex.size(); j++)
					for(int k = 0; k < (int)findex.size(); k++)
					{
						if(k == j)	IdMatrix[j][k] = 1;
						else 		IdMatrix[j][k] = 0;
					}


				TMatrixT<double> BatemanMatrixDL = TMatrixT<double>(findex.size(),findex.size());   // Order 0 Term from the DL : Id
				TMatrixT<double> BatemanMatrixDLTermN = TMatrixT<double>(findex.size(),findex.size());  // Addind it;

				{
					BatemanMatrix *= TStepMax ;
					BatemanMatrixDLTermN = IdMatrix;
					BatemanMatrixDL = BatemanMatrixDLTermN;
					int j = 1;
					double NormN;

					do
					{
						TMatrixT<double> BatemanMatrixDLTermtmp = TMatrixT<double>(findex.size(),findex.size());  // Adding it;
						BatemanMatrixDLTermtmp = BatemanMatrixDLTermN;

						BatemanMatrixDLTermN.Mult(BatemanMatrixDLTermtmp, BatemanMatrix );

						BatemanMatrixDLTermN *= 1./j;
						BatemanMatrixDL += BatemanMatrixDLTermN;

						NormN = 0;
						for(int m = 0; m < (int)findex.size(); m++)
							for(int n = 0; n < (int)findex.size(); n++)
								NormN += BatemanMatrixDLTermN[m][n]*BatemanMatrixDLTermN[m][n];
						j++;
						
					} while ( NormN != 0 );
				}
				NEvolutionMatrix = BatemanMatrixDL * NEvolutionMatrix ;
			}
			NMatrix.push_back(NEvolutionMatrix);
			
		}

		timevector[i+1] = timevector[i] + TStepMax;


	}
	FissionXSMatrix.push_back(GetFissionXsMatrix(EvolutionDataStep, DBTimeStep[NStep])); //Feel the reaction Matrix
	CaptureXSMatrix.push_back(GetCaptureXsMatrix(EvolutionDataStep, DBTimeStep[NStep])); //Feel the reaction Matrix
	n2nXSMatrix.push_back(Getn2nXsMatrix(EvolutionDataStep, DBTimeStep[NStep])); //Feel the reaction Matrix


	EvolutionData GeneratedDB = EvolutionData(GetLog());


	double ESigmaN = 0;
	for (int j = 0; j < (int)findex.size() ; j++)
		ESigmaN -= FissionXSMatrix[NStep-1][j][j]*NMatrix.back()[j][0]*1.6e-19*FissionEnergy[j][0];
	Flux[NStep-1] = Power/ESigmaN;

	GeneratedDB.SetFlux( new TGraph(NStep, timevector, Flux)  );
	
	for(int i = 0; i < (int)findex.size(); i++)
	{
		double ZAIQuantity[NMatrix.size()];
		double FissionXS[NStep];
		double CaptureXS[NStep];
		double n2nXS[NStep];
		for(int j = 0; j < (int)NMatrix.size(); j++)
			ZAIQuantity[j] = (NMatrix[j])[i][0];

		for(int j = 0; j < NStep; j++)
		{
			FissionXS[j]	= FissionXSMatrix[j][i][i];
			CaptureXS[j]	= CaptureXSMatrix[j][i][i];
			n2nXS[j]	= n2nXSMatrix[j][i][i];
		}
		
		GeneratedDB.NucleiInsert(pair<ZAI, TGraph*> (findex.find(i)->second, new TGraph(NMatrix.size(), timevector, ZAIQuantity)));
		GeneratedDB.FissionXSInsert(pair<ZAI, TGraph*> (findex.find(i)->second, new TGraph(NStep, timevector, FissionXS)));
		GeneratedDB.CaptureXSInsert(pair<ZAI, TGraph*> (findex.find(i)->second, new TGraph(NStep, timevector, CaptureXS)));
		GeneratedDB.n2nXSInsert(pair<ZAI, TGraph*> (findex.find(i)->second, new TGraph(NStep, timevector, n2nXS)));
		}
	
	GeneratedDB.SetPower(Power );
	GeneratedDB.SetFuelType(fFuelType );
	GeneratedDB.SetReactorType(ReactorType );
	GeneratedDB.SetCycleTime(cycletime);
	
	fDataBankCalculated.insert( pair< IsotopicVector, EvolutionData > ( GeneratedDB.GetIsotopicVectorAt(0.), GeneratedDB) );

	//exit(1);
	ResetTheMatrix();
	ResetTheNucleiVector();
	return GeneratedDB;
	
}


//________________________________________________________________________

template<>
TMatrixT<double> DataBank<IsotopicVector>::GetFissionXsMatrix(EvolutionData EvolutionDataStep,double TStep)
{
	
	map<ZAI ,TGraph* >::iterator it;
	TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
	
	// ----------------  A(n,.) X+Y
	
	map<ZAI ,TGraph* > FissionXS = EvolutionDataStep.GetFissionXS();
	
	for(it = FissionXS.begin() ; it != FissionXS.end(); it++)
	{
		map<ZAI, int>::iterator findex_inver_it = findex_inver.find( (*it).first );
		if( findex_inver_it != findex_inver.end() )
		{
			double y = (*it).second->Eval(TStep);
			BatemanMatrix[ findex_inver_it->second ][ findex_inver_it->second ] += -y* 1e-24;
			BatemanMatrix[1][ findex_inver_it->second ] += 2*y* 1e-24;
		}
		
	}
	
	return BatemanMatrix;
	
}
//________________________________________________________________________
template<>
TMatrixT<double> DataBank<IsotopicVector>::GetCaptureXsMatrix(EvolutionData EvolutionDataStep,double TStep)
{
	
	
	map<ZAI ,TGraph* >::iterator it;
	TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
	
	map<ZAI, map<ZAI, double> > Capture;
	{	// 241Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(96,242,0) , 0.8733*0.827) ); //directly cut the Am242 as in MURE
		toAdd.insert(pair<ZAI, double> ( ZAI(94,242,0) , 0.8733*0.173) ); //directly cut the Am242 as in MURE
		toAdd.insert(pair<ZAI, double> ( ZAI(95,242,1) , 0.1267) );
		Capture.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,241,0), toAdd ) );
	}
	{	// 242Am*
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(95,243,0) , 1) );
		Capture.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,242,1), toAdd ) );
	}
	
	
	// ----------------  A(n,.)A+1
	map<ZAI ,TGraph* > CaptureXS = EvolutionDataStep.GetCaptureXS();
	for(it = CaptureXS.begin(); it != CaptureXS.end(); it++)
	{
		map<ZAI, int>::iterator Index_it = findex_inver.find( (*it).first );
		if( Index_it != findex_inver.end() )
		{
			double y;
			y = (*it).second->Eval(TStep);
			
			BatemanMatrix[Index_it->second][ Index_it->second ] += -y* 1e-24 ;
			
			map<ZAI, map<ZAI, double> >::iterator it3 = Capture.find( (*it).first );
			
			if( it3 == Capture.end() )
			{
				map<ZAI, int >::iterator it6 = findex_inver.find( ZAI( (*it).first.Z(), (*it).first.A()+1, (*it).first.I()) );
				
				if( it6 != findex_inver.end() )
				{
					BatemanMatrix[(*it6).second][Index_it->second] += y* 1e-24  ;
				}
				else
				{
					map<ZAI, map<ZAI, double> >::iterator it4 = fFastDecay.find(  ZAI( (*it).first.Z(), (*it).first.A()+1, (*it).first.I()) );
					
					if( it4 == fFastDecay.end() )
					{
						cout << "Capture Problem in FastDecay for nuclei " << (*it).first.Z() << " " << (*it).first.A()+1 << " " << (*it).first.I() << endl;
						
						BatemanMatrix[0][Index_it->second] += y* 1e-24  ;
					}
					else
					{
						
						map< ZAI, double >::iterator it5;
						map< ZAI, double > decaylist2 = (*it4).second;
						for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
						{
							it6 = findex_inver.find( (*it5).first );
							if( it6 == findex_inver.end() )
								BatemanMatrix[0][Index_it->second] += y* 1e-24 * (*it5).second;
							else
								BatemanMatrix[(*it6).second][Index_it->second] += y* 1e-24 * (*it5).second;
						}
					}
				}
			}
			else
			{
				map<ZAI, double>::iterator it4;
				map<ZAI, double> CaptureList = (*it3).second;
				for(it4 = CaptureList.begin(); it4 != CaptureList.end() ; it4++)
				{
					
					map<ZAI, int >::iterator it6 = findex_inver.find( (*it4).first );
					if( it6 != findex_inver.end() )
						BatemanMatrix[(*it6).second][Index_it->second] += y* 1e-24 * (*it4).second ;
					else
					{
						map<ZAI, map<ZAI, double> >::iterator it7 = fFastDecay.find( (*it4).first );
						
						if( it7 == fFastDecay.end() )
						{
							cout << "CaptureList Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
							exit(1);
						}
						
						map< ZAI, double >::iterator it5;
						map< ZAI, double > decaylist2 = (*it7).second;
						for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
						{
							
							it6 = findex_inver.find( (*it5).first );
							if( it6 == findex_inver.end() )
							{
								cout << "CaptureList Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
								exit(1);
							}
							
							BatemanMatrix[(*it6).second][Index_it->second] += y * 1e-24 * (*it5).second * (*it4).second;
						}
					}
					
				}
			}
			
			
		}
	}
	return BatemanMatrix;
	
}


//________________________________________________________________________

template<>
TMatrixT<double> DataBank<IsotopicVector>::Getn2nXsMatrix(EvolutionData EvolutionDataStep,double TStep)
{
	
	
	map<ZAI ,TGraph* >::iterator it;
	TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
	
	map<ZAI, map<ZAI, double> > n2n;
	{	// 237Np
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(93,236,0) , 0.2) );
		toAdd.insert(pair<ZAI, double> ( ZAI(93,236,1) , 0.8) );
		n2n.insert( pair< ZAI, map<ZAI, double> > ( ZAI(93,237,0), toAdd ) );
	}
	{	// 242Am*
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(95,241,0) , 1) );
		n2n.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,242,1), toAdd ) );
	}
	
	// ----------------  A(n,2n)A-1
	map<ZAI ,TGraph* > n2nXS = EvolutionDataStep.Getn2nXS();
	for(it = n2nXS.begin(); it != n2nXS.end(); it++)
	{
		map<ZAI, int>::iterator Index_it = findex_inver.find( (*it).first );
		if( Index_it != findex_inver.end() )
		{
			double y;
			y = (*it).second->Eval(TStep);
			
			BatemanMatrix[Index_it->second][ Index_it->second ] += -y* 1e-24 ;
			
			map<ZAI, map<ZAI, double> >::iterator it3 = n2n.find( (*it).first );
			
			if( it3 == n2n.end() )
			{
				map<ZAI, int >::iterator it6 = findex_inver.find( ZAI( (*it).first.Z(), (*it).first.A()-1, (*it).first.I()) );
				
				if( it6 != findex_inver.end() )
				{
					BatemanMatrix[(*it6).second][Index_it->second] += y* 1e-24  ;
				}
				else
				{
					map<ZAI, map<ZAI, double> >::iterator it4 = fFastDecay.find(  ZAI( (*it).first.Z(), (*it).first.A()-1, (*it).first.I()) );
					
					if( it4 == fFastDecay.end() )
					{
						cout << "n2n Problem in FastDecay for nuclei " << (*it).first.Z() << " " << (*it).first.A()-1 << " " << (*it).first.I() << endl;
						
						BatemanMatrix[0][Index_it->second] += y* 1e-24  ;
					}
					else
					{
						
						map< ZAI, double >::iterator it5;
						map< ZAI, double > decaylist2 = (*it4).second;
						for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
						{
							it6 = findex_inver.find( (*it5).first );
							if( it6 == findex_inver.end() )
								BatemanMatrix[0][Index_it->second] += y* 1e-24 * (*it5).second;
							else
								BatemanMatrix[(*it6).second][Index_it->second] += y* 1e-24 * (*it5).second;
						}
					}
				}
			}
			else
			{
				map<ZAI, double>::iterator it4;
				map<ZAI, double> n2nList = (*it3).second;
				for(it4 = n2nList.begin(); it4 != n2nList.end() ; it4++)
				{
					
					map<ZAI, int >::iterator it6 = findex_inver.find( (*it4).first );
					if( it6 != findex_inver.end() )
						BatemanMatrix[(*it6).second][Index_it->second] += y* 1e-24 * (*it4).second ;
					else
					{
						map<ZAI, map<ZAI, double> >::iterator it7 = fFastDecay.find( (*it4).first );
						
						if( it7 == fFastDecay.end() )
						{
							cout << "n2nList Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
							exit(1);
						}
						
						map< ZAI, double >::iterator it5;
						map< ZAI, double > decaylist2 = (*it7).second;
						for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
						{
							
							it6 = findex_inver.find( (*it5).first );
							if( it6 == findex_inver.end() )
							{
								cout << "n2nList Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
								exit(1);
							}
							
							BatemanMatrix[(*it6).second][Index_it->second] += y * 1e-24 * (*it5).second * (*it4).second;
						}
					}
					
				}
			}
			
			
		}
	}
	return BatemanMatrix;
}


template<>
TMatrixT<double> DataBank<IsotopicVector>::ExtractXS(EvolutionData EvolutionDataStep,double TStep)
{
	
	
	map<ZAI ,TGraph* >::iterator it;
	// ----------------  A(n,.) X+Y
	
	map<ZAI ,TGraph* > FissionXS = EvolutionDataStep.GetFissionXS();
	
	TMatrixT<double> SigmaPhi = TMatrixT<double>(findex.size()*3+1,1);
	for(it = FissionXS.begin() ; it != FissionXS.end(); it++)
	{
		
		if( findex_inver.find( (*it).first ) != findex_inver.end() )
		{
			double y;
			y = (*it).second->Eval(TStep);
			SigmaPhi[findex_inver.find( (*it).first )->second][0] = y ;
		}
		
	}
	
	// ----------------  A(n,.)A+1
	map<ZAI ,TGraph* > CaptureXS = EvolutionDataStep.GetCaptureXS();
	for(it = CaptureXS.begin(); it != CaptureXS.end(); it++)
	{
		if( findex_inver.find( (*it).first ) != findex_inver.end() )
		{
			double y;
			y = (*it).second->Eval(TStep);
			SigmaPhi[findex_inver.find( (*it).first )->second + findex.size() ][0] = y ;
			
		}
	}
	
	// ----------------  A(n,2n)A-1
	map<ZAI ,TGraph* > n2nXS = EvolutionDataStep.Getn2nXS();
	for(it = n2nXS.begin() ; it != n2nXS.end(); it++)
	{
		if( findex_inver.find( (*it).first ) != findex_inver.end() )
		{
			double y;
			y = (*it).second->Eval(TStep);
			SigmaPhi[findex_inver.find( (*it).first )->second + findex.size() + findex.size()][0] = y ;
			
		}
	}
	return SigmaPhi;
}
//________________________________________________________________________
//________________________________________________________________________
template<>
string DataBank<IsotopicVector>::GetDecay(string DecayModes, double &BR,int &Iso, int &StartPos)
{
	string header;
	
	BR=0;
	//extraction of the decay mode and the BR
	string DecayBR=StringLine::NextWord(DecayModes,StartPos,',');
	//extraction of the decay
	int ss=0;
	string Decay=StringLine::NextWord(DecayBR,ss,':');
	//extraction of the BR if exist (i.e. for non stable isotop)
	if(ss<int(DecayBR.size()))
		BR=atof(DecayBR.substr(ss+1).c_str());
	//BR in % -> BR
	BR/=100.;
	//find the Isomeric state of Daughter
	Iso=0;
	if(Decay.find("/",0)<string::npos)
	{
		Iso=atoi(Decay.substr(Decay.find("/")+1).c_str());
		Decay=Decay.substr(0,Decay.find("/"));
	}
	return Decay;
}



template<>
void DataBank<IsotopicVector>::BuildDecayMatrix()
{
	// List of Decay Time and Properties
	map<ZAI, pair<double, map< ZAI, double > > > ZAIDecay;
	{	// TMP
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(-3,-3,-3), pair<double, map< ZAI, double > > ( 1e28 ,toAdd )) ) ;
	}
	{	// PF
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-2,-2,-2), 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(-2,-2,-2), pair<double, map< ZAI, double > > ( 1e28 ,toAdd )) ) ;
	}
	
	string DataFullPathName = GetDataDirectoryName()+ GetDataFileName();
	ifstream infile(DataFullPathName.c_str());
	
	if(!infile)
	{
		cout << "!!Warning!! !!!DataBank!!! \n Can't open \"" << DataFullPathName << "\"\n" << endl;
		GetLog()->fLog << "!!Warning!! !!!DataBank!!! \n Can't open \"" << DataFullPathName<< "\"\n" << endl;
	}
	
	
	
	
	
	do
	{
		int A = -1;
		int Z = -1;
		int I = -1;
		string zainame;
		string Iname;
		int unkown;
		double HalfLife;
		string DecayModes;
		
		infile >> A >> Z >> zainame >> Iname >> unkown >> HalfLife >> DecayModes;
		if(Z >= 90 )
		{
			// Get Isomeric State;
			
			if(Iname=="gs")
				I = 0;
			else
				if(Iname[0]=='m')
				{
					if( atoi( Iname.substr(1).c_str() )==0 )
						I = 1;
					else
						I = atoi( Iname.substr(1).c_str() );
					
				}
			
			
			int start=0;
			double branch_test=0;
			double branch_test_f=0;
			
			ZAI ParentZAI = ZAI(Z,A,I);
			map< ZAI, double > DaughtersMap;
			bool stable = true;

			while(start<int(DecayModes.size()))
			{
				ZAI DaughterZAI;
				double BR;
				int daughter_A=0;
				int daughter_N=0;
				int daughter_Z=0;
				int Iso=0;
				int DM=-1;
				//				FPDistribution *FP=0;
				string decay_name = GetDecay(DecayModes, BR, Iso, start);
				
				if (decay_name == "s")	{DM=0;}
				if (decay_name == "b-")	{DM=1;stable=false;	daughter_N=(A-Z)-1;	daughter_Z=Z+1;}
				if (decay_name == "n")	{DM=2;stable=false;	daughter_N=(A-Z)-1;	daughter_Z=Z;}
				if (decay_name == "nn")	{DM=3;stable=false;	daughter_N=(A-Z)-2;	daughter_Z=Z;}
				if (decay_name == "b-n"){DM=4;stable=false;	daughter_N=(A-Z)-2;	daughter_Z=Z+1;}
				if (decay_name == "p")	{DM=5;stable=false;	daughter_N=(A-Z);	daughter_Z=Z-1;}
				if (decay_name == "b-a"){DM=6;stable=false;	daughter_N=(A-Z)-3;	daughter_Z=Z-1;}
				if (decay_name == "pp")	{DM=7;stable=false;	daughter_N=(A-Z);	daughter_Z=Z-2;}
				if (decay_name == "ce")	{DM=8;stable=false;	daughter_N=(A-Z)+1;	daughter_Z=Z-1;}
				if (decay_name == "a")	{DM=9;stable=false;	daughter_N=(A-Z)-2;	daughter_Z=Z-2;}
				if (decay_name == "cen"){DM=10;stable=false;	daughter_N=(A-Z);	daughter_Z=Z-1;}
				if (decay_name == "cep"){DM=11;stable=false;	daughter_N=(A-Z)+1;	daughter_Z=Z-2;}
				if (decay_name == "it")	{DM=12;stable=false;	daughter_N=(A-Z);	daughter_Z=Z;Iso = I-1;}
				if (decay_name == "db-"){DM=13;stable=false;	daughter_N=(A-Z)-2;	daughter_Z=Z+2;}
				if (decay_name == "db+"){DM=14;stable=false;	daughter_N=(A-Z)+2;	daughter_Z=Z-2;}
				if (decay_name == "ita"){DM=15;stable=false;	daughter_N=(A-Z)-2;	daughter_Z=Z-2;}
				if (decay_name == "sf")	{DM=16;stable=false;	daughter_N=0;		daughter_Z=-2;	Iso = -2;}
				if (decay_name == "cesf"){DM=17;stable=false;	daughter_N=0;		daughter_Z=-2;	Iso = -2;}
				if (decay_name == "b-sf"){DM=18;stable=false;	daughter_N=0;		daughter_Z=-2;	Iso = -2;}

				daughter_A = daughter_Z + daughter_N;
				{
					if( daughter_Z < 90 && daughter_Z!=-2 )
						daughter_A = daughter_Z = Iso = -3;
					// not spontaneous fission
					ZAI DaughterZAI = ZAI(daughter_Z,daughter_A,Iso);
					if((BR>1e-10) && (!stable))
					{
						if(DM <= 15)
						{
							pair<map<ZAI, double>::iterator, bool> IResult;
							IResult = DaughtersMap.insert(pair<ZAI, double> (DaughterZAI , BR) );
							if( !IResult.second)
								(*IResult.first).second += BR;

							branch_test+=BR;
						}
						else if( DM <= 18)
						{
							pair<map<ZAI, double>::iterator, bool> IResult;
							IResult = DaughtersMap.insert(pair<ZAI, double> (ZAI(-2,-2,-2) , 2*BR) );
							if( !IResult.second)
								(*IResult.first).second += 2*BR;
							branch_test_f += BR;
						}
						
					}
					
				}
				if (DM !=0)
					stable = false;
				// End of While loop
			}
			
			double btest = fabs(branch_test + branch_test_f-1.0);
			if ( btest > 1e-8 && !stable )
			{
				
				map< ZAI, double >::iterator DM_it = DaughtersMap.begin();
				if(branch_test+branch_test_f>0)
					for(  DM_it = DaughtersMap.begin();  DM_it != DaughtersMap.end(); DM_it++)
					{
						if ( (*DM_it).first != ZAI(-2,-2,-2) )
							(*DM_it).second *= 1./(branch_test+branch_test_f);
						else
							(*DM_it).second *= 1./(branch_test+branch_test_f);
					}
			}
			
			
			
			if (HalfLife < fShorstestHalflife)
			{
				fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ParentZAI, DaughtersMap ) );
			}
			else
			{
				ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >
								( ParentZAI, pair<double, map< ZAI, double > >
								  ( HalfLife, DaughtersMap) ) );
			}


		}
		
	} while (!infile.eof());
	
	{
		int i = 0;
		map<ZAI, pair<double, map< ZAI, double > > >::iterator it;
		for(it = ZAIDecay.begin() ; it != ZAIDecay.end(); it++)
		{
			findex.insert( pair<int, ZAI > ( i, (*it).first ) );
			findex_inver.insert( pair<ZAI, int > ( (*it).first , i ));
			i++;
		}
	}

	// Fill the Decay Part of the Bateman Matrix Always the same !
	bool FastDecayValidation  = false;
	while (!FastDecayValidation)
	{
		map<ZAI, map<ZAI, double> > FastDecayCopy = fFastDecay;
		
		map<ZAI, map<ZAI, double> >::iterator	FD_it;
		map<ZAI, map<ZAI, double> >::iterator	FD_it_Origin;
		
		
		for(FD_it = FastDecayCopy.begin(); FD_it != FastDecayCopy.end(); FD_it++)
		{
			
			FD_it_Origin = fFastDecay.find(FD_it->first);
			
			
			map<ZAI, double> BR = (*FD_it).second;
			map<ZAI, double>::iterator BR_it;
			
			for(BR_it = BR.begin(); BR_it != BR.end(); BR_it++)
			{
				
				map<ZAI, int>::iterator it_index = findex_inver.find( (*BR_it).first );
				
				if( it_index == findex_inver.end() )
				{
					map<ZAI, map<ZAI, double> >::iterator FD2_it = FastDecayCopy.find((*BR_it).first);
					if( FD2_it == fFastDecay.end() )
					{
						(*FD_it_Origin).second.erase((*BR_it).first);
						pair<map<ZAI, double>::iterator, bool> IResult;
						IResult = (*FD_it_Origin).second.insert( pair<ZAI, double> ( ZAI(-3,-3,-3), (*BR_it).second) );
						if( !IResult.second)
							(*IResult.first).second += (*BR_it).second ;
					}
					else
					{
						map<ZAI, double>::iterator BR2_it;
						(*FD_it_Origin).second.erase((*BR_it).first);
						
						for (BR2_it = (*FD2_it).second.begin(); BR2_it != (*FD2_it).second.end(); BR2_it++)
						{
							pair<map<ZAI, double>::iterator, bool> IResult;
							IResult = (*FD_it_Origin).second.insert( pair<ZAI, double> ( (*BR2_it).first, (*BR_it).second * (*BR2_it).second ) );
							if( !IResult.second)
								(*IResult.first).second += (*BR_it).second * (*BR2_it).second ;
						}
						
					}
					
				}
			}
			
		}
		
		FastDecayValidation =true;
		for(FD_it = fFastDecay.begin(); FD_it != fFastDecay.end(); FD_it++)
		{
			map<ZAI, double>::iterator BR_it;
			for (BR_it = (*FD_it).second.begin(); BR_it != (*FD_it).second.end(); BR_it++)
			{
				map<ZAI, int>::iterator Index_it = findex_inver.find( (*BR_it).first );
				map<ZAI, map<ZAI, double> >::iterator FD2_it = fFastDecay.find( (*BR_it).first );
				if(Index_it == findex_inver.end() && FD2_it == fFastDecay.end())
					FastDecayValidation = false;
			}
		}
	}


	fDecayMatrix.ResizeTo(findex.size(),findex.size());
	for(int i = 0; i < (int)findex.size(); i++)
		for(int j = 0; j < (int)findex.size(); j++)
			fDecayMatrix[i][j] = 0;




	{
		int i = 0;
		map<ZAI, pair<double, map< ZAI, double > > >::iterator it;
		for(it = ZAIDecay.begin() ; it != ZAIDecay.end(); it++)
		{
			map< ZAI, double >::iterator it2;
			map< ZAI, double > decaylist = (*it).second.second;
			for(it2 = decaylist.begin(); it2!= decaylist.end(); it2++)
			{
				
				map<ZAI, int >::iterator it3 = findex_inver.find( (*it2).first );
				if( it3 != findex_inver.end() )
				{
					fDecayMatrix[(*it3).second][i] += log(2.)/(*it).second.first * (*it2).second;
				}
				else
				{
					map<ZAI, map<ZAI, double> >::iterator it4 = fFastDecay.find( (*it2).first );
					
					if( it4 == fFastDecay.end() )
					{


						fDecayMatrix[0][i] += log(2.)/(*it).second.first * (*it2).second;
					}
					else
					{
						map< ZAI, double >::iterator it5;
						map< ZAI, double > decaylist2 = (*it4).second;
						for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
						{
							it3 = findex_inver.find( (*it5).first );
							if( it3 == findex_inver.end() )
								fDecayMatrix[0][i] += log(2.)/(*it).second.first * (*it2).second * (*it5).second;
							
							else
							{
								fDecayMatrix[(*it3).second][i] += log(2.)/(*it).second.first * (*it2).second * (*it5).second;

							}
						}
					}
					
				}
			}
			fDecayMatrix[i][i] += -log(2.)/(*it).second.first;

			i++;
			
			
		}
	}
	//exit(1);
	
}


template<>
void DataBank<IsotopicVector>::SetTheMatrixToZero()
{
	//ResetTheMatrix();
	
	fNVar = findex.size();
	fTheMatrix = new double*[fNVar];
	
#pragma omp parallel for
	for(int i= 0; i < fNVar; i++)
		fTheMatrix[i] = new double[fNVar];
	
	for(int i = 0; i < fNVar; i++)
	{
		for(int k = 0; k < fNVar; k++)
		{
			fTheMatrix[i][k]=0.0;
		}
	}

}

template<>
void DataBank<IsotopicVector>::ResetTheMatrix()
{
	
	if(fTheMatrix)
	{
		for(int i= 0; i<fNVar; i++)
			delete [] fTheMatrix[i];
		delete [] fTheMatrix;
	}
	fTheMatrix = 0;
}


template<>
void DataBank<IsotopicVector>::SetTheNucleiVectorToZero()
{
	ResetTheNucleiVector();
	fTheNucleiVector = new double[fNVar];

#pragma omp parallel for
	for(int i = 0; i < fNVar; i++)
			fTheNucleiVector[i]=0.0;
	
}

template<>
void DataBank<IsotopicVector>::ResetTheNucleiVector()
{
	if(fTheNucleiVector)
		delete [] fTheNucleiVector;
	fTheNucleiVector = 0;
}


template<>
void DataBank<IsotopicVector>::BuildEqns(double t, double *N, double *dNdt)
{
	double sum=0;
	// pragma omp parallel for reduction(+:sum)
	for(int i = 0; i < fNVar; i++)
	{
		sum=0;
		for(int k = 0; k < fNVar; k++)
		{
			sum += fTheMatrix[i][k]*N[k];
		}
		dNdt[i] = sum;
	}
}

template<>
void DataBank<IsotopicVector>::SetTheMatrix(TMatrixT<double> BatemanMatrix)
{
	for (int k = 0; k < (int)fNVar; k++)
		for (int l = 0; l < (int)findex_inver.size(); l++)
			fTheMatrix[l][k] = BatemanMatrix[l][k];
}

template<>
TMatrixT<double> DataBank<IsotopicVector>::GetTheMatrix()
{
	TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
	for (int k = 0; k < (int)fNVar; k++)
		for (int l = 0; l < (int)findex_inver.size(); l++)
			BatemanMatrix[l][k] = fTheMatrix[l][k];

	return BatemanMatrix;
}


template<>
void DataBank<IsotopicVector>::SetTheNucleiVector(TMatrixT<double> NEvolutionMatrix)
{
	for (int k = 0; k < (int)fNVar; k++)
		fTheNucleiVector[k] = NEvolutionMatrix[k][0];
}

template<>
TMatrixT<double> DataBank<IsotopicVector>::GetTheNucleiVector()
{
	TMatrixT<double> NEvolutionMatrix = TMatrixT<double>(findex.size(),1);
	for (int k = 0; k < (int)fNVar; k++)
		NEvolutionMatrix[k][0] = fTheNucleiVector[k];

	return NEvolutionMatrix;
}


//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________

/*
 
 
 OLD READ
 
 
 */

//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________

































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































































template<>
void DataBank<IsotopicVector>::OldBuildDecayMatrix()
{
	// List of Decay Time and Properties
	map<ZAI, pair<double, map< ZAI, double > > > ZAIDecay;
	
	{	// TMP
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(-3,-3,-3), pair<double, map< ZAI, double > > ( 1e28 ,toAdd )) ) ;
	}
	{	// PF
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-2,-2,-2), 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(-2,-2,-2), pair<double, map< ZAI, double > > ( 1e28 ,toAdd )) ) ;
	}
	{	// 232Th
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(90,232,0), pair<double, map< ZAI, double > > ( 4.41806400000000000e+17 , toAdd ) ) );
	}
	{	// 233U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,233,0), pair<double, map< ZAI, double > > ( 5.02396992000000000e+12, toAdd) ) );
	}
	{	// 234U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,234,0), pair<double, map< ZAI, double > > ( 7.74739080000000000e+12, toAdd) ) );
	}
	{	// 235U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,235,0), pair<double, map< ZAI, double > > ( 2.22165504000000000e+16, toAdd) ) );
	}
	{	// 236U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(90,232,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,236,0), pair<double, map< ZAI, double > > ( 7.39078992000000000e+14, toAdd) ) );
	}
	{	// 238U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,238,0), pair<double, map< ZAI, double > > ( 1.40999356800000000e+17, toAdd) ) );
	}
	{	// 237Np
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(91,233,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(93,237,0), pair<double, map< ZAI, double > > ( 6.76594944000000000e+13, toAdd) ) );
	}
	{	// 238Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,234,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,238,0), pair<double, map< ZAI, double > > ( 2.76760152000000000e+09, toAdd) ) );
	}
	{	// 239Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,235,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,239,0), pair<double, map< ZAI, double > > ( 7.60853736000000000e+11, toAdd) ) );
	}
	{	// 240Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,236,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,240,0), pair<double, map< ZAI, double > > ( 2.07049413600000000e+11, toAdd) ) );
	}
	{	// 241Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(95,241,0) , (1-2.5e-5)) );
		toAdd.insert(pair<ZAI, double> ( ZAI(92,237,0) , 2.5e-5) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,241,0), pair<double, map< ZAI, double > > ( 4.52062620000000000e+08, toAdd) ) );
	}
	{	// 242Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,238,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,242,0), pair<double, map< ZAI, double > > ( 1.18341000000000000e+13, toAdd) ) );
	}
	{	// 241Am
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(93,237,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(95,241,0), pair<double, map< ZAI, double > > ( 1.36518177600000000e+10, toAdd) ) );
	}
	{	// 242Am*
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(93,238,0) , 0.00459) );
		toAdd.insert(pair<ZAI, double> ( ZAI(95,242,0) , 0.99541) );
		
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(95,242,1), pair<double, map< ZAI, double > > ( 4.44962160000000000e+09, toAdd) ) );
	}
	{	// 243Am
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,239,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(95,243,0), pair<double, map< ZAI, double > > ( 2.32579512000000000e+11, toAdd) ) );
	}
	{	// 242Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,238,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,242,0), pair<double, map< ZAI, double > > ( 1.40659200000000000e+07 , toAdd) ) );
	}
	{	// 243Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,239,0) , 0.9971) );
		toAdd.insert(pair<ZAI, double> ( ZAI(95,243,0) , 0.0029) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,243,0), pair<double, map< ZAI, double > > ( 9.18326160000000000e+08, toAdd) ) );
	}
	{	// 244Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,240,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,244,0), pair<double, map< ZAI, double > > ( 5.71192560000000000e+08, toAdd) ) );
	}
	{	// 245Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,241,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,245,0), pair<double, map< ZAI, double > > ( 2.65809664800000000e+11, toAdd) ) );
	}
	{	// 246Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,242,0) , 0.9997) );
		toAdd.insert(pair<ZAI, double> ( ZAI(-2,-2,-2) , 0.0003) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,246,0), pair<double, map< ZAI, double > > ( 1.48510065600000000e+11, toAdd) ) );
	}
	{	// 247Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,243,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,247,0), pair<double, map< ZAI, double > > ( 4.92298560000000000e+14, toAdd) ) );
	}
	{	// 248Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 0.9161) );
		toAdd.insert(pair<ZAI, double> ( ZAI(-2,-2,-2) , 0.0839) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,248,0), pair<double, map< ZAI, double > > ( 1.09820448000000000e+13, toAdd) ) );
	}
	
	{	// 231Th
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(90,231,0), toAdd ) );
	}
	{	// 233Th
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,233,0) , 1) );
		
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(90,233,0), toAdd ) );
	}
	{	// 233Pa
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,233,0) , 1) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(91,233,0), toAdd ) );
	}
	{	// 237U
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(93,237,0) , 1) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(92,237,0), toAdd ) );
	}
	{	// 239U
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,239,0) , 1) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(92,239,0), toAdd ) );
	}
	{	// 238Np
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,238,0) , 1) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(93,238,0), toAdd ) );
	}
	{	// 239Np
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,239,0) , 1) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(93,239,0), toAdd ) );
	}
	{	// 240Np
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,240,0) , 1) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(93,240,0), toAdd ) );
	}
	{	// 241Np
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,241,0) , 1) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(93,241,0), toAdd ) );
	}
	{	// 237Pu
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(93,237,0) , 1) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(94,237,0), toAdd ) );
	}
	{	// 243Pu
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(95,243,0) , 1) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(94,243,0), toAdd ) );
	}
	{	// 240Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,240,0) , 1) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,240,0), toAdd ) );
	}
	{	// 242Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(96,242,0) , 0.827) );
		toAdd.insert(pair<ZAI, double> ( ZAI(94,242,0) , 0.173) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,242,0), toAdd ) );
	}
	{	// 244Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(96,244,0) , 1) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,244,0), toAdd ) );
	}
	{	// 245Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(96,245,0) , 1) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,245,0), toAdd ) );
	}
	{	// 249Cm
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(96,249,0), toAdd ) );
	}
	
	
	// Build Index avd reverse index, correspondance ZAI-nt (and int-ZAI)
	{
		int i = 0;
		map<ZAI, pair<double, map< ZAI, double > > >::iterator it;
		for(it = ZAIDecay.begin() ; it != ZAIDecay.end(); it++)
		{
			findex.insert( pair<int, ZAI > ( i, (*it).first ) );
			findex_inver.insert( pair<ZAI, int > ( (*it).first , i ));
			i++;
		}
	}
	
	fDecayMatrix.ResizeTo(findex.size(),findex.size());
	for(int i = 0; i < (int)findex.size(); i++)
		for(int j = 0; j < (int)findex.size(); j++)
			fDecayMatrix[i][j] = 0;
	
	
	// Fill the Decay Part of the Bateman Matrix Always the same !
	{
		int i = 0;
		map<ZAI, pair<double, map< ZAI, double > > >::iterator it;
		for(it = ZAIDecay.begin() ; it != ZAIDecay.end(); it++)
		{
			map< ZAI, double >::iterator it2;
			map< ZAI, double > decaylist = (*it).second.second;
			for(it2 = decaylist.begin(); it2!= decaylist.end(); it2++)
			{
				
				map<ZAI, int >::iterator it3 = findex_inver.find( (*it2).first );
				if( it3 != findex_inver.end() )
					fDecayMatrix[(*it3).second][i] = log(2.)/(*it).second.first * (*it2).second;
				else
				{
					map<ZAI, map<ZAI, double> >::iterator it4 = fFastDecay.find( (*it2).first );
					
					if( it4 == fFastDecay.end() )
					{
						cout << "Problem in FastDecay for nuclei " << (*it2).first.Z() << " " << (*it2).first.A() << " " << (*it2).first.I() << endl;
						exit(1);
					}
					
					map< ZAI, double >::iterator it5;
					map< ZAI, double > decaylist2 = (*it4).second;
					for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
					{
						it3 = findex_inver.find( (*it5).first );
						if( it3 == findex_inver.end() )
						{
							cout << "Problem in FastDecay for nuclei " << (*it2).first.Z() << " " << (*it2).first.A() << " " << (*it2).first.I() << endl;
							exit(1);
						}
						fDecayMatrix[(*it3).second][i] = log(2.)/(*it).second.first * (*it2).second * (*it5).second;
					}
					
				}
			}
			fDecayMatrix[i][i] += -log(2.)/(*it).second.first;
			i++;
			
			
		}
	}
	
	
	
	
	
}



template<>
EvolutionData DataBank<IsotopicVector>::OldGenerateEvolutionData(IsotopicVector isotopicvector, double cycletime, double Power)
{
	
	string ReactorType;
	map<ZAI, pair<double, map< ZAI, double > > > ZAIDecay;
	
	{	// TMP
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(-3,-3,-3), pair<double, map< ZAI, double > > ( 1e28 ,toAdd )) ) ;
	}
	{	// PF
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-2,-2,-2), 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(-2,-2,-2), pair<double, map< ZAI, double > > ( 1e28 ,toAdd )) ) ;
	}
	{	// 232Th
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(90,232,0), pair<double, map< ZAI, double > > ( 2.37944304000000000e+18 , toAdd ) ) );
	}
	{	// 233U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,233,0), pair<double, map< ZAI, double > > ( 5.02396992000000000e+12, toAdd) ) );
	}
	{	// 234U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,234,0), pair<double, map< ZAI, double > > ( 7.74739080000000000e+12, toAdd) ) );
	}
	{	// 235U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,235,0), pair<double, map< ZAI, double > > ( 2.22165504000000000e+16, toAdd) ) );
	}
	{	// 236U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(90,232,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,236,0), pair<double, map< ZAI, double > > ( 7.39078992000000000e+14, toAdd) ) );
	}
	{	// 238U
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(92,238,0), pair<double, map< ZAI, double > > ( 1.40999356800000000e+17, toAdd) ) );
	}
	{	// 237Np
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(91,233,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(93,237,0), pair<double, map< ZAI, double > > ( 6.76594944000000000e+13, toAdd) ) );
	}
	{	// 238Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,234,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,238,0), pair<double, map< ZAI, double > > ( 2.76760152000000000e+09, toAdd) ) );
	}
	{	// 239Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,235,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,239,0), pair<double, map< ZAI, double > > ( 7.60853736000000000e+11, toAdd) ) );
	}
	{	// 240Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,236,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,240,0), pair<double, map< ZAI, double > > ( 2.07049413600000000e+11, toAdd) ) );
	}
	{	// 241Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(95,241,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,241,0), pair<double, map< ZAI, double > > ( 4.52062620000000000e+08, toAdd) ) );
	}
	{	// 242Pu
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,238,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(94,242,0), pair<double, map< ZAI, double > > ( 1.18341000000000000e+13, toAdd) ) );
	}
	{	// 241Am
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(93,237,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(95,241,0), pair<double, map< ZAI, double > > ( 1.36518177600000000e+10, toAdd) ) );
	}
	{	// 242Am*
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(93,238,0) , 0.00459) );
		toAdd.insert(pair<ZAI, double> ( ZAI(95,242,0) , 0.99541) );
		
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(95,242,1), pair<double, map< ZAI, double > > ( 4.44962160000000000e+09, toAdd) ) );
	}
	{	// 243Am
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,239,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(95,243,0), pair<double, map< ZAI, double > > ( 2.32579512000000000e+11, toAdd) ) );
	}
	{	// 242Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,238,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,242,0), pair<double, map< ZAI, double > > ( 1.40659200000000000e+07 , toAdd) ) );
	}
	{	// 243Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,239,0) , 0.9971) );
		toAdd.insert(pair<ZAI, double> ( ZAI(95,243,0) , 0.0029) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,243,0), pair<double, map< ZAI, double > > ( 9.18326160000000000e+08, toAdd) ) );
	}
	{	// 244Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,240,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,244,0), pair<double, map< ZAI, double > > ( 5.71192560000000000e+08, toAdd) ) );
	}
	{	// 245Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,241,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,245,0), pair<double, map< ZAI, double > > ( 2.65809664800000000e+11, toAdd) ) );
	}
	{	// 246Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,242,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,246,0), pair<double, map< ZAI, double > > ( 1.48510065600000000e+11, toAdd) ) );
	}
	{	// 247Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,243,0) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,247,0), pair<double, map< ZAI, double > > ( 4.92298560000000000e+14, toAdd) ) );
	}
	{	// 248Cm
		map< ZAI, double > toAdd;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >( ZAI(96,248,0), pair<double, map< ZAI, double > > ( 1.09820448000000000e+13, toAdd) ) );
	}
	
	map<ZAI, map<ZAI, double> > FastDecay;
	{	// 231Th
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(90,231,0), toAdd ) );
	}
	{	// 233Th
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,233,0) , 1) );
		
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(90,233,0), toAdd ) );
	}
	{	// 233Pa
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(92,233,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(91,233,0), toAdd ) );
	}
	{	// 237U
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(93,237,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(92,237,0), toAdd ) );
	}
	{	// 239U
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,239,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(92,239,0), toAdd ) );
	}
	{	// 238Np
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,238,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(93,238,0), toAdd ) );
	}
	{	// 239Np
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,239,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(93,239,0), toAdd ) );
	}
	{	// 240Np
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,240,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(93,240,0), toAdd ) );
	}
	{	// 241Np
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,241,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(93,241,0), toAdd ) );
	}
	{	// 237Pu
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(93,237,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(94,237,0), toAdd ) );
	}
	{	// 243Pu
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(95,243,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(94,243,0), toAdd ) );
	}
	{	// 240Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(94,240,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,240,0), toAdd ) );
	}
	{	// 242Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(96,242,0) , 0.827) );
		toAdd.insert(pair<ZAI, double> ( ZAI(94,242,0) , 0.173) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,242,0), toAdd ) );
	}
	{	// 244Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(96,244,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,244,0), toAdd ) );
	}
	{	// 245Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(96,245,0) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,245,0), toAdd ) );
	}
	{	// 249Cm
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(-3,-3,-3) , 1) );
		FastDecay.insert( pair< ZAI, map<ZAI, double> > ( ZAI(96,249,0), toAdd ) );
	}
	
	
	map<ZAI, map<ZAI, double> > Capture;
	{	// 241Am
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(95,242,0) , 0.086) );
		toAdd.insert(pair<ZAI, double> ( ZAI(95,242,1) , 0.914) );
		Capture.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,241,0), toAdd ) );
	}
	{	// 242Am*
		map<ZAI, double> toAdd ;
		toAdd.insert(pair<ZAI, double> ( ZAI(95,243,0) , 1) );
		Capture.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,242,1), toAdd ) );
	}
	map<ZAI, int> index_inver;
	map<int, ZAI> index;
	{
		int i = 0;
		map<ZAI, pair<double, map< ZAI, double > > >::iterator it;
		for(it = ZAIDecay.begin() ; it != ZAIDecay.end(); it++)
		{
			index.insert( pair<int, ZAI > ( i, (*it).first ) );
			index_inver.insert( pair<ZAI, int > ( (*it).first , i ));
			i++;
		}
	}
	
	TMatrixT<double> DecayMatrix = TMatrixT<double>(index.size(),index.size());
	for(int i = 0; i < (int)index.size(); i++)
		for(int j = 0; j < (int)index.size(); j++)
			DecayMatrix[i][j] = 0;
	
	
	// Fill the Decay Part of the Bateman Matrix
	{
		int i = 0;
		map<ZAI, pair<double, map< ZAI, double > > >::iterator it;
		for(it = ZAIDecay.begin() ; it != ZAIDecay.end(); it++)
		{
			map< ZAI, double >::iterator it2;
			map< ZAI, double > decaylist = (*it).second.second;
			for(it2 = decaylist.begin(); it2!= decaylist.end(); it2++)
			{
				
				map<ZAI, int >::iterator it3 = index_inver.find( (*it2).first );
				if( it3 != index_inver.end() )
					DecayMatrix[(*it3).second][i] = log(2.)/(*it).second.first * (*it2).second;
				else
				{
					map<ZAI, map<ZAI, double> >::iterator it4 = FastDecay.find( (*it2).first );
					
					if( it4 == FastDecay.end() )
					{
						cout << "Problem in FastDecay for nuclei " << (*it2).first.Z() << " " << (*it2).first.A() << " " << (*it2).first.I() << endl;
						exit(1);
					}
					
					map< ZAI, double >::iterator it5;
					map< ZAI, double > decaylist2 = (*it4).second;
					for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
					{
						it3 = index_inver.find( (*it5).first );
						if( it3 == index_inver.end() )
						{
							cout << "Problem in FastDecay for nuclei " << (*it2).first.Z() << " " << (*it2).first.A() << " " << (*it2).first.I() << endl;
							exit(1);
						}
						DecayMatrix[(*it3).second][i] = log(2.)/(*it).second.first * (*it2).second * (*it5).second;
					}
					
				}
			}
			DecayMatrix[i][i] += -log(2.)/(*it).second.first;
			i++;
			
			
		}
	}
	
	
	
	
	//-------------------------//
	//--- Perform Evolution ---//
	//-------------------------//
	double timevector[17];
	timevector[0] = 0.;
	vector< TMatrixT<double> > NMatrix ;//  TMatrixT<double>(decayindex.size(),1))
	{	// Filling the t=0 State;
		map<ZAI, double > isotopicquantity = isotopicvector.GetIsotopicQuantity();
		TMatrixT<double>  N_0Matrix =  TMatrixT<double>( index.size(),1) ;
		
		map<ZAI, double >::iterator it ;
		for(int i = 0; i < (int)index.size(); i++)
			N_0Matrix[i] = 0;
		
		for(it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		{
			
			map<ZAI, int >::iterator it2;
			
			if( (*it).first.Z() < 90 )
				it2 = index_inver.find( ZAI(-2,-2,-2) );
			else it2 = index_inver.find( (*it).first );
			
			if(it2 == index_inver.end() )				//If not in index should be TMP, can't be fast decay for new Fuel !!!
				it2 = index_inver.find( ZAI(-3,-3,-3) );
			
			N_0Matrix[ (*it2).second ][0] = (*it).second ;
			
			
		}
		NMatrix.push_back(N_0Matrix);
	}
	
	
	TMatrixT<double> SigmaPhi = TMatrixT<double>(index.size()*3+1,16);
	
	EvolutionData EvolutionDataStep = GetClosest(isotopicvector.GetActinidesComposition(), 0.);	//GetCLosest at the begining of evolution
	
	
	ReactorType = EvolutionDataStep.GetReactorType();
	
	for(int i = 0; i < 16; i++)
	{
		
		
		double TStep = cycletime/16*i;
		
		TMatrixT<double> BatemanMatrix = TMatrixT<double>(index.size(),index.size());
		BatemanMatrix = DecayMatrix ;
		
		
		IsotopicVector IVStep;
		for(int k = 0; k < (int)index.size(); k++)
			IVStep += index.find(k)->second * NMatrix.back()[k][0];
		
		if(i != 0);		//GetCLosest at the each of evolution step (begining already done...)
		EvolutionDataStep = GetClosest(IVStep, TStep);
		
		double NormFactor = 1;
		{
			IsotopicVector WantedHMIV = isotopicvector.GetSpeciesComposition(90)
			+ isotopicvector.GetSpeciesComposition(92)
			+ isotopicvector.GetSpeciesComposition(93)
			+ isotopicvector.GetSpeciesComposition(94)
			+ isotopicvector.GetSpeciesComposition(95)
			+ isotopicvector.GetSpeciesComposition(96);
			
			IsotopicVector DBHMIV = EvolutionDataStep.GetIsotopicVectorAt(0).GetSpeciesComposition(90)
			+ EvolutionDataStep.GetIsotopicVectorAt(0).GetSpeciesComposition(92)
			+ EvolutionDataStep.GetIsotopicVectorAt(0).GetSpeciesComposition(93)
			+ EvolutionDataStep.GetIsotopicVectorAt(0).GetSpeciesComposition(94)
			+ EvolutionDataStep.GetIsotopicVectorAt(0).GetSpeciesComposition(95)
			+ EvolutionDataStep.GetIsotopicVectorAt(0).GetSpeciesComposition(96);
			
			NormFactor = Norme(WantedHMIV)/ Norme(DBHMIV);
		}
		
		double Flux = EvolutionDataStep.GetFlux()->Eval(TStep)*Power/(EvolutionDataStep.GetPower()*NormFactor);
		SigmaPhi[index.size()*3][i] = Flux;
		
		
		map<ZAI ,TGraph* >::iterator it;
		// ----------------  A(n,.) X+Y
		
		map<ZAI ,TGraph* > FissionXS = EvolutionDataStep.GetFissionXS();
		
		for(it = FissionXS.begin() ; it != FissionXS.end(); it++)
		{
			
			if( index_inver.find( (*it).first ) != index_inver.end() )
			{
				double y;
				y = (*it).second->Eval(TStep);
				
				BatemanMatrix[ index_inver.find( (*it).first )->second ][index_inver.find( (*it).first )->second] += -y* 1e-24 *Flux;
				BatemanMatrix[1][ index_inver.find( (*it).first )->second] += 2*y* 1e-24 *Flux;
				
				SigmaPhi[index_inver.find( (*it).first )->second][i] = y ;
			}
			
		}
		
		// ----------------  A(n,.)A+1
		map<ZAI ,TGraph* > CaptureXS = EvolutionDataStep.GetCaptureXS();
		for(it = CaptureXS.begin(); it != CaptureXS.end(); it++)
		{
			if( index_inver.find( (*it).first ) != index_inver.end() )
			{
				double y;
				y = (*it).second->Eval(TStep);
				
				BatemanMatrix[index_inver.find( (*it).first )->second][ index_inver.find( (*it).first )->second ] += -y* 1e-24 *Flux;
				SigmaPhi[index_inver.find( (*it).first )->second + index.size() ][i] = y ;
				
				map<ZAI, map<ZAI, double> >::iterator it3 = Capture.find( (*it).first );
				
				if( it3 == Capture.end() )
				{
					map<ZAI, int >::iterator it6 = index_inver.find( ZAI( (*it).first.Z(), (*it).first.A()+1, (*it).first.I()) );
					
					if( it6 != index_inver.end() )
					{
						BatemanMatrix[(*it6).second][index_inver.find( (*it).first )->second] += y* 1e-24 *Flux ;
					}
					else
					{
						map<ZAI, map<ZAI, double> >::iterator it4 = FastDecay.find(  ZAI( (*it).first.Z(), (*it).first.A()+1, (*it).first.I()) );
						
						if( it4 == FastDecay.end() )
						{
							cout << "Problem in FastDecay for nuclei " << (*it).first.Z() << " " << (*it).first.A()+1 << " " << (*it).first.I() << endl;
							exit(1);
						}
						
						map< ZAI, double >::iterator it5;
						map< ZAI, double > decaylist2 = (*it4).second;
						for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
						{
							it6 = index_inver.find( (*it5).first );
							if( it6 == index_inver.end() )
							{
								cout << "Problem in FastDecay for nuclei " << (*it).first.Z() << " " << (*it).first.A() << " " << (*it).first.I() << endl;
								exit(1);
							}
							BatemanMatrix[(*it6).second][index_inver.find( (*it).first )->second] += y* 1e-24 *Flux * (*it5).second;
						}
					}
				}
				else
				{
					map<ZAI, double>::iterator it4;
					map<ZAI, double> CaptureList = (*it3).second;
					for(it4 = CaptureList.begin(); it4 != CaptureList.end() ; it4++)
					{
						
						map<ZAI, int >::iterator it6 = index_inver.find( (*it4).first );
						if( it6 != index_inver.end() )
							BatemanMatrix[(*it6).second][index_inver.find( (*it).first )->second] += y* 1e-24 *Flux * (*it4).second ;
						else
						{
							map<ZAI, map<ZAI, double> >::iterator it7 = FastDecay.find( (*it4).first );
							
							if( it7 == FastDecay.end() )
							{
								cout << "Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
								exit(1);
							}
							
							map< ZAI, double >::iterator it5;
							map< ZAI, double > decaylist2 = (*it7).second;
							for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
							{
								
								it6 = index_inver.find( (*it5).first );
								if( it6 == index_inver.end() )
								{
									cout << "Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
									exit(1);
								}
								BatemanMatrix[(*it6).second][index_inver.find( (*it).first )->second] += y * 1e-24 * Flux * (*it5).second * (*it4).second;
							}
						}
						
					}
				}
				
				
			}
		}
		
		// ----------------  A(n,2n)A-1
		map<ZAI ,TGraph* > n2nXS = EvolutionDataStep.Getn2nXS();
		for(it = n2nXS.begin() ; it != n2nXS.end(); it++)
		{
			if( index_inver.find( (*it).first ) != index_inver.end() )
			{
				double y;
				y = (*it).second->Eval(TStep);
				BatemanMatrix[ index_inver.find( (*it).first )->second ][index_inver.find( (*it).first )->second] += -y* 1e-24 *Flux;
				SigmaPhi[index_inver.find( (*it).first )->second + index.size() + index.size()][i] = y ;
				
				
				map<ZAI, int>::iterator it3 = index_inver.find( ZAI( (*it).first.Z(), (*it).first.A()-1, 0) );
				
				if( it3 != index_inver.end() )
					BatemanMatrix[(*it3).second][index_inver.find( (*it).first )->second] += y* 1e-24 *Flux;
				else
				{
					
					map<ZAI, map<ZAI, double> >::iterator it4 = FastDecay.find( ZAI( (*it).first.Z(), (*it).first.A()-1, 0) );
					
					if( it4 == FastDecay.end() )
					{
						it3 = index_inver.find( ZAI( -3, -3, -3 ) );
						BatemanMatrix[(*it3).second][index_inver.find( (*it).first )->second] += y* 1e-24 *Flux;
					}
					else
					{
						map< ZAI, double >::iterator it5;
						map< ZAI, double > decaylist2 = (*it4).second;
						for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
						{
							
							it3 = index_inver.find( (*it5).first );
							if( it3 == index_inver.end() )
							{
								cout << "Problem in FastDecay for nuclei " << (*it4).first.Z() << " " << (*it4).first.A() << " " << (*it4).first.I() << endl;
								exit(1);
							}
							BatemanMatrix[(*it3).second][index_inver.find( (*it).first )->second] += y* 1e-24 *Flux * (*it5).second ;
						}
					}
				}
			}
		}
		
		// ----------------   Evolution
		TMatrixT<double> NEvolutionMatrix = TMatrixT<double>(index.size(),1);
		
		double TStepMax = cycletime/16.;
		timevector[i+1] = timevector[i] + TStepMax;
		
		BatemanMatrix *= TStepMax;
		TMatrixT<double> IdMatrix = TMatrixT<double>(index.size(),index.size());
		for(int j = 0; j < (int)index.size(); j++)
			for(int k = 0; k < (int)index.size(); k++)
			{
				if(k == j)	IdMatrix[j][k] = 1;
				else 		IdMatrix[j][k] = 0;
			}
		
		
		TMatrixT<double> BatemanMatrixDL = TMatrixT<double>(index.size(),index.size());   // Order 0 Term from the DL : Id
		TMatrixT<double> BatemanMatrixDLTermN = TMatrixT<double>(index.size(),index.size());  // Addind it;
		
		{
			BatemanMatrixDLTermN = IdMatrix;
			BatemanMatrixDL = BatemanMatrixDLTermN;
			
			
			int j = 1;
			double NormN = 0;
			
			do
			{
				
				NormN = 0;
				TMatrixT<double> BatemanMatrixDLTermtmp = TMatrixT<double>(index.size(),index.size());  // Adding it;
				BatemanMatrixDLTermtmp = BatemanMatrixDLTermN;
				BatemanMatrixDLTermN.Mult(BatemanMatrixDLTermtmp, BatemanMatrix );
				
				BatemanMatrixDLTermN *= 1./j;
				BatemanMatrixDL += BatemanMatrixDLTermN;
				
				NormN = 0;
				for(int m = 0; m < (int)index.size(); m++)
					for(int n = 0; n < (int)index.size(); n++)
						NormN += BatemanMatrixDLTermN[m][n]*BatemanMatrixDLTermN[m][n];
				j++;
			} while ( NormN != 0 );
		}
		
		NEvolutionMatrix = BatemanMatrixDL * NMatrix.back() ;
		NMatrix.push_back(NEvolutionMatrix);
	}
	
	
	EvolutionData GeneratedDB = EvolutionData(GetLog());
	double Flux[16];
	for(int j = 0; j < 16; j++)
		Flux[j] = SigmaPhi[index.size()*3][j];
	GeneratedDB.SetFlux( new TGraph(16, timevector, Flux)  );
	
	for(int i = 0; i < (int)index.size(); i++)
	{
		double ZAIQuantity[NMatrix.size()];
		double FissionXS[16];
		double CaptureXS[16];
		double n2nXS[16];
		for(int j = 0; j < (int)NMatrix.size(); j++)
			ZAIQuantity[j] = (NMatrix[j])[i][0];
		
		for(int j = 0; j < 16; j++)
		{
			FissionXS[j]	= SigmaPhi[i][j];
			CaptureXS[j]	= SigmaPhi[i + index.size()][j];
			n2nXS[j]	= SigmaPhi[i + index.size() + index.size()][j];
		}
		
		GeneratedDB.NucleiInsert(pair<ZAI, TGraph*> (index.find(i)->second, new TGraph(NMatrix.size(), timevector, ZAIQuantity) ) );
		GeneratedDB.FissionXSInsert(pair<ZAI, TGraph*> (index.find(i)->second, new TGraph(16, timevector, FissionXS) ) );
		GeneratedDB.CaptureXSInsert(pair<ZAI, TGraph*> (index.find(i)->second, new TGraph(16, timevector, CaptureXS) ) );
		GeneratedDB.n2nXSInsert(pair<ZAI, TGraph*> (index.find(i)->second, new TGraph(16, timevector, n2nXS) ) );
	}
	
	GeneratedDB.SetPower(Power );
	GeneratedDB.SetFuelType(fFuelType );
	GeneratedDB.SetReactorType(ReactorType );
	GeneratedDB.SetCycleTime(cycletime);
	
	
	return GeneratedDB;
	
}
