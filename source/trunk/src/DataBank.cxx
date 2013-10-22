#include "DataBank.hxx"

#include "IsotopicVector.hxx"
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






	//________________________________________________________________________
	//________________________________________________________________________
	//________________________________________________________________________
	//________________________________________________________________________
	//________________________________________________________________________
template<>
void DataBank<IsotopicVector>::ReadDataBase();
template<>
TMatrixT<double> DataBank<IsotopicVector>::GetFissionXsMatrix(EvolutionData EvolutionDataStep,double TStep);
template<>
TMatrixT<double> DataBank<IsotopicVector>::GetCaptureXsMatrix(EvolutionData EvolutionDataStep,double TStep);
template<>
TMatrixT<double> DataBank<IsotopicVector>::Getn2nXsMatrix(EvolutionData EvolutionDataStep,double TStep);
template<>
TMatrixT<double> DataBank<IsotopicVector>::ExtractXS(EvolutionData EvolutionDataStep,double TStep);
template<>
void DataBank<IsotopicVector>::BuildDecayMatrix();
template<>
EvolutionData DataBank<IsotopicVector>::OldGenerateEvolutionData(IsotopicVector isotopicvector, double cycletime, double Power);

	//________________________________________________________________________


template<>
DataBank<IsotopicVector>::DataBank(LogFile* Log, string DB_index_file, bool setlog, bool olfreadmethod)
{
	
	SetLog(Log);
	IsLog(setlog);

	fDataBaseIndex = DB_index_file;
	fUpdateReferenceDBatEachStep = false;
	fOldReadMethod = olfreadmethod;
	fUseOldGeneration = false;
	fDistanceType = 0;
	BuildDecayMatrix();
	ReadDataBase();
	
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
				   / Norme( evolutiondb.begin()->second.GetIsotopicVectorAt(t).GetActinidesComposition() )
				   * Norme(isotopicvector.GetActinidesComposition()),
				   fDistanceType, fDistanceParameter);
	
	EvolutionData CloseEvolData = evolutiondb.begin()->second ;
	
	map<IsotopicVector, EvolutionData >::iterator it;
	for( it = evolutiondb.begin(); it != evolutiondb.end(); it++ )
	{
		pair<map<double, EvolutionData>::iterator, bool> IResult;
		double D = Distance(isotopicvector.GetActinidesComposition(),
				    (*it).second.GetIsotopicVectorAt(t).GetActinidesComposition()
				    / Norme( (*it).second.GetIsotopicVectorAt(t).GetActinidesComposition() )
				    * Norme(isotopicvector.GetActinidesComposition()),
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
	
	
	/*{
		map<IsotopicVector, EvolutionData>::iterator it;
		for( it = fDataBankCalculated.begin(); it != fDataBankCalculated.end(); it++)
		{
			if (Distance(isotopicvector, (*it).first) == 0
			    && cycletime == (*it).second.GetCycleTime()
			    && Power == (*it).second.GetPower() ) {
				return (*it).second;
			}
		}
	}*/
	
	if(fUseOldGeneration)
	{
		EvolutionData GeneratedDB = OldGenerateEvolutionData( isotopicvector,  cycletime, Power);
		fDataBankCalculated.insert( pair< IsotopicVector, EvolutionData > ( GeneratedDB.GetIsotopicVectorAt(0.), GeneratedDB) );
		return GeneratedDB;
	
	}
	string ReactorType;
	
		//-------------------------//
		//--- Perform Evolution ---//
		//-------------------------//
	double timevector[17];
	timevector[0] = 0.;
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
	
	
	TMatrixT<double> SigmaPhi = TMatrixT<double>(findex.size()*3+1,16); // Store the XS and the flux trought the evolution calculation.
	
	EvolutionData EvolutionDataStep = GetClosest(isotopicvector.GetActinidesComposition(), 0.);	//GetCLosest at the begining of evolution
	
	ReactorType = EvolutionDataStep.GetReactorType();
	TMatrixT<double> FissionEnergy = TMatrixT<double>(findex.size(),1);
	{
		map< ZAI, int >::iterator it;
		for(it = findex_inver.begin(); it != findex_inver.end(); it++)
		{
			map< ZAI, double >::iterator it2 = fFissionEnergy.find(it->first);
			if(it2 == fFissionEnergy.end())
				FissionEnergy[it->second][0] = 200.e6;
			else
				FissionEnergy[it->second][0] = it2->second;
		}
	}
	
	vector< TMatrixT<double> > FissionXSMatrix; //The Fisison XS Matrix
	vector< TMatrixT<double> > CaptureXSMatrix; //The Capture XS Matrix
	vector< TMatrixT<double> > n2nXSMatrix;     //The n2N XS Matrix
	double  Flux[16];
	for(int i = 0; i < 16; i++)
	{
		
		double TStep = cycletime/16*i;
		
		if(fUpdateReferenceDBatEachStep)	//Updated At each Step
		{
			IsotopicVector IVStep;
			for(int k = 0; k < (int)findex.size(); k++)
				IVStep += findex.find(k)->second * NMatrix.back()[k][0];
			
			if(i != 0)			// Update closest
				EvolutionDataStep = GetClosest(IVStep, TStep);
			
				//Extraxt XS (if the refecence is not updated will get it from the Evolution DataBase
			TMatrixT<double> SigmaPhiTmp = ExtractXS(EvolutionDataStep, TStep);
			
			for(int j=0; j < (int)findex.size()*3; j++)	// Store it
				SigmaPhi[j][i] = SigmaPhiTmp[j][0];
			
		}
		
		TMatrixT<double> NEvolutionMatrix = TMatrixT<double>(findex.size(),1);
		NEvolutionMatrix = NMatrix.back();
		

		FissionXSMatrix.push_back(GetFissionXsMatrix(EvolutionDataStep, TStep)); //Feel the reaction Matrix
		CaptureXSMatrix.push_back(Getn2nXsMatrix(EvolutionDataStep, TStep)); //Feel the reaction Matrix
		n2nXSMatrix.push_back(GetCaptureXsMatrix(EvolutionDataStep, TStep)); //Feel the reaction Matrix

		double ESigmaN = 0;
		
		for (int j = 0; j < (int)findex.size() ; j++)
			ESigmaN -= FissionXSMatrix[i][j][j]*NEvolutionMatrix[j][0]*1.60217653e-19*FissionEnergy[j][0];
		
			// Update Flux
		Flux[i] = Power/ESigmaN;
		
		
		// ----------------   Evolution
		TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
		
		TMatrixT<double> BatemanReactionMatrix = TMatrixT<double>(findex.size(),findex.size());
		BatemanReactionMatrix = FissionXSMatrix[i];
		BatemanReactionMatrix += CaptureXSMatrix[i];
		BatemanReactionMatrix += n2nXSMatrix[i];
		
		BatemanMatrix = BatemanReactionMatrix;
		BatemanMatrix *= Flux[i];

		BatemanMatrix += fDecayMatrix ;

		double TStepMax = cycletime/16.;
		timevector[i+1] = timevector[i] + TStepMax;
		
		BatemanMatrix *= TStepMax;
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
			BatemanMatrixDLTermN = IdMatrix;
			BatemanMatrixDL = BatemanMatrixDLTermN;
			int j = 1;
			double NormN = 0;
			
			do
			{

				NormN = 0;
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

		NEvolutionMatrix = BatemanMatrixDL * NMatrix.back() ;
		NMatrix.push_back(NEvolutionMatrix);
	}
	
	
	EvolutionData GeneratedDB = EvolutionData(GetLog());
	
	GeneratedDB.SetFlux( new TGraph(16, timevector, Flux)  );
	
	for(int i = 0; i < (int)findex.size(); i++)
	{
		double ZAIQuantity[NMatrix.size()];
		double FissionXS[16];
		double CaptureXS[16];
		double n2nXS[16];
		for(int j = 0; j < (int)NMatrix.size(); j++)
			ZAIQuantity[j] = (NMatrix[j])[i][0];
		
		for(int j = 0; j < 16; j++)
		{
			FissionXS[j]	= FissionXSMatrix[j][i][i];
			CaptureXS[j]	= CaptureXSMatrix[j][i][i];
			n2nXS[j]	= n2nXSMatrix[j][i][i];
		}
		
		GeneratedDB.NucleiInsert(pair<ZAI, TGraph*> (findex.find(i)->second, new TGraph(NMatrix.size(), timevector, ZAIQuantity)));
		if(fUpdateReferenceDBatEachStep)
		{
			GeneratedDB.FissionXSInsert(pair<ZAI, TGraph*> (findex.find(i)->second, new TGraph(16, timevector, FissionXS)));
			GeneratedDB.CaptureXSInsert(pair<ZAI, TGraph*> (findex.find(i)->second, new TGraph(16, timevector, CaptureXS)));
			GeneratedDB.n2nXSInsert(pair<ZAI, TGraph*> (findex.find(i)->second, new TGraph(16, timevector, n2nXS)));
		}
		else
		{
			GeneratedDB.SetFissionXS(EvolutionDataStep.GetFissionXS());
			GeneratedDB.SetCaptureXS(EvolutionDataStep.GetCaptureXS());
			GeneratedDB.Setn2nXS(EvolutionDataStep.Getn2nXS());

		}
	}
	
	GeneratedDB.SetPower(Power );
	GeneratedDB.SetFuelType(fFuelType );
	GeneratedDB.SetReactorType(ReactorType );
	GeneratedDB.SetCycleTime(cycletime);
	
	fDataBankCalculated.insert( pair< IsotopicVector, EvolutionData > ( GeneratedDB.GetIsotopicVectorAt(0.), GeneratedDB) );
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
		
		if( findex_inver.find( (*it).first ) != findex_inver.end() )
		{
			double y;
			y = (*it).second->Eval(TStep);
			
			BatemanMatrix[ findex_inver.find( (*it).first )->second ][findex_inver.find( (*it).first )->second] += -y* 1e-24;
			BatemanMatrix[1][ findex_inver.find( (*it).first )->second] += 2*y* 1e-24;
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
		toAdd.insert(pair<ZAI, double> ( ZAI(95,242,0) , 0.086) );
		toAdd.insert(pair<ZAI, double> ( ZAI(95,242,1) , 0.914) );
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
		if( findex_inver.find( (*it).first ) != findex_inver.end() )
		{
			double y;
			y = (*it).second->Eval(TStep);
			
			BatemanMatrix[findex_inver.find( (*it).first )->second][ findex_inver.find( (*it).first )->second ] += -y* 1e-24 ;
			
			map<ZAI, map<ZAI, double> >::iterator it3 = Capture.find( (*it).first );
			
			if( it3 == Capture.end() )
			{
				map<ZAI, int >::iterator it6 = findex_inver.find( ZAI( (*it).first.Z(), (*it).first.A()+1, (*it).first.I()) );
				
				if( it6 != findex_inver.end() )
				{
					BatemanMatrix[(*it6).second][findex_inver.find( (*it).first )->second] += y* 1e-24  ;
				}
				else
				{
					map<ZAI, map<ZAI, double> >::iterator it4 = fFastDecay.find(  ZAI( (*it).first.Z(), (*it).first.A()+1, (*it).first.I()) );
					
					if( it4 == fFastDecay.end() )
					{
						cout << "Problem in FastDecay for nuclei " << (*it).first.Z() << " " << (*it).first.A()+1 << " " << (*it).first.I() << endl;
						exit(1);
					}
					
					map< ZAI, double >::iterator it5;
					map< ZAI, double > decaylist2 = (*it4).second;
					for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
					{
						it6 = findex_inver.find( (*it5).first );
						if( it6 == findex_inver.end() )
						{
							cout << "Problem in FastDecay for nuclei " << (*it).first.Z() << " " << (*it).first.A() << " " << (*it).first.I() << endl;
							exit(1);
						}
						BatemanMatrix[(*it6).second][findex_inver.find( (*it).first )->second] += y* 1e-24  * (*it5).second;
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
						BatemanMatrix[(*it6).second][findex_inver.find( (*it).first )->second] += y* 1e-24  * (*it4).second ;
					else
					{
						map<ZAI, map<ZAI, double> >::iterator it7 = fFastDecay.find( (*it4).first );
						
						if( it7 == fFastDecay.end() )
						{
							cout << "Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
							exit(1);
						}
						
						map< ZAI, double >::iterator it5;
						map< ZAI, double > decaylist2 = (*it7).second;
						for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
						{
							
							it6 = findex_inver.find( (*it5).first );
							if( it6 == findex_inver.end() )
							{
								cout << "Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
								exit(1);
							}

							BatemanMatrix[(*it6).second][findex_inver.find( (*it).first )->second] += y * 1e-24 * (*it5).second * (*it4).second;
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

		// ----------------  A(n,2n)A-1
	map<ZAI ,TGraph* > n2nXS = EvolutionDataStep.Getn2nXS();
	for(it = n2nXS.begin() ; it != n2nXS.end(); it++)
	{
		if( findex_inver.find( (*it).first ) != findex_inver.end() )
		{
			double y;
			y = (*it).second->Eval(TStep);
			BatemanMatrix[ findex_inver.find( (*it).first )->second ][findex_inver.find( (*it).first )->second] += -y* 1e-24 ;
			
			
			map<ZAI, int>::iterator it3 = findex_inver.find( ZAI( (*it).first.Z(), (*it).first.A()-1, 0) );
			
			if( it3 != findex_inver.end() )
				BatemanMatrix[(*it3).second][findex_inver.find( (*it).first )->second] += y* 1e-24 ;
			else
			{
				
				map<ZAI, map<ZAI, double> >::iterator it4 = fFastDecay.find( ZAI( (*it).first.Z(), (*it).first.A()-1, 0) );
				
				if( it4 == fFastDecay.end() )
				{
					it3 = findex_inver.find( ZAI( -3, -3, -3 ) );
					BatemanMatrix[(*it3).second][findex_inver.find( (*it).first )->second] += y* 1e-24;
				}
				else
				{
					map< ZAI, double >::iterator it5;
					map< ZAI, double > decaylist2 = (*it4).second;
					for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
					{
						
						it3 = findex_inver.find( (*it5).first );
						if( it3 == findex_inver.end() )
						{
							cout << "Problem in FastDecay for nuclei " << (*it4).first.Z() << " " << (*it4).first.A() << " " << (*it4).first.I() << endl;
							exit(1);
						}
						BatemanMatrix[(*it3).second][findex_inver.find( (*it).first )->second] += y* 1e-24 * (*it5).second ;
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
	
		//return EvolutionDataStep.GenerateDBFor(isotopicvector);
	
	ReactorType = EvolutionDataStep.GetReactorType();
	
	for(int i = 0; i < 16; i++)
	{
		
		
		double TStep = cycletime/16*i;
		
		TMatrixT<double> BatemanMatrix = TMatrixT<double>(index.size(),index.size());
		BatemanMatrix = DecayMatrix ;
		
		
		IsotopicVector IVStep;
		for(int k = 0; k < (int)index.size(); k++)
			IVStep += index.find(k)->second * NMatrix.back()[k][0];
		
		if(fUpdateReferenceDBatEachStep && i != 0);		//GetCLosest at the each of evolution step (begining already done...)
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
						//if( (*it3).first.Z() == 90 && (*it3).first.A() == 232) cout << y* 1e-24 *Flux << endl;
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
									//if( (*it6).first.Z() == 92 && (*it6).first.A() == 233) cout << y* 1e-24 *Flux * (*it5).second << endl;
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
