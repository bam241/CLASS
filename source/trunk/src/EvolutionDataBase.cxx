#include "EvolutionDataBase.hxx"

#include "IsotopicVector.hxx"
#include "EvolutiveProduct.hxx"
#include "LogFile.hxx"
#include "Defines.hxx"
#include "StringLine.hxx"

#include <TGraph.h>
#include <TString.h>


#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <map>
#include <cmath>


using namespace std;
//________________________________________________________________________
//
//		EvolutionDataBase
//
//
//
//
//________________________________________________________________________

template<>
void EvolutionDataBase<IsotopicVector>::ReadDataBase();

template<>
EvolutionDataBase<ZAI>::EvolutionDataBase(LogFile* Log, string DB_index_file)
{
	DBGL;
	fLog = Log;
	fDataBaseIndex = DB_index_file;
	
	
	// Warning
	
	cout	<< "!!Info!! !!!EvolutionDataBase<ZAI>!!! A EvolutionData<ZAI> has been define :" << endl;
	cout	<< "\t His index is : \"" << DB_index_file << "\"" << endl << endl;
	
	fLog->fLog 	<< "!!Info!! !!!EvolutionDataBase<ZAI>!!! A EvolutionData<ZAI> has been define :" << endl;
	fLog->fLog	<< "\t His index is : \"" << DB_index_file << "\"" << endl << endl;
	
	DBGL;
}

//________________________________________________________________________
template<>
EvolutionDataBase<ZAI>::~EvolutionDataBase()
{
	DBGL;
	DBGL;
}

template<>
IsotopicVector	EvolutionDataBase<ZAI>::Evolution(const ZAI& zai, double dt)
{
	DBGL;
	IsotopicVector	returnIV;
	
	map<ZAI ,EvolutiveProduct >::iterator it = fEvolutionDataBase.find(zai);
	
	if (it == fEvolutionDataBase.end() )
	{
		ifstream DB_index(fDataBaseIndex.c_str());
		if( !DB_index)
		{
			cout << "!!!EVOLUTIVE DB!!!! Can't open \"" << fDataBaseIndex << "\"" << endl;
			fLog->fLog << "!!!EVOLUTIVE DB!!!! Can't open \"" << fDataBaseIndex << "\"" << endl;
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
				EvolutiveProduct evolutionproduct = EvolutiveProduct(fLog, file_name);
#pragma omp critical(DBupdate)
				{fEvolutionDataBase.insert( pair<ZAI ,EvolutiveProduct >(zai, evolutionproduct) );}
				returnIV = evolutionproduct.GetIsotopicVectorAt(dt);
				zaifind = true;
			}
		}
		
		if(zaifind == false)
		{
			{
				fLog->fLog << "!!Warning!! !!!EVOLUTIVE DB!!! Oups... Can't Find the ZAI : "
				<< zai.Z() << " " << zai.A() << " "	<< zai.I() << "!!! It will be considered as stable !!" << endl;
				EvolutiveProduct evolutionproduct = EvolutiveProduct(fLog," " , true, zai);
				{fEvolutionDataBase.insert( pair<ZAI, EvolutiveProduct >(zai, evolutionproduct) );}
				returnIV = evolutionproduct.GetIsotopicVectorAt(dt);
				
			}
		}
		
		
	}
	else	returnIV = (*it).second.GetIsotopicVectorAt(dt);
	
	return returnIV;
}

template<>
bool EvolutionDataBase<ZAI>::IsDefine(const ZAI& zai) const
{
	DBGL;
	map<ZAI ,EvolutiveProduct > evolutiondb = (*this).GetEvolutionDataBase();
	if (evolutiondb.find(zai) != evolutiondb.end())
		return true;
	else
		return false;
	
}
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________

template<>
EvolutionDataBase<IsotopicVector>::~EvolutionDataBase()
{
	DBGL;
	DBGL;
}

template<>
EvolutionDataBase<IsotopicVector>::EvolutionDataBase(LogFile* Log, string DB_index_file)
{
	DBGL;
	fLog = Log;
	fDataBaseIndex = DB_index_file;
	
	ReadDataBase();
	
	
	// Warning
	cout	<< "!!Info!! !!!EvolutionDataBase<IsotopicVector>!!! A EvolutionData<ZAI> has been define :" << endl;
	cout	<< "\t His index is : \"" << DB_index_file << "\"" << endl;
	cout	<< "\t " << fEvolutionDataBase.size() << " EvolutiveProduct have been read."<< endl << endl;
	
	fLog->fLog 	<< "!!Info!! !!!EvolutionDataBase<IsotopicVector>!!! A EvolutionData<ZAI> has been define :" << endl;
	fLog->fLog	<< "\t His index is : \"" << DB_index_file << "\"" << endl;
	fLog->fLog	<< "\t " << fEvolutionDataBase.size() << " EvolutiveProduct have been read."<< endl << endl;
	
	DBGL;
}

template<>
void EvolutionDataBase<IsotopicVector>::ReadDataBase()
{
	DBGL;
	ifstream DataDB(fDataBaseIndex.c_str());							// Open the File
	if(!DataDB)
	{
		cout << "!!Warning!! !!!EvolutionDataBase!!! \n Can't open \"" << fDataBaseIndex << "\"\n" << endl;
		fLog->fLog << "!!Warning!! !!!EvolutionDataBase!!! \n Can't open \"" << fDataBaseIndex << "\"\n" << endl;
	}
	vector<double> vTime;
	vector<double> vTimeErr;
	
	string line;
	int start = 0;
	
	
	
	// First Get Fuel Type
	getline(DataDB, line);
	if( StringLine::NextWord(line, start, ' ') != "TYPE")
	{
		cout << "!!Bad Trouble!! !!!EvolutionDataBase!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the type of the DataBase"<< endl;
		fLog->fLog << "!!Bad Trouble!! !!!EvolutionDataBase!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the type of the DataBase"<< endl;
		exit (1);
	}
	fFuelType = StringLine::NextWord(line, start, ' ');
	// First Get Fuel Parameter
	getline(DataDB, line);
	start = 0;
	if( StringLine::NextWord(line, start, ' ') != "PARAM")
	{
		cout << "!!Bad Trouble!! !!!EvolutionDataBase!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the Parameter of the DataBase"<< endl;
		fLog->fLog << "!!Bad Trouble!! !!!EvolutionDataBase!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the Parameter of the DataBase"<< endl;
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
			EvolutiveProduct* evolutionproduct = new EvolutiveProduct(fLog, line);
			IsotopicVector ivtmp  = evolutionproduct->GetIsotopicVectorAt(0.).GetActinidesComposition();
			fEvolutionDataBase.insert( pair<IsotopicVector, EvolutiveProduct >(ivtmp , (*evolutionproduct) ));
		}
	}
	DBGL;
}
//________________________________________________________________________



template<>
map<double, EvolutiveProduct> EvolutionDataBase<IsotopicVector>::GetDistancesTo(IsotopicVector isotopicvector, double t) const
{
	DBGL;
	map<double, EvolutiveProduct> distances;
	
	map<IsotopicVector, EvolutiveProduct > evolutiondb = fEvolutionDataBase;
	
	map<IsotopicVector, EvolutiveProduct >::iterator it;
	for( it = evolutiondb.begin(); it != evolutiondb.end(); it++ )
	{
		pair<map<double, EvolutiveProduct>::iterator, bool> IResult;
		double D = Distance(isotopicvector.GetActinidesComposition(), (*it).second.GetIsotopicVectorAt(t).GetActinidesComposition()/ Norme( (*it).second.GetIsotopicVectorAt(t).GetActinidesComposition() )*Norme(isotopicvector.GetActinidesComposition()) );
		//		double D = RelativDistance(isotopicvector, (*it).second.GetIsotopicVectorAt(t)/ Norme( (*it).second.GetIsotopicVectorAt(t) )*Norme(isotopicvector) );
		IResult = distances.insert( pair<double, EvolutiveProduct>( D , (*it).second ) );
	}
	
	return distances;
	DBGL;
}

template<>
EvolutiveProduct EvolutionDataBase<IsotopicVector>::GenerateDB(IsotopicVector isotopicvector, double cycletime, double Power)
{
	DBGL;
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
	
	for(int i = 0; i < 16; i++)
	{
		
		
		double TStep = cycletime/16*i;
		
		TMatrixT<double> BatemanMatrix = TMatrixT<double>(index.size(),index.size());
		BatemanMatrix = DecayMatrix ;
		
		
		IsotopicVector IVStep;
		for(int k = 0; k < (int)index.size(); k++)
			IVStep += index.find(k)->second * NMatrix.back()[k][0];
		
		
		EvolutiveProduct evolutiveproductStep = GetDistancesTo(IVStep, TStep).begin()->second;
		
		double NormFactor = 1;
		{
			IsotopicVector WantedHMIV = 	  isotopicvector.GetSpeciesComposition(90)
			+ isotopicvector.GetSpeciesComposition(92)
			+ isotopicvector.GetSpeciesComposition(93)
			+ isotopicvector.GetSpeciesComposition(94)
			+ isotopicvector.GetSpeciesComposition(95)
			+ isotopicvector.GetSpeciesComposition(96);
			
			IsotopicVector DBHMIV =   evolutiveproductStep.GetIsotopicVectorAt(0).GetSpeciesComposition(90)
			+ evolutiveproductStep.GetIsotopicVectorAt(0).GetSpeciesComposition(92)
			+ evolutiveproductStep.GetIsotopicVectorAt(0).GetSpeciesComposition(93)
			+ evolutiveproductStep.GetIsotopicVectorAt(0).GetSpeciesComposition(94)
			+ evolutiveproductStep.GetIsotopicVectorAt(0).GetSpeciesComposition(95)
			+ evolutiveproductStep.GetIsotopicVectorAt(0).GetSpeciesComposition(96);
			
			NormFactor = Norme(WantedHMIV)/ Norme(DBHMIV);
		}
		
		
		double Flux = evolutiveproductStep.GetFlux()->Eval(TStep)*Power/(evolutiveproductStep.GetPower()*NormFactor);
		
		map<ZAI ,TGraph* >::iterator it;
		// ----------------  A(n,.) X+Y
		
		map<ZAI ,TGraph* > FissionXS = evolutiveproductStep.GetFissionXS();
		
		for(it = FissionXS.begin() ; it != FissionXS.end(); it++)
		{
			
			if( index_inver.find( (*it).first ) != index_inver.end() )
			{
				double y;
				y = (*it).second->Eval(TStep);
				
				BatemanMatrix[ index_inver.find( (*it).first )->second ][index_inver.find( (*it).first )->second] += -y* 1e-24 *Flux;
				BatemanMatrix[1][ index_inver.find( (*it).first )->second] += 2*y* 1e-24 *Flux;
			}
			
		}
		
		// ----------------  A(n,.)A+1
		map<ZAI ,TGraph* > CaptureXS = evolutiveproductStep.GetCaptureXS();
		for(it = CaptureXS.begin(); it != CaptureXS.end(); it++)
		{
			
			if( index_inver.find( (*it).first ) != index_inver.end() )
			{
				double y;
				y = (*it).second->Eval(TStep);
				
				BatemanMatrix[index_inver.find( (*it).first )->second][ index_inver.find( (*it).first )->second ] += -y* 1e-24 *Flux;
				
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
		map<ZAI ,TGraph* > n2nXS = evolutiveproductStep.Getn2nXS();
		for(it = n2nXS.begin() ; it != n2nXS.end(); it++)
		{
			if( index_inver.find( (*it).first ) != index_inver.end() )
			{
				double y;
				y = (*it).second->Eval(TStep);
				BatemanMatrix[ index_inver.find( (*it).first )->second ][index_inver.find( (*it).first )->second] += -y* 1e-24 *Flux;
				
				
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
			} while ( NormN != 0);
		}
		
		NEvolutionMatrix = BatemanMatrixDL * NMatrix.back() ;
		NMatrix.push_back(NEvolutionMatrix);
	}
	
	
	EvolutiveProduct GeneratedDB = EvolutiveProduct(fLog);
	
	for(int i = 0; i < (int)index.size(); i++)
	{
		double ZAIQuantity[NMatrix.size()];
		for(int j = 0; j < (int)NMatrix.size(); j++)
			ZAIQuantity[j] = (NMatrix[j])[i][0];
		
		GeneratedDB.Insert(pair<ZAI, TGraph*> (index.find(i)->second, new TGraph(NMatrix.size(), timevector, ZAIQuantity) ) );
	}
	
	GeneratedDB.SetPower(Power );
	GeneratedDB.SetFuelType(fFuelType );
		//GeneratedDB.SetReactorType(fReactorType );
		//GeneratedDB.SetHMMass(fHMMass*NormFactor );
	
	return GeneratedDB;
	DBGL;
}






