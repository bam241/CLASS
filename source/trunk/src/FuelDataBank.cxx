//
//  FuelDataBank.cxx
//  CLASSSource
//
//  Created by BaM on 04/05/2014.
//  Copyright (c) 2014 BaM. All rights reserved.
//

#include "FuelDataBank.hxx"

#include "IsotopicVector.hxx"
#include "CLASSHeaders.hxx"
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
//		FuelDataBank
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
			XS += DB.GetXSForAt(0., (*it).first, i);
		
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



//________________________________________________________________________
FuelDataBank::FuelDataBank():DynamicalSystem()
{
	fTheNucleiVector = 0;
	fTheMatrix = 0;
	
	fWeightedDistance = false;
	fEvolutionDataInterpolation = false;
	
	
	
	fOldReadMethod = true;
	fUseRK4EvolutionMethod = true;
	fDistanceType = 0;
	fShorstestHalflife = 3600.*24*2.;
	fZAIThreshold = 90;
	
	
	fDataDirectoryName = getenv("CLASS_PATH");
	fDataDirectoryName += "/source/data/";
	fDataFileName = "chart.JEF3T";
	
	SetForbidNegativeValue();
	
}


FuelDataBank::FuelDataBank(LogFile* Log, string DB_index_file, bool setlog, bool olfreadmethod):DynamicalSystem()
{
	SetLog(Log);
	fWeightedDistance = false;
	fEvolutionDataInterpolation = false;
	
	
	fTheNucleiVector = 0;
	fTheMatrix = 0;
	
	fDataDirectoryName = getenv("CLASS_PATH");
	fDataDirectoryName += "/source/data/";
	fDataFileName = "chart.JEF3T";
	
	fDataBaseIndex = DB_index_file;
	fOldReadMethod = olfreadmethod;
	fUseRK4EvolutionMethod = true;
	
	fDistanceType = 0;
	
	fShorstestHalflife = 3600.*24*2.;
	fZAIThreshold = 90;
	
	ReadDataBase();
	
	SetForbidNegativeValue();
	
	if(IsLog())
	{
		// Warning
		cout	<< "!!INFO!! !!!FuelDataBank!!! A EvolutionData has been define :" << endl;
		cout	<< "\t His index is : \"" << DB_index_file << "\"" << endl;
		cout	<< "\t " << fFuelDataBank.size() << " EvolutionData have been read."<< endl << endl;
		
		GetLog()->fLog 	<< "!!INFO!! !!!FuelDataBank!!! A EvolutionData has been define :" << endl;
		GetLog()->fLog	<< "\t His index is : \"" << DB_index_file << "\"" << endl;
		GetLog()->fLog	<< "\t " << fFuelDataBank.size() << " EvolutionData have been read."<< endl << endl;
	}
	
	
}

FuelDataBank::~FuelDataBank()
{
	if(fTheMatrix)
	{
		for(int i= 0; i<fNVar; i++)
			delete [] fTheMatrix[i];
		delete [] fTheMatrix;
	}
	if(fTheNucleiVector)
	{
		delete fTheNucleiVector;
	}
	
	map<IsotopicVector ,EvolutionData >::iterator it_del;
	for( it_del = fFuelDataBank.begin(); it_del != fFuelDataBank.end(); it_del++)
		(*it_del).second.DeleteEvolutionData();
	fFuelDataBank.clear();
	
	for( it_del = fFuelDataBankCalculated.begin(); it_del != fFuelDataBankCalculated.end(); it_del++)
		(*it_del).second.DeleteEvolutionData();
	fFuelDataBankCalculated.clear();
	
	fFissionEnergy.clear();
	fFastDecay.clear();
	fSpontaneusYield.clear();
	fReactionYield.clear();
	findex_inver.clear();
	findex.clear();
	fDecayMatrix.Clear();
	
}

void FuelDataBank::Clear()
{
	if(fTheMatrix)
	{
		for(int i= 0; i<fNVar; i++)
			delete [] fTheMatrix[i];
		delete [] fTheMatrix;
	}
	if(fTheNucleiVector)
	{
		delete fTheNucleiVector;
	}
	
	fFuelDataBank.clear();
	fFuelDataBankCalculated.clear();
	fFissionEnergy.clear();
	fFastDecay.clear();
	fSpontaneusYield.clear();
	fReactionYield.clear();
	findex_inver.clear();
	findex.clear();
	
	
	fTheNucleiVector = 0;
	fTheMatrix = 0;
	fNVar = 0;
	
	fOldReadMethod = true;
	fUseRK4EvolutionMethod = true;
	fDistanceType = 0;
	fShorstestHalflife = 3600.*24*2.;
	fZAIThreshold = 90;
	
	
	fDataDirectoryName = getenv("CLASS_PATH");
	fDataDirectoryName += "/source/data/";
	fDataFileName = "chart.JEF3T";
	
	SetForbidNegativeValue();
	
	
	fDecayMatrix.Clear();
	
	
	
}


//________________________________________________________________________
void FuelDataBank::ReadDataBase()
{
	
	if(fFuelDataBank.size() != 0)
	{
		map<IsotopicVector ,EvolutionData >::iterator it_del;
		for( it_del = fFuelDataBank.begin(); it_del != fFuelDataBank.end(); it_del++)
			(*it_del).second.DeleteEvolutionData();
		fFuelDataBank.clear();
	}
	
	if(fFuelDataBankCalculated.size() != 0)
	{
		map<IsotopicVector ,EvolutionData >::iterator it_del;
		for( it_del = fFuelDataBankCalculated.begin(); it_del != fFuelDataBankCalculated.end(); it_del++)
			(*it_del).second.DeleteEvolutionData();
		fFuelDataBankCalculated.clear();
	}
	
	
	
	
	ifstream DataDB(fDataBaseIndex.c_str());							// Open the File
	if(!DataDB)
	{
		cout << "!!Warning!! !!!FuelDataBank!!! \n Can't open \"" << fDataBaseIndex << "\"\n" << endl;
		GetLog()->fLog << "!!Warning!! !!!FuelDataBank!!! \n Can't open \"" << fDataBaseIndex << "\"\n" << endl;
	}
	vector<double> vTime;
	vector<double> vTimeErr;
	
	string line;
	int start = 0;
	
	
	
	// First Get Fuel Type
	getline(DataDB, line);
	if( StringLine::NextWord(line, start, ' ') != "TYPE")
	{
		cout << "!!Bad Trouble!! !!!FuelDataBank!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the type of the DataBase"<< endl;
		GetLog()->fLog << "!!Bad Trouble!! !!!FuelDataBank!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the type of the DataBase"<< endl;
		exit (1);
	}
	fFuelType = StringLine::NextWord(line, start, ' ');
	// First Get Fuel Parameter
	getline(DataDB, line);
	start = 0;
	if( StringLine::NextWord(line, start, ' ') != "PARAM")
	{
		cout << "!!Bad Trouble!! !!!FuelDataBank!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the Parameter of the DataBase"<< endl;
		GetLog()->fLog << "!!Bad Trouble!! !!!FuelDataBank!!! Bad Database file : " <<  fDataBaseIndex << " Can't find the Parameter of the DataBase"<< endl;
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
			fFuelDataBank.insert( pair<IsotopicVector, EvolutionData >(ivtmp , (*evolutionproduct) ));
		}
	}
	
}
//________________________________________________________________________




//________________________________________________________________________
//________________________________________________________________________
/*				Physics				*/
//________________________________________________________________________
//________________________________________________________________________



//________________________________________________________________________
/*				Decay Stuff			*/
//________________________________________________________________________
string FuelDataBank::GetDecay(string DecayModes, double &BR,int &Iso, int &StartPos)
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

//________________________________________________________________________
void FuelDataBank::BuildDecayMatrix()
{
	fDecayMatrix.Clear();
	
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
		cout << "!!Warning!! !!!FuelDataBank!!! \n Can't open \"" << DataFullPathName << "\"\n" << endl;
		GetLog()->fLog << "!!Warning!! !!!FuelDataBank!!! \n Can't open \"" << DataFullPathName<< "\"\n" << endl;
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
		if(Z >= fZAIThreshold )
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
			IsotopicVector DaughtersMap;
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
				
				if (decay_name == "s")	{DM=0;daughter_N=(A-Z);	daughter_Z=Z;BR=1;}
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
					if( daughter_Z < fZAIThreshold && daughter_Z!=-2 )
						daughter_A = daughter_Z = Iso = -3;
					// not spontaneous fission
					ZAI DaughterZAI = ZAI(daughter_Z,daughter_A,Iso);
					
					if((BR>1e-10) && (!stable))
					{
						if(DM <= 15)
						{
							DaughtersMap += BR * DaughterZAI;
							branch_test+=BR;
						}
						else if( DM <= 18)
						{
							if(fSpontaneusYield.size() == 0 || fReactionYield.size() == 0 || DM == 17 || DM == 18)
							{
								DaughtersMap += 2*BR * ZAI(-2,-2,-2);
								
								branch_test_f += BR;
							}
							else
							{
								
								map<ZAI, IsotopicVector>::iterator it_yield = fSpontaneusYield.find(ParentZAI);
								if(it_yield != fSpontaneusYield.end())
								{
									DaughtersMap += (BR* (*it_yield).second );
									branch_test_f += BR* (*it_yield).second.GetSumOfAll() / 2.;
									
								}
								else
								{
									DaughtersMap += 2*BR * ZAI(-2,-2,-2);
									branch_test_f += BR;
								}
							}
						}
						
					}
					
				}
				if (DM !=0)
					stable = false;
				// End of While loop
			}
			
			double btest = fabs(branch_test + branch_test_f-1.0);
			if ( btest > 1e-8 && !stable )
				if(branch_test+branch_test_f > 0)
					DaughtersMap = DaughtersMap /(branch_test+branch_test_f);
			
			
			
			if (HalfLife < fShorstestHalflife && !stable)
				fFastDecay.insert( pair< ZAI, map<ZAI, double> > ( ParentZAI, DaughtersMap.GetIsotopicQuantity() ) );
			else if (stable)
			{
				IsotopicVector StableIV = ParentZAI *1;
				ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >
						(ParentZAI, pair<double, map< ZAI, double > >
						 ( 1e36, StableIV.GetIsotopicQuantity()) ) );
			}
			else
				ZAIDecay.insert( pair< ZAI, pair<double, map< ZAI, double > > >
						(ParentZAI, pair<double, map< ZAI, double > >
						 ( HalfLife, DaughtersMap.GetIsotopicQuantity()) ) );
			
			
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
					if( FD2_it != FastDecayCopy.end() )
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
					else
					{
						(*FD_it_Origin).second.erase( (*BR_it).first );
						pair< map<ZAI, double>::iterator, bool> IResult;
						IResult = (*FD_it_Origin).second.insert( pair<ZAI, double> ( ZAI(-3,-3,-3), (*BR_it).second) );
						if( !IResult.second)
							(*IResult.first).second += (*BR_it).second ;
					}
					
					
					
				}
			}
			
		}
		
		
		FastDecayValidation = true;
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
	
}

//________________________________________________________________________
/*				Fission Stuff			*/
//________________________________________________________________________
void FuelDataBank::SetFissionEnergy(ZAI zai, double E)
{
	pair<map<ZAI, double>::iterator, bool> IResult;
	IResult = fFissionEnergy.insert( pair<ZAI ,double>(zai, E));
	if(IResult.second == false)
		IResult.first->second = E;
	
}

//________________________________________________________________________
void FuelDataBank::SetFissionEnergy(string FissionEnergyFile)
{
	ifstream FissionFile(FissionEnergyFile.c_str());	// Open the File
	if(!FissionFile)				//check if file is correctly open
	{
		cout << "!!Warning!! !!!FuelDataBank!!! \n Can't open \"" << FissionFile << "\"\n" << endl;
		GetLog()->fLog << "!!Warning!! !!!FuelDataBank!!! \n Can't open \"" << FissionFile << "\"\n" << endl;
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
map< ZAI,IsotopicVector > FuelDataBank::ReadFPYield(string Yield)
{
	IsotopicVector EmptyIV;
	map< ZAI,IsotopicVector >  Yield_map;
	
	ifstream infile(Yield.c_str());
	if(!infile)
	{
		cout << "!!Warning!! !!!FuelDataBank!!! \n Can't open \"" << Yield << "\"\n" << endl;
		GetLog()->fLog << "!!Warning!! !!!FuelDataBank!!! \n Can't open \"" << Yield<< "\"\n" << endl;
	}
	
	
	string line;
	int start = 0;
	
	getline(infile, line);
	
	do
	{
		int Z = atof(StringLine::NextWord(line, start, ' ').c_str());
		int A = atof(StringLine::NextWord(line, start, ' ').c_str());
		int I = 0;
		
		//		if(Z!=0 && A!=0)
		{
			pair<map<ZAI, IsotopicVector>::iterator, bool> IResult;
			IResult = Yield_map.insert(pair<ZAI,IsotopicVector>(ZAI(Z,A,I),EmptyIV) );
			if (IResult.second == false)
			{
				cout << "!!Error!! !!!FuelDataBank!!! Many accurance of ZAI " << Z << " " << A;
				cout << " in " << Yield << " file!! Please Check it !!!" << endl;
				exit(1);
				
			}
		}
	}while(start < (int)line.size()-1);
	
	do
	{
		start = 0;
		
		getline(infile, line);
		int Z = atof(StringLine::NextWord(line, start, ' ').c_str());
		int A = atof(StringLine::NextWord(line, start, ' ').c_str());
		int I = atof(StringLine::NextWord(line, start, ' ').c_str());
		map<ZAI, IsotopicVector>::iterator it = Yield_map.begin();
		do
		{
			if (it == Yield_map.end())
			{
				cout << "!!Error!! !!!FuelDataBank!!! Many accurance of the PF " << Z << " " << A;
				cout << " in " << Yield << " file!! Please Check it";
				cout << "(Number of yield does not match the number of ZAI that fission !!!" << endl;
				exit(1);
				
			}
			
			double Yield_values = atof(StringLine::NextWord(line, start, ' ').c_str());
			(*it).second +=  Yield_values * ZAI(Z,A,I);
			
			it++;
		}while(start < (int)line.size()-1);
		
		
		
		
	} while (!infile.eof());
	return Yield_map;
}

//________________________________________________________________________
void FuelDataBank::LoadFPYield(string SponfaneusYield, string ReactionYield)
{
	
	fSpontaneusYield = ReadFPYield(SponfaneusYield);
	fReactionYield = ReadFPYield(ReactionYield);
	fZAIThreshold = 0;
}



//________________________________________________________________________
/*				Reaction Stuff			*/
//________________________________________________________________________
TMatrixT<double> FuelDataBank::GetFissionXsMatrix(EvolutionData EvolutionDataStep,double TStep)
{
	
	map<ZAI ,TGraph* >::iterator it;
	TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
	for(int i = 0; i < (int)findex.size(); i++)
		for(int j = 0; j < (int)findex.size(); j++)
			BatemanMatrix[i][j] = 0;
	
	// ----------------  A(n,.) X+Y
	
	map<ZAI ,TGraph* > FissionXS = EvolutionDataStep.GetFissionXS();
	
	for(it = FissionXS.begin() ; it != FissionXS.end(); it++)
	{
		map<ZAI, int>::iterator findex_inver_it = findex_inver.find( (*it).first );
		if( findex_inver_it != findex_inver.end() )
		{
			double y = (*it).second->Eval(TStep);
			BatemanMatrix[ findex_inver_it->second ][ findex_inver_it->second ] += -y* 1e-24;
			
			if(fSpontaneusYield.size() == 0 || fReactionYield.size() == 0)
				BatemanMatrix[1][ findex_inver_it->second ] += 2*y* 1e-24;
			else
			{
				map<ZAI, IsotopicVector>::iterator it_yield = fReactionYield.find( (*it).first );
				
				if( it_yield != fReactionYield.end())
				{
					map<ZAI ,double>::iterator it_IVQ;
					map<ZAI ,double> IVQ = (*it_yield).second.GetIsotopicQuantity();
					
					for( it_IVQ = IVQ.begin(); it_IVQ != IVQ.end(); it_IVQ++ )
					{
						map<ZAI, int>::iterator findex_it_PF = findex_inver.find( (*it_IVQ).first );
						
						if(findex_it_PF != findex_inver.end() )
							BatemanMatrix[(*findex_it_PF).second][ (*findex_inver_it).second ] += (*it_IVQ).second*y* 1e-24;
						else
						{
							map<ZAI, map<ZAI, double> >::iterator it_FD = fFastDecay.find( (*it_IVQ).first);
							
							if( it_FD == fFastDecay.end() )
							{
								BatemanMatrix[1][ (*findex_inver_it).second ] += (*it_IVQ).second * y * 1e-24  ;
							}
							else
							{
								
								map< ZAI, double >::iterator it5;
								map< ZAI, double > decaylist2 = (*it_FD).second;
								for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
								{
									findex_it_PF = findex_inver.find( (*it5).first );
									if( findex_it_PF == findex_inver.end() )
										BatemanMatrix[0][findex_inver_it->second] +=
										(*it_IVQ).second * y * 1e-24 * (*it5).second;
									else
										BatemanMatrix[(*findex_it_PF).second][findex_inver_it->second]+=
										(*it_IVQ).second * y * 1e-24 * (*it5).second;
								}
							}
							
						}
						
					}
				}
				else
					BatemanMatrix[1][ findex_inver_it->second ] += 2*y* 1e-24;
				
			}
		}
		
	}
	
	return BatemanMatrix;
	
}

//________________________________________________________________________
TMatrixT<double> FuelDataBank::GetCaptureXsMatrix(EvolutionData EvolutionDataStep,double TStep)
{
	
	
	map<ZAI ,TGraph* >::iterator it;
	TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
	for(int i = 0; i < (int)findex.size(); i++)
		for(int j = 0; j < (int)findex.size(); j++)
			BatemanMatrix[i][j] = 0;
	
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
TMatrixT<double> FuelDataBank::Getn2nXsMatrix(EvolutionData EvolutionDataStep,double TStep)
{
	
	
	map<ZAI ,TGraph* >::iterator it;
	TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
	for(int i = 0; i < (int)findex.size(); i++)
		for(int j = 0; j < (int)findex.size(); j++)
			BatemanMatrix[i][j] = 0;
	
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



//________________________________________________________________________
//________________________________________________________________________
/*			Distance Calculation			*/
//________________________________________________________________________
//________________________________________________________________________
map<double, EvolutionData> FuelDataBank::GetDistancesTo(IsotopicVector isotopicvector, double t) const
{
	
	map<double, EvolutionData> distances;
	
	map<IsotopicVector, EvolutionData > evolutiondb = fFuelDataBank;
	
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

//________________________________________________________________________
EvolutionData FuelDataBank::GetClosest(IsotopicVector isotopicvector, double t) const
{
	
	map<IsotopicVector, EvolutionData > evolutiondb = fFuelDataBank;
	double distance = 0;
	
	map<IsotopicVector, EvolutionData >::iterator it_close = evolutiondb.begin();
	
	
	map<IsotopicVector, EvolutionData >::iterator it;
	
	
	if(fWeightedDistance)
	{
		distance = Distance(isotopicvector.GetActinidesComposition()
				    * evolutiondb.begin()->second.GetIsotopicVectorAt(t).GetActinidesComposition().GetSumOfAll()
				    / isotopicvector.GetActinidesComposition().GetSumOfAll(),
				    evolutiondb.begin()->second);
		
		
		for( it = evolutiondb.begin(); it != evolutiondb.end(); it++ )
		{
			double D = 0;
			D = Distance(isotopicvector.GetActinidesComposition()
				     * (*it).second.GetIsotopicVectorAt(t).GetActinidesComposition().GetSumOfAll()
				     / isotopicvector.GetActinidesComposition().GetSumOfAll(),
				     (*it).second);
			
			if (D< distance)
			{
				distance = D;
				it_close = it;
			}
		}
		
		return (*it_close).second;
	}
	else if (fEvolutionDataInterpolation)
	{
		
		map<double, EvolutionData> distance_map = GetDistancesTo(isotopicvector, t);
		map<double, EvolutionData>::iterator it_distance;
		int NClose = 64;
		int Nstep = 0;
		EvolutionData EvolInterpolate;
		double SumOfDistance = 0;
		for( it_distance = distance_map.begin(); it_distance != distance_map.end(); it_distance++)
		{
			
			if((*it_distance).first == 0 )
			{
				EvolInterpolate = Multiply(1,(*it_distance).second);
				return EvolInterpolate;
			}
			if(Nstep == 0)
				EvolInterpolate = Multiply(1./(*it_distance).first, (*it_distance).second);
			else
			{
				
				EvolutionData Evoltmp = EvolInterpolate;
				EvolutionData Evoltmp2 = Multiply(1./(*it_distance).first, (*it_distance).second);
				
				EvolInterpolate = Sum(Evoltmp,  Evoltmp2);
				Evoltmp.DeleteEvolutionData();
				Evoltmp2.DeleteEvolutionData();
				
				
			}
			
			SumOfDistance += 1./(*it_distance).first;
			Nstep++;
			if(Nstep == NClose) break;
			
		}
		
		EvolutionData Evoltmp = EvolInterpolate;
		EvolInterpolate.Clear();
		
		EvolInterpolate = Multiply(1/SumOfDistance, Evoltmp);
		
		Evoltmp.DeleteEvolutionData();
		return EvolInterpolate;
		
	}
	else
	{
		distance = Distance(isotopicvector.GetActinidesComposition(),
				    evolutiondb.begin()->second.GetIsotopicVectorAt(t).GetActinidesComposition()
				    / evolutiondb.begin()->second.GetIsotopicVectorAt(t).GetActinidesComposition().GetSumOfAll()
				    * isotopicvector.GetActinidesComposition().GetSumOfAll(),
				    fDistanceType, fDistanceParameter);
		for( it = evolutiondb.begin(); it != evolutiondb.end(); it++ )
		{
			
			double D = 0;
			
			
			D = Distance(isotopicvector.GetActinidesComposition(),
				     (*it).second.GetIsotopicVectorAt(t).GetActinidesComposition()
				     / (*it).second.GetIsotopicVectorAt(t).GetActinidesComposition().GetSumOfAll()
				     * isotopicvector.GetActinidesComposition().GetSumOfAll(),
				     fDistanceType, fDistanceParameter);
			
			if (D< distance)
			{
				distance = D;
				it_close = it;
			}
		}
		return (*it_close).second;
		
	}
	
	
}

//________________________________________________________________________

void FuelDataBank::CalculateDistanceParameter()
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
	map<IsotopicVector ,EvolutionData > FuelDataBank = (*this).GetFuelDataBank();
	int NevolutionDatainFuelDataBank = 0;
	
	for( it = FuelDataBank.begin(); it != FuelDataBank.end(); it++ )
	{
		NevolutionDatainFuelDataBank++;
		map<ZAI ,double>::iterator itit;
		map<ZAI ,double> isovector=(*it).first.GetIsotopicQuantity();
		for(itit=isovector.begin(); itit != isovector.end(); itit++) //Boucle sur ZAI
		{
			double TmpXS=0;
			
			for( int i=1; i<4; i++ ) //Loop on Reactions 1==fission, 2==capture, 3==n2n
				TmpXS+=	(*it).second.GetXSForAt(0, (*itit).first, i);
			
			fDistanceParameter.Add((*itit).first,TmpXS);
		}
		
		
	}
	fDistanceParameter.Multiply( (double)1.0/NevolutionDatainFuelDataBank );
	
	if(GetLog()){
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
	
	
}

//________________________________________________________________________
void FuelDataBank::SetDistanceParameter(IsotopicVector DistanceParameter)
{
	
	fDistanceParameter = DistanceParameter;
	
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
void FuelDataBank::SetDistanceType(int DistanceType)
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



//________________________________________________________________________
//________________________________________________________________________
/*				Evolution 			*/
//________________________________________________________________________
//________________________________________________________________________

//________________________________________________________________________
/*				RK4 Stuff			*/
//________________________________________________________________________
//________________________________________________________________________
void FuelDataBank::ResetTheMatrix()
{
	
	if(fTheMatrix)
	{
		for(int i= 0; i<fNVar; i++)
			delete [] fTheMatrix[i];
		delete [] fTheMatrix;
	}
	fTheMatrix = 0;
}

void FuelDataBank::SetTheMatrixToZero()
{
	ResetTheMatrix();
	
	fNVar = findex.size();
	fTheMatrix = new double*[fNVar];
	
#pragma omp parallel for
	for(int i= 0; i < fNVar; i++)
		fTheMatrix[i] = new double[fNVar];
	
	for(int i = 0; i < fNVar; i++)
		for(int k = 0; k < fNVar; k++)
		{
			fTheMatrix[i][k]=0.0;
		}
	
}

//________________________________________________________________________
void FuelDataBank::ResetTheNucleiVector()
{
	if(fTheNucleiVector)
		delete [] fTheNucleiVector;
	fTheNucleiVector = 0;
}

//________________________________________________________________________
void FuelDataBank::SetTheNucleiVectorToZero()
{
	ResetTheNucleiVector();
	fTheNucleiVector = new double[fNVar];
	
#pragma omp parallel for
	for(int i = 0; i < fNVar; i++)
		fTheNucleiVector[i]=0.0;
	
}

//________________________________________________________________________
void FuelDataBank::BuildEqns(double t, double *N, double *dNdt)
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

//________________________________________________________________________
void FuelDataBank::SetTheMatrix(TMatrixT<double> BatemanMatrix)
{
	for (int k = 0; k < (int)fNVar; k++)
		for (int l = 0; l < (int)findex_inver.size(); l++)
			fTheMatrix[l][k] = BatemanMatrix[l][k];
}

//________________________________________________________________________
TMatrixT<double> FuelDataBank::GetTheMatrix()
{
	TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
	for (int k = 0; k < (int)fNVar; k++)
		for (int l = 0; l < (int)findex_inver.size(); l++)
			BatemanMatrix[l][k] = fTheMatrix[l][k];
	
	return BatemanMatrix;
}

//________________________________________________________________________
void FuelDataBank::SetTheNucleiVector(TMatrixT<double> NEvolutionMatrix)
{
	for (int k = 0; k < (int)fNVar; k++)
		fTheNucleiVector[k] = NEvolutionMatrix[k][0];
}

//________________________________________________________________________
TMatrixT<double> FuelDataBank::GetTheNucleiVector()
{
	TMatrixT<double> NEvolutionMatrix = TMatrixT<double>(findex.size(),1);
	for (int k = 0; k < (int)fNVar; k++)
		NEvolutionMatrix[k][0] = fTheNucleiVector[k];
	
	return NEvolutionMatrix;
}


//________________________________________________________________________
/*			Evolution Calculation			*/
//________________________________________________________________________
EvolutionData FuelDataBank::GenerateEvolutionData(IsotopicVector isotopicvector, double cycletime, double Power)
{
	
	if(fFastDecay.size() == 0)
	{
		BuildDecayMatrix();
		fNVar = findex_inver.size();
	}
	
	SetTheMatrixToZero();
	SetTheNucleiVectorToZero();
	
	string ReactorType;
	
	
	vector< TMatrixT<double> > NMatrix ;//  TMatrixT<double>(decayindex.size(),1))
	{	// Filling the t=0 State;
		map<ZAI, double > isotopicquantity = isotopicvector.GetIsotopicQuantity();
		TMatrixT<double>  N_0Matrix =  TMatrixT<double>( findex.size(),1) ;
		for(int i = 0; i < (int)findex.size(); i++)
			N_0Matrix[i] = 0;
		
		map<ZAI, double >::iterator it ;
		for(int i = 0; i < (int)findex.size(); i++)
			N_0Matrix[i] = 0;
		
		for(it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		{
			/// Need TO change with FP managment
			map<ZAI, int >::iterator it2;
			
			if( (*it).first.Z() < fZAIThreshold )
				it2 = findex_inver.find( ZAI(-2,-2,-2) );
			else it2 = findex_inver.find( (*it).first );
			
			if(it2 == findex_inver.end() )		//If not in index should be TMP, can't be fast decay for new Fuel !!!
				it2 = findex_inver.find( ZAI(-3,-3,-3) );
			N_0Matrix[ (*it2).second ][0] = (*it).second ;
			
			
		}
		
		isotopicquantity.clear();
		
		NMatrix.push_back(N_0Matrix);
		N_0Matrix.Clear();
		
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
		isotopicquantity.clear();
		
	}
	
	int DBTimeStepN = EvolutionDataStep.GetFissionXS().begin()->second->GetN();
	double* DBTimeStep = EvolutionDataStep.GetFissionXS().begin()->second->GetX();
	
	int InsideStep = 10;
	
	int NStep = (DBTimeStepN);
	double timevector[NStep];
	timevector[0] = 0;
	
	double  Flux[NStep];
	
	TMatrixT<double> SigmaPhi = TMatrixT<double>(findex.size()*3+1,NStep); // Store the XS and the flux trought the evolution calculation.
	for(int i = 0; i < (int)findex.size()*3+1; i++)
		for(int j = 0; j < (int)NStep; j++)
			SigmaPhi[i][j] = 0;
	
	TMatrixT<double> FissionEnergy = TMatrixT<double>(findex.size(),1);
	for(int i = 0; i < (int)findex.size(); i++)
		FissionEnergy[i] = 0;
	
	{
		map< ZAI, int >::iterator it;
		for(it = findex_inver.begin(); it != findex_inver.end(); it++)
		{
			map< ZAI, double >::iterator it2 = fFissionEnergy.find(it->first);
			if(it2 == fFissionEnergy.end())
			{
				if(it->first.Z() > fZAIThreshold)
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
		
		BatemanMatrix.Clear();
		BatemanReactionMatrix.Clear();
		NEvolutionMatrix.Clear();
		
		
	}
	FissionXSMatrix.push_back(GetFissionXsMatrix(EvolutionDataStep, DBTimeStep[NStep-1])); //Feel the reaction Matrix
	CaptureXSMatrix.push_back(GetCaptureXsMatrix(EvolutionDataStep, DBTimeStep[NStep-1])); //Feel the reaction Matrix
	n2nXSMatrix.push_back(Getn2nXsMatrix(EvolutionDataStep, DBTimeStep[NStep-1])); //Feel the reaction Matrix
	
	
	EvolutionData GeneratedDB = EvolutionData(GetLog());
	
	double ESigmaN = 0;
	for (int j = 0; j < (int)findex.size() ; j++)
		ESigmaN -= FissionXSMatrix.back()[j][j]*NMatrix.back()[j][0]*1.6e-19*FissionEnergy[j][0];
	
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
	
	//	fFuelDataBankCalculated.insert( pair< IsotopicVector, EvolutionData > ( GeneratedDB.GetIsotopicVectorAt(0.), GeneratedDB) );
	
	ResetTheMatrix();
	ResetTheNucleiVector();
	
	for (int i = 0; i < (int) FissionXSMatrix.size(); i++)
	{
		FissionXSMatrix[i].Clear();
		CaptureXSMatrix[i].Clear();
		n2nXSMatrix[i].Clear();
	}
	FissionXSMatrix.clear();
	CaptureXSMatrix.clear();
	n2nXSMatrix.clear();
	
	if(fEvolutionDataInterpolation)
		EvolutionDataStep.DeleteEvolutionData();
	
	return GeneratedDB;
	
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
