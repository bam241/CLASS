//
//  IrradiationModel.cxx
//  CLASSSource
//
//  Created by BaM on 04/05/2014.
//  Copyright (c) 2014 BaM. All rights reserved.
//

#include "IrradiationModel.hxx"

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
//________________________________________________________________________
/*				Physics				*/
//________________________________________________________________________
//________________________________________________________________________



//________________________________________________________________________
/*				Decay Stuff			*/
//________________________________________________________________________
string IrradiationModel::GetDecay(string DecayModes, double &BR,int &Iso, int &StartPos)
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
void IrradiationModel::BuildDecayMatrix()
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
		cout << "!!Warning!! !!!IrradiationModel!!! \n Can't open \"" << DataFullPathName << "\"\n" << endl;
		GetLog()->fLog << "!!Warning!! !!!IrradiationModel!!! \n Can't open \"" << DataFullPathName<< "\"\n" << endl;
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
void IrradiationModel::SetFissionEnergy(ZAI zai, double E)
{
	pair<map<ZAI, double>::iterator, bool> IResult;
	IResult = fFissionEnergy.insert( pair<ZAI ,double>(zai, E));
	if(!IResult.second)
		IResult.first->second = E;

}

//________________________________________________________________________
void IrradiationModel::SetFissionEnergy(string FissionEnergyFile)
{
	ifstream FissionFile(FissionEnergyFile.c_str());	// Open the File
	if(!FissionFile)				//check if file is correctly open
	{
		cout << "!!Warning!! !!!IrradiationModel!!! \n Can't open \"" << FissionFile << "\"\n" << endl;
		GetLog()->fLog << "!!Warning!! !!!IrradiationModel!!! \n Can't open \"" << FissionFile << "\"\n" << endl;
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
map< ZAI,IsotopicVector > IrradiationModel::ReadFPYield(string Yield)
{
	IsotopicVector EmptyIV;
	map< ZAI,IsotopicVector >  Yield_map;

	ifstream infile(Yield.c_str());
	if(!infile)
	{
		cout << "!!Warning!! !!!IrradiationModel!!! \n Can't open \"" << Yield << "\"\n" << endl;
		GetLog()->fLog << "!!Warning!! !!!IrradiationModel!!! \n Can't open \"" << Yield<< "\"\n" << endl;
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
			if(!IResult.second)
			{
				cout << "!!Error!! !!!IrradiationModel!!! Many accurance of ZAI " << Z << " " << A;
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
				cout << "!!Error!! !!!IrradiationModel!!! Many accurance of the PF " << Z << " " << A;
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
void IrradiationModel::LoadFPYield(string SponfaneusYield, string ReactionYield)
{

	fSpontaneusYield = ReadFPYield(SponfaneusYield);
	fReactionYield = ReadFPYield(ReactionYield);
	fZAIThreshold = 0;
}
