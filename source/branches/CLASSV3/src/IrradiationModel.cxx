//
//  IrradiationModel.cxx
//  CLASSSource
//
//  Created by BaM on 04/05/2014.
//  Copyright (c) 2014 BaM. All rights reserved.
//

#include "IrradiationModel.hxx"

#include "IsotopicVector.hxx"
#include "CLASSLogger.hxx"
#include "StringLine.hxx"

#include <TGraph.h>
#include <TString.h>


#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>


using namespace std;

IrradiationModel::IrradiationModel():CLASSObject()
{
	fShorstestHalflife = 3600.*24*2.;
	fZAIThreshold = 90;


	fDataDirectoryName = getenv("CLASS_PATH");
	fDataDirectoryName += "/data/";
	fDataFileName = "chart.JEF3T";
}

IrradiationModel::IrradiationModel(CLASSLogger* log):CLASSObject(log)
{
	fShorstestHalflife = 3600.*24*2.;
	fZAIThreshold = 90;


	fDataDirectoryName = getenv("CLASS_PATH");
	fDataDirectoryName += "/data/";
	fDataFileName = "chart.JEF3T";
}


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
void IrradiationModel::LoadDecay()
{
DBGL

	// Add TMP and PF as a stable nuclei in the normal decay lsit
	fNormalDecay.Add(ZAI(-3,-3,-3), IsotopicVector() );		// TMP
	fNormalDecay.Add(ZAI(-2,-2,-2), IsotopicVector() );		// PF



	string DataFullPathName = GetDataDirectoryName()+ GetDataFileName();
	ifstream infile(DataFullPathName.c_str());

	if(!infile)
		WARNING << " Can't open \"" << DataFullPathName<< "\"" << endl;


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
			IsotopicVector DaughtersIV;
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
				if (decay_name == "it")	{DM=12;stable=false;	daughter_N=(A-Z);	daughter_Z=Z;}
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
							DaughtersIV += BR * DaughterZAI;
							branch_test+=BR;
						}
						else if( DM <= 18)
						{
							if(fSpontaneusYield.size() == 0 )
							{
								DaughtersIV += 2*BR * ZAI(-2,-2,-2);

								branch_test_f += 2*BR;
							}
							else
							{

								map<ZAI, IsotopicVector>::iterator it_yield = fSpontaneusYield.find(ParentZAI);
								if(it_yield != fSpontaneusYield.end())
								{
									DaughtersIV += BR* (*it_yield).second ;
									branch_test_f += BR* (*it_yield).second.GetSumOfAll();

								}
								else
								{
									DaughtersIV += 2*BR * ZAI(-2,-2,-2);
									branch_test_f += 2*BR;
								}
							}
						}

					}

				}
				if (DM !=0)
					stable = false;
				// End of While loop
			}

			double btest = fabs(branch_test + branch_test_f/2.-1.0);
			if ( btest > 1e-8 && !stable )
				DaughtersIV = DaughtersIV / (branch_test+branch_test_f/2.);
			
				
					
			if (HalfLife < fShorstestHalflife && !stable)
				fFastDecay.Add( ParentZAI, DaughtersIV );	// Fill the FastDecay by the Daughter IV
			else if (stable)
			{
				fNormalDecay.Add(ParentZAI, IsotopicVector() ); // Fill the NormalDecay with a empty IV (mother is stable)
			}
			else
				fNormalDecay.Add( ParentZAI,  ln(2)/HalfLife * DaughtersIV ); // FIll the NormalDecay with the daughter IV scaled by the decay constante.


		}

	} while (!infile.eof());
	
	
	//Build the Matrix index :
	fReverseMatrixIndex = fNormalDecay.GetZAIList();
	for(int i = 0; i< (int)fReverseMatrixIndex.size(); i++)
		fMatrixIndex.insert(pair<ZAI, int> (fReverseMatrixIndex[i], i) );
	
	
	fFastDecay.SelfFiliationCleanUp(fMatrixIndex);
	fNormalDecay.FiliationCleanUp(fMatrixIndex, fFastDecay);


DBGL
}


void IrradiationModel::BuildDecayMatrix()
{
	
	fDecayMatrix.Clear();
	
	fDecayMatrix.ResizeTo(findex.size(),findex.size());
	for(int i = 0; i < (int)findex.size(); i++)
		for(int j = 0; j < (int)findex.size(); j++)
			fDecayMatrix[i][j] = 0;
	
	for(int i = 0; i < (int)fReverseMatrixIndex.size(); i++)
	{
		
		IsotopicVector DaughterIV = fNormalDecay.GetFiliation(fReverseMatrixIndex[i]);
		vector<ZAI> DaughterZAIList = DaughterIV.GetZAIList();
		
		for(int j = 0; j < (int)DaughterZAIList.size(); j++)
		{
			if(fMatrixIndex.find(DaughterZAIList[j]) != fMatrixIndex.end() )
				fDecayMatrix[fMatrixIndex[ DaughterZAIList[j] ]][i] += DaughterIV.GetQuantity(DaughterZAIList[j]);
			else
				fDecayMatrix[0][i] += DaughterIV.GetQuantity(DaughterZAIList[j]);
			
		}
		
	}
	
	DBGL
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
		WARNING << " Can't open \"" << FissionEnergyFile << "\"\n" << endl;

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
		WARNING << " Can't open \"" << Yield<< "\"\n" << endl;


	string line;
	int start = 0;

	getline(infile, line);
	vector<ZAI> Fissile;
	vector<IsotopicVector> FPYields;
	do
	{
		int Z = atof(StringLine::NextWord(line, start, ' ').c_str());
		int A = atof(StringLine::NextWord(line, start, ' ').c_str());
		int I = atof(StringLine::NextWord(line, start, ' ').c_str());
		Fissile.push_back(ZAI(Z,A,I));
		FPYields.push_back(EmptyIV);
		
	}while(start < (int)line.size()-1);

	getline(infile, line);
	do
	{
		start = 0;

		int Z = atof(StringLine::NextWord(line, start, ' ').c_str());
		int A = atof(StringLine::NextWord(line, start, ' ').c_str());
		int I = atof(StringLine::NextWord(line, start, ' ').c_str());
		int i=0;
		bool NucleiFound=false;
		do
		{
			double Yield_values = atof(StringLine::NextWord(line, start, ' ').c_str());
			//	cout<<"****Searching for ZAI : "<<Z<<" "<<A<<" "<<I<<endl;
			while(!NucleiFound && Z<200)
			{	//Is in the non cut nuclei chart ?
				map<ZAI, int>::iterator findNonCut_it = findex_inver.find( ZAI(Z,A,I) );
				if( findNonCut_it != findex_inver.end() )
					NucleiFound=true;
					//cout<<"Find in non cut"<<endl;
					
				else
				{	//Is in the cutted nuclei chart ?
					map<ZAI, map<ZAI,double> >::iterator findFastDecay_it = fFastDecay.find(ZAI(Z,A,I));
					if( findFastDecay_it != fFastDecay.end() )
						NucleiFound=true;
					
     				//if after a SF the FP does not exist, make either beta- (for gs) or Isomeric transition
     				// to find an existing one (in the chart)
					else
					{	//cout<<"Artificial decay ... "<<endl;
						if(I>0)//Isomeric Decay
							I--;
						else//Beta- decay
							Z++;	
						//cout<<"new ZAI "<<Z<<" "<<A<<" "<<I<<" "<<endl;				
					}
				}	
			}
			if(!NucleiFound)
			{	WARNING<<"Nuclei not found in chart "<<endl;}

			else
				FPYields[i] +=  Yield_values * ZAI(Z,A,I);

			i++;

		}while(start < (int)line.size()-1);

		getline(infile, line);

	} while (!infile.eof());


	for(int i=0 ; i<(int) Fissile.size() ; i++)
	{
		pair<map<ZAI, IsotopicVector>::iterator, bool> IResult;
		IResult = Yield_map.insert(pair<ZAI,IsotopicVector>(Fissile[i],FPYields[i]));
		if(!IResult.second)
		{
		ERROR << " Many occurances of ZAI  file!! " << Fissile[i].Z()<<" "<<Fissile[i].A()<<" "<<Fissile[i].I()<<" "<<endl;
		exit(1);
		}
	}

	return Yield_map;
}

//________________________________________________________________________
void IrradiationModel::LoadFPYield(string SpontaneusYield, string ReactionYield)
{
	fSpontaneusYieldFile = SpontaneusYield;
	fReactionYieldFile   = ReactionYield;	
	fZAIThreshold = 0;
}
//________________________________________________________________________
void IrradiationModel::NuclearDataInitialization()
{

	LoadDecay();
	
	BuildDecayMatrix();

	if(fSpontaneusYieldFile!="")
		fSpontaneusYield = ReadFPYield(fSpontaneusYieldFile);

	if(fReactionYieldFile!="")
		fReactionYield = ReadFPYield(fReactionYieldFile);
}
//________________________________________________________________________
/*				Reaction Stuff			*/
//________________________________________________________________________
TMatrixT<double> IrradiationModel::GetFissionXsMatrix(EvolutionData EvolutionDataStep,double TStep)
{
DBGL
	map<ZAI ,TGraph* >::iterator it_XS;
	TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
	for(int i = 0; i < (int)findex.size(); i++)
		for(int j = 0; j < (int)findex.size(); j++)
			BatemanMatrix[i][j] = 0;

	// ----------------  A(n,.) X+Y

	map<ZAI ,TGraph* > FissionXS = EvolutionDataStep.GetFissionXS();

	for(it_XS = FissionXS.begin() ; it_XS != FissionXS.end(); it_XS++)	//loop on fissionable nuclei
	{	
	//	cout<<"***************"<<(*it_XS).first.Z()<<" "<<(*it_XS).first.A()<<" "<<(*it_XS).first.I()<<endl;

		map<ZAI, int>::iterator findex_inver_it = findex_inver.find( (*it_XS).first );
		if( findex_inver_it != findex_inver.end() )
		{
			double XS_Value = (*it_XS).second->Eval(TStep);
			BatemanMatrix[ findex_inver_it->second ][ findex_inver_it->second ] += -XS_Value* 1e-24;

			if(fReactionYield.size() == 0)
				BatemanMatrix[1][ findex_inver_it->second ] += 2*XS_Value* 1e-24;
			else
			{
				map<ZAI, IsotopicVector>::iterator it_yield = fReactionYield.find( (*it_XS).first );

				if( it_yield != fReactionYield.end())
				{
					map<ZAI ,double>::iterator it_FissionProductMap;
					map<ZAI ,double> FissionProductMap = (*it_yield).second.GetIsotopicQuantity();

					for( it_FissionProductMap = FissionProductMap.begin(); it_FissionProductMap != FissionProductMap.end(); it_FissionProductMap++ )//loop on fission product
					{	//cout<<(*it_FissionProductMap).first.Z()<<" "<<(*it_FissionProductMap).first.A()<<" "<<(*it_FissionProductMap).first.I()<<" "<<(*it_FissionProductMap).second<<endl;						
						map<ZAI, int>::iterator findex_it_PF = findex_inver.find( (*it_FissionProductMap).first );

						if(findex_it_PF != findex_inver.end() )
						{	BatemanMatrix[(*findex_it_PF).second][ (*findex_inver_it).second ] += (*it_FissionProductMap).second*XS_Value* 1e-24;
							
						}	
						else
						{
							map<ZAI, map<ZAI, double> >::iterator it_FD = fFastDecay.find( (*it_FissionProductMap).first);

							if( it_FD == fFastDecay.end() )
							{
								BatemanMatrix[1][ (*findex_inver_it).second ] += (*it_FissionProductMap).second * XS_Value * 1e-24  ;
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
										(*it_FissionProductMap).second * XS_Value * 1e-24 * (*it5).second;
									else
									{	BatemanMatrix[(*findex_it_PF).second][findex_inver_it->second]+=
										(*it_FissionProductMap).second * XS_Value * 1e-24 * (*it5).second;
									}	
								}
							}

						}

					}
				}
				else
					BatemanMatrix[1][ findex_inver_it->second ] += 2*XS_Value* 1e-24;


			}
		}

	}

DBGL
	return BatemanMatrix;
}

//________________________________________________________________________
TMatrixT<double> IrradiationModel::GetCaptureXsMatrix(EvolutionData EvolutionDataStep,double TStep)
{
DBGL
	map<ZAI ,TGraph* >::iterator it;
	TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
	for(int i = 0; i < (int)findex.size(); i++)
		for(int j = 0; j < (int)findex.size(); j++)
			BatemanMatrix[i][j] = 0;

	map<ZAI, map<ZAI, double> > Capture;
	{	// 241Am
		map<ZAI, double> Am1Case ;
		Am1Case.insert(pair<ZAI, double> ( ZAI(96,242,0) , 0.8733*0.827) ); //directly cut the Am242 as in MURE
		Am1Case.insert(pair<ZAI, double> ( ZAI(94,242,0) , 0.8733*0.173) ); //directly cut the Am242 as in MURE
		Am1Case.insert(pair<ZAI, double> ( ZAI(95,242,1) , 0.1267) );
		Capture.insert( pair< ZAI, map<ZAI, double> > ( ZAI(95,241,0), Am1Case ) );
		// 165Ho
		map<ZAI, double> Ho165Case ;
		Ho165Case.insert(pair<ZAI, double> ( ZAI(68,166,0) , 0.9490) ); 
		Ho165Case.insert(pair<ZAI, double> ( ZAI(67,166,1) , 1-0.9490 ) ); 
		Capture.insert( pair< ZAI, map<ZAI, double> > ( ZAI(67,165,0), Ho165Case ) );
		// 147Pm
		map<ZAI, double> Pm147Case ;
		Pm147Case.insert(pair<ZAI, double> ( ZAI(61,148,0) , 0.5330) ); 
		Pm147Case.insert(pair<ZAI, double> ( ZAI(61,148,1) , 1-0.5330 ) ); 
		Capture.insert( pair< ZAI, map<ZAI, double> > ( ZAI(61,147,0), Pm147Case ) );
		// 109Ag
		map<ZAI, double> Ag109Case ;
		Ag109Case.insert(pair<ZAI, double> ( ZAI(48,110,0) , 0.9970*0.9508) ); 
		Ag109Case.insert(pair<ZAI, double> ( ZAI(46,110,0) , 0.0030*0.9508) ); 
		Ag109Case.insert(pair<ZAI, double> ( ZAI(47,110,1) , 1-0.9508 ) ); 
		Capture.insert( pair< ZAI, map<ZAI, double> > ( ZAI(47,109,0), Ag109Case ) );
		// 107Ag
		map<ZAI, double> Ag107Case ;
		Ag107Case.insert(pair<ZAI, double> ( ZAI(48,108,0) , 0.9715*0.9895) ); 
		Ag107Case.insert(pair<ZAI, double> ( ZAI(46,108,0) , 0.0285*0.9895) ); 
		Ag107Case.insert(pair<ZAI, double> ( ZAI(47,108,1) , 1-0.9895) ); 
		Capture.insert( pair< ZAI, map<ZAI, double> > ( ZAI(47,107,0), Ag107Case ) );

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
				map<ZAI, int >::iterator it6 = findex_inver.find( ZAI( (*it).first.Z(), (*it).first.A()+1, 0) );

				if( it6 != findex_inver.end() )
				{
					BatemanMatrix[(*it6).second][Index_it->second] += y* 1e-24  ;
				}
				else
				{
					map<ZAI, map<ZAI, double> >::iterator it4 = fFastDecay.find(  ZAI( (*it).first.Z(), (*it).first.A()+1, 0) );

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
							ERROR << " CaptureList Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
							exit(1);
						}

						map< ZAI, double >::iterator it5;
						map< ZAI, double > decaylist2 = (*it7).second;
						for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
						{

							it6 = findex_inver.find( (*it5).first );
							if( it6 == findex_inver.end() )
							{
								ERROR << " CaptureList Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
								exit(1);
							}

							BatemanMatrix[(*it6).second][Index_it->second] += y * 1e-24 * (*it5).second * (*it4).second;
						}
					}

				}
			}


		}
	}
DBGL
	return BatemanMatrix;
}


//________________________________________________________________________
TMatrixT<double> IrradiationModel::Getn2nXsMatrix(EvolutionData EvolutionDataStep,double TStep)
{
DBGL

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
				map<ZAI, int >::iterator it6 = findex_inver.find( ZAI( (*it).first.Z(), (*it).first.A()-1, 0) );

				if( it6 != findex_inver.end() )
				{
					BatemanMatrix[(*it6).second][Index_it->second] += y* 1e-24  ;
				}
				else
				{
					map<ZAI, map<ZAI, double> >::iterator it4 = fFastDecay.find(  ZAI( (*it).first.Z(), (*it).first.A()-1, 0) );

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
							ERROR << " n2nList Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
							exit(1);
						}

						map< ZAI, double >::iterator it5;
						map< ZAI, double > decaylist2 = (*it7).second;
						for(it5 = decaylist2.begin(); it5!= decaylist2.end(); it5++)
						{

							it6 = findex_inver.find( (*it5).first );
							if( it6 == findex_inver.end() )
							{
								ERROR << " n2nList Problem in FastDecay for nuclei " << (*it7).first.Z() << " " << (*it7).first.A() << " " << (*it7).first.I() << endl;
								exit(1);
							}

							BatemanMatrix[(*it6).second][Index_it->second] += y * 1e-24 * (*it5).second * (*it4).second;
						}
					}

				}
			}


		}
	}
DBGL
	return BatemanMatrix;
}

