//
//  IrradiationModel.cxx
//  CLASSSource
//
//  Created by BaM on 04/05/2014.
//  Copyright (c) 2014 BaM. All rights reserved.
//

#include "IrradiationModel.hxx"

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

	fNormalDecay = CLASSNucleiFiliation( log );
	fFastDecay   = CLASSNucleiFiliation( log );

	fCaptureReaction= CLASSNucleiFiliation( log );
	fn2nReaction	= CLASSNucleiFiliation( log );
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
	
	BR = 0;
	string DecayBR = StringLine::NextWord(DecayModes,StartPos,',');	//extraction of the decay mode and the BR
	
	int ss = 0;
	string Decay = StringLine::NextWord(DecayBR,ss,':'); 	//extraction of the decay
	
	
	if( ss < (int)DecayBR.size() ) 	//extraction of the BR if exist (i.e. for non stable isotop)
		BR = atof(DecayBR.substr(ss+1).c_str());
	
	BR /= 100.;	//BR in % -> BR
	
	
	Iso = 0;
	if( Decay.find("/",0) < string::npos )	//find the Isomeric state of Daughter
	{
		Iso = atoi( Decay.substr( Decay.find("/") + 1 ).c_str() );
		Decay = Decay.substr( 0, Decay.find("/") );
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
			{
				if(Iname[0]=='m')
				{
					if( atoi( Iname.substr(1).c_str() )==0 )
						I = 1;
					else
						I = atoi( Iname.substr(1).c_str() );
				}
			}
			
			int start = 0;
			double branch_test = 0;
			double branch_test_f = 0;
			
			ZAI ParentZAI = ZAI(Z,A,I);
			IsotopicVector DaughtersIV;
			bool stable = true;
			
			while(start<int(DecayModes.size()))
			{
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
					
					if((BR>1e-10) && (!stable))
					{
						if(DM <= 15)
						{
							ZAI DaughterZAI = ZAI(daughter_Z,daughter_A,Iso);
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
								
								IsotopicVector SpontanuesFissionProduct = fSpontaneusYield.GetFiliation( ParentZAI );
								if(SpontanuesFissionProduct.GetQuantity(-1,-1,-1) == 0)
								{
									DaughtersIV += BR* SpontanuesFissionProduct ;
									branch_test_f += BR* SpontanuesFissionProduct.GetSumOfAll();
									
								}
								else
								{
									WARNING << " Unknwon Spontanues yield for ZAI : " << ParentZAI.Z() << " " << ParentZAI.A() << " " << ParentZAI.I() << endl;
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
			{
				fNormalDecay.Add( ParentZAI,  DaughtersIV ); // FIll the NormalDecay with the daughter IV scaled by the decay constante.
				fDecayConstante += ParentZAI*log(2)/HalfLife;
			}
		}
		
	} while (!infile.eof());
	DBGL
	
	//Build the Matrix index :
	fReverseMatrixIndex = fNormalDecay.GetZAIList();
	for(int i = 0; i< (int)fReverseMatrixIndex.size(); i++)
		fMatrixIndex.insert(pair<ZAI, int> (fReverseMatrixIndex[i], i) );
	
	DBGL
	fFastDecay.SelfFiliationCleanUp(fMatrixIndex);
	DBGL
	fNormalDecay.FiliationCleanUp(fMatrixIndex, fFastDecay);
	
	
	DBGL
}


void IrradiationModel::BuildDecayMatrix()
{
	DBGL
	
	fDecayMatrix.Clear();
	
	fDecayMatrix.ResizeTo( fMatrixIndex.size(), fMatrixIndex.size() );
	for(int i = 0; i < (int)fMatrixIndex.size(); i++)
		for(int j = 0; j < (int)fMatrixIndex.size(); j++)
			fDecayMatrix[i][j] = 0;
	
	for(int i = 0; i < (int)fReverseMatrixIndex.size(); i++)
	{
		
		IsotopicVector DaughterIV = fNormalDecay.GetFiliation(fReverseMatrixIndex[i]);
		(*this).GetNuclearProcessMatrix(fDecayMatrix, fReverseMatrixIndex[i], DaughterIV, fDecayConstante.GetQuantity(fReverseMatrixIndex[i]) );
		
	}
	
	DBGL
}


void IrradiationModel::GetNuclearProcessMatrix(TMatrixT<double> &NuclearProcessMatrix, ZAI Mother, IsotopicVector ProductedIV, double XSValue)
{
	DBGL
		
	vector<ZAI> ProductedZAIList = ProductedIV.GetZAIList();
	
	if(fMatrixIndex.find(Mother) != fMatrixIndex.end())
	{
		
		
		int i = fMatrixIndex[Mother];
		
		NuclearProcessMatrix[i][i] -= XSValue;
		
		for(int j = 0; j < (int)ProductedZAIList.size(); j++)
		{
			if(fMatrixIndex.find(ProductedZAIList[j]) != fMatrixIndex.end() )
			{	NuclearProcessMatrix[fMatrixIndex[ ProductedZAIList[j] ]][i] += ProductedIV.GetQuantity(ProductedZAIList[j])*XSValue;
				//cout<<ProductedIV.GetQuantity(ProductedZAIList[j])*XSValue<<endl;
			}	
			else
				NuclearProcessMatrix[0][i] += ProductedIV.GetQuantity(ProductedZAIList[j])*XSValue;
		}
	}
	else
		WARNING << " Can't have nuclear process on this nucleus, ZAI : " << Mother.Z() << " " << Mother.A() << " " << Mother.I() << " its halflife seems to be below the threshold!" << endl;
		
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
	ifstream FissionFile( FissionEnergyFile.c_str() );	// Open the File
	if(!FissionFile)				//check if file is correctly open
		WARNING << " Can't open \"" << FissionEnergyFile << "\" File" << endl;
	
	do
	{
		int Z = 0;
		int A = 0;
		int I = 0;
		double E = 0;
		FissionFile >> Z >> A >> I >> E;
		SetFissionEnergy(Z, A, I, E);
	} while (!FissionFile.eof());
}

//________________________________________________________________________
CLASSNucleiFiliation IrradiationModel::ReadFPYield(string Yield)
{
	DBGL
	
	CLASSNucleiFiliation MyYield = CLASSNucleiFiliation( fLog );
	
	ifstream infile(Yield.c_str());
	if(!infile)
		WARNING << " Can't open \"" << Yield<< "\" File !!!" << endl;
	
	
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
		FPYields.push_back( IsotopicVector() );
		
	}while(start < (int)line.size()-1);
	
	getline(infile, line);
	do
	{
		start = 0;
		
		int Z = atof(StringLine::NextWord(line, start, ' ').c_str());
		int A = atof(StringLine::NextWord(line, start, ' ').c_str());
		int I = atof(StringLine::NextWord(line, start, ' ').c_str());
		int i=0;

		do
		{
			double Yield_values = atof(StringLine::NextWord(line, start, ' ').c_str());
			
			FPYields[i] +=  Yield_values * ZAI(Z,A,I);
			
			i++;
			
		}while(start < (int)line.size()-1);
		
		getline(infile, line);
		
	} while (!infile.eof());
	
	
	for(int i=0 ; i<(int) Fissile.size() ; i++)			// Fill the CLASSNucleiFiliation
		MyYield.Add( Fissile[i], FPYields[i] );
	
	
	DBGL
	return MyYield;
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
	DBGL
	
	if(fSpontaneusYieldFile!="")
		fSpontaneusYield = ReadFPYield(fSpontaneusYieldFile);
	if(fReactionYieldFile!="")
		fReactionYield = ReadFPYield(fReactionYieldFile);

	LoadDecay();
	
	
	if(fSpontaneusYieldFile!="")
		fSpontaneusYield.FiliationCleanUp(fMatrixIndex, fFastDecay);		// remove the cutted nuclei....
	if(fReactionYieldFile!="")
		fReactionYield.FiliationCleanUp(fMatrixIndex, fFastDecay);		// remove the cutted nuclei....

	
	BuildDecayMatrix();
	BuildReactionFiliation();
	
	DBGL
}


void 	IrradiationModel::BuildReactionFiliation()
{
	DBGL
	
	// (n,Gamma) Special Reaction.....
	{
		// 241Am(n,Gamma)
		{
			fCaptureReaction.Add( ZAI(95,241,0), ZAI(96,242,0) * 0.8733*0.827 ); //directly cut the Am242 as in MURE
			fCaptureReaction.Add( ZAI(95,241,0), ZAI(94,242,0) * 0.8733*0.173 ); //directly cut the Am242 as in MURE
			fCaptureReaction.Add( ZAI(95,241,0), ZAI(95,242,1) * 0.1267 );
		}
		
		// 165Ho(n,Gamma)
		{
			fCaptureReaction.Add( ZAI(67,165,0), ZAI(68,166,0) * 0.9490 ); //
			fCaptureReaction.Add( ZAI(67,165,0), ZAI(67,166,1) * 0.0510 ); //
		}
		
		// 147Pm(n,Gamma)
		{
			fCaptureReaction.Add( ZAI(61,147,0), ZAI(61,148,0) * 0.5330 );
			fCaptureReaction.Add( ZAI(61,147,0), ZAI(61,148,1) * 0.4670 );
			
		}
		
		// 109Ag(n, Gamma)
		{
			fCaptureReaction.Add( ZAI(47,109,0), ZAI(48,110,0) * 0.9970*0.9508);
			fCaptureReaction.Add( ZAI(47,109,0), ZAI(46,110,0) * 0.0030*0.9508);
			fCaptureReaction.Add( ZAI(47,109,0), ZAI(47,110,1)         *0.0492);
		}
		
		
		// 107Ag(n, Gamma)
		{
			fCaptureReaction.Add( ZAI(47,107,0), ZAI(48,108,0) * 0.9715*0.9895 );
			fCaptureReaction.Add( ZAI(47,107,0), ZAI(46,108,0) * 0.0285*0.9895 );
			fCaptureReaction.Add( ZAI(47,107,0), ZAI(47,108,1)         *0.0105 );
		}
	}
	
	
	
	// (n,2n) Special Reaction.....
	{
		// 237Np(n,2n)
		{
			fn2nReaction.Add(ZAI(93,237,0), ZAI(93,236,0) * 0.2 );
			fn2nReaction.Add(ZAI(93,237,0), ZAI(93,236,1) * 0.8 );
		}
		
	}
	
	
	
	for(int i = 2; i < (int)fReverseMatrixIndex.size(); i++)	// Start at 2 to skeep "TMP" ZAI and "PF" ZAI
	{
		
		int Z = fReverseMatrixIndex[i].Z();
		int A = fReverseMatrixIndex[i].A();

		if(fCaptureReaction.GetFiliation(fReverseMatrixIndex[i]).GetQuantity(ZAI(-1,-1,-1))  == 1 )
		{
			fCaptureReaction.Add(fReverseMatrixIndex[i], ZAI(Z,A+1)*1);
		}
		if(fn2nReaction.GetFiliation(fReverseMatrixIndex[i]).GetQuantity(ZAI(-1,-1,-1))  == 1 )
		{
			if(A>1)
				fn2nReaction.Add(fReverseMatrixIndex[i], ZAI(Z,A-1)*1);
		}
	}
	
	
	fCaptureReaction.FiliationCleanUp(fMatrixIndex, fFastDecay);	// clean the filiation link
	fCaptureReaction.NormalizeBranchingRatio();			// normalize it
	
	
	fn2nReaction.FiliationCleanUp(fMatrixIndex, fFastDecay);	// clean the filiation link
	fn2nReaction.NormalizeBranchingRatio();				// normalize it
	
	
	DBGL
}


//________________________________________________________________________
/*				Reaction Stuff			*/
//________________________________________________________________________
TMatrixT<double> IrradiationModel::GetFissionXsMatrix(EvolutionData EvolutionDataStep,double TStep)
{
	DBGL
	TMatrixT<double> FissionMatrix = TMatrixT<double>( fReverseMatrixIndex.size(),fReverseMatrixIndex.size() );
	for(int i = 0; i < (int)fReverseMatrixIndex.size(); i++)
		for(int j = 0; j < (int)fReverseMatrixIndex.size(); j++)
			FissionMatrix[i][j] = 0;
	
	// ----------------  A(n,.) X+Y
	
	map<ZAI ,TGraph* > FissionXS = EvolutionDataStep.GetFissionXS();
	map<ZAI ,TGraph* >::iterator it_XS;
	
	for(it_XS = FissionXS.begin() ; it_XS != FissionXS.end(); it_XS++)	//loop on fissionable nuclei
	{
		ZAI Mother = (*it_XS).first;					// Note the Mother ZAI (not necessary but help for reading the code)
		double XS_Value = (*it_XS).second->Eval(TStep) * 1e-24;		// Get Cross section values
		
		IsotopicVector FissionProductIV = fReactionYield.GetFiliation(Mother);		// Get the Isotopicvector produced by the reaction
		
		if(FissionProductIV.GetQuantity(ZAI(-1,-1,-1)) != 1)						// Check if ZAI is dealed
			GetNuclearProcessMatrix( FissionMatrix, Mother, FissionProductIV,  XS_Value );	// add the Nuclear process in the Reaction Matrix
		else
		{
			WARNING << "Don't have fission Yield for this nuclei, ZAI : " << Mother.Z() << " " << Mother.A() << " " << Mother.I() << endl;
			GetNuclearProcessMatrix(FissionMatrix, Mother, ZAI(-2, -2, -2) * 2 , XS_Value );	// add the Nuclear process in the Reaction Matrix
		}
		
		
	}
	
	DBGL
	return FissionMatrix;
}

//________________________________________________________________________
TMatrixT<double> IrradiationModel::GetCaptureXsMatrix(EvolutionData EvolutionDataStep,double TStep)
{
	DBGL
	TMatrixT<double> CaptureMatrix = TMatrixT<double>( fReverseMatrixIndex.size(),fReverseMatrixIndex.size() );
	for(int i = 0; i < (int)fReverseMatrixIndex.size(); i++)
		for(int j = 0; j < (int)fReverseMatrixIndex.size(); j++)
			CaptureMatrix[i][j] = 0;
	
	// ----------------  A(n,Gamma) A+1
	
	map<ZAI ,TGraph* > CaptureXS = EvolutionDataStep.GetCaptureXS();
	map<ZAI ,TGraph* >::iterator it_XS;
	
	for(it_XS = CaptureXS.begin() ; it_XS != CaptureXS.end(); it_XS++)	//loop on nuclei
	{
		ZAI Mother = (*it_XS).first;					// Note the Mother ZAI (not necessary but help for reading the code)
		double XS_Value = (*it_XS).second->Eval(TStep) * 1e-24;		// Get Cross section values
		
		IsotopicVector CaptureProductIV = fCaptureReaction.GetFiliation(Mother);		// Get the Isotopicvector produced by the reaction
		
		if(CaptureProductIV.GetQuantity(ZAI(-1,-1,-1)) != 1)						// Check if ZAI is dealed
			GetNuclearProcessMatrix(CaptureMatrix, Mother, CaptureProductIV,  XS_Value );	// add the Nuclear process in the Reaction Matrix
		else
			WARNING << "Can't have capture reaction on this nuclei, ZAI : " << Mother.Z() << " " << Mother.A() << " " << Mother.I() << endl;
			
		
	
		
		
	}
	
	DBGL
	return CaptureMatrix;
}


//________________________________________________________________________
TMatrixT<double> IrradiationModel::Getn2nXsMatrix(EvolutionData EvolutionDataStep,double TStep)
{
	DBGL
	
	TMatrixT<double> n2nMatrix = TMatrixT<double>( fReverseMatrixIndex.size(),fReverseMatrixIndex.size() );
	for(int i = 0; i < (int)fReverseMatrixIndex.size(); i++)
		for(int j = 0; j < (int)fReverseMatrixIndex.size(); j++)
			n2nMatrix[i][j] = 0;
	
	// ----------------  A(n,2n) A-1
	
	map<ZAI ,TGraph* > CaptureXS = EvolutionDataStep.Getn2nXS();
	map<ZAI ,TGraph* >::iterator it_XS;
	
	for(it_XS = CaptureXS.begin() ; it_XS != CaptureXS.end(); it_XS++)	//loop on nuclei
	{
		ZAI Mother = (*it_XS).first;					// Note the Mother ZAI (not necessary but help for reading the code)
		double XS_Value = (*it_XS).second->Eval(TStep) * 1e-24;		// Get Cross section values
		
		IsotopicVector n2nProductIV = fn2nReaction.GetFiliation(Mother);		// Get the Isotopicvector produced by the reaction
		
		if(n2nProductIV.GetQuantity(ZAI(-1,-1,-1)) != 1)						// Check if ZAI is dealed
			GetNuclearProcessMatrix(n2nMatrix, Mother, n2nProductIV,  XS_Value );	// add the Nuclear process in the Reaction Matrix
		else
			WARNING << "Can't have n,2n reaction on this nuclei, ZAI : " << Mother.Z() << " " << Mother.A() << " " << Mother.I() << endl;
		
		
	}
	
	DBGL
	return n2nMatrix;
}

