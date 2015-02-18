//
//  IM_Matrix_Decay.cxx
//  CLASSSource
//
//  Created by BaM on 04/05/2014.
//  Copyright (c) 2014 BaM. All rights reserved.
//

#include "IM_Matrix_Decay.hxx"

#include "IsotopicVector.hxx"
#include "CLASSConstante.hxx"
#include "CLASSLogger.hxx"
#include "StringLine.hxx"

#include <TGraph.h>
#include <TString.h>
#include <TMatrixT.h>

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>



using namespace std;


//________________________________________________________________________
IM_Matrix_Decay::IM_Matrix_Decay(IsotopicVector IVList):IrradiationModel(new CLASSLogger("IM_Matrix_Decay.log"))
{
	
	fShorstestHalflife = 0;
	fZAIThreshold = 0;

	fIVList = IVList;

	cout << "Loading Decay !! Could take few minutes... Please wait !!" << endl;
	NuclearDataInitialization();
	
	fExponentialDecayMatrix.Clear();
	
	fExponentialDecayMatrix.ResizeTo( fMatrixIndex.size(), fMatrixIndex.size() );

	fExponentialDecayMatrix = ExponentialCalculation(fDecayMatrix);

}


IM_Matrix_Decay::IM_Matrix_Decay(CLASSLogger* log, IsotopicVector IVList):IrradiationModel(log)
{
	fShorstestHalflife = 0;
	fZAIThreshold = 0;
	fIVList = IVList;

	
	
	cout << "laoding Decay !! Could take fiew minutes...." << endl;
	NuclearDataInitialization();
	
	fExponentialDecayMatrix.Clear();
	
	fExponentialDecayMatrix.ResizeTo( fMatrixIndex.size(), fMatrixIndex.size() );
	
	fExponentialDecayMatrix = ExponentialCalculation(fDecayMatrix);


}


TMatrixT<double> IM_Matrix_Decay::ExponentialCalculation(TMatrixT<double> myMatrix)
{
	TMatrixT<double> IdMatrix = TMatrixT<double>(fReverseMatrixIndex.size(),fReverseMatrixIndex.size());
	for(int j = 0; j < (int)fReverseMatrixIndex.size(); j++)
		for(int k = 0; k < (int)fReverseMatrixIndex.size(); k++)
		{
			if(k == j)	IdMatrix[j][k] = 1;
			else 		IdMatrix[j][k] = 0;
		}
	

	TMatrixT<double> MatrixDL = TMatrixT<double>(fReverseMatrixIndex.size(),fReverseMatrixIndex.size());   // Order 0 Term from the DL : Id
	TMatrixT<double> MatrixDLTermN = TMatrixT<double>(fReverseMatrixIndex.size(),fReverseMatrixIndex.size());  // Addind it;
	
	{
		MatrixDLTermN = IdMatrix;
		MatrixDL = MatrixDLTermN;
		int j = 1;
		double NormN;

		do
		{

			TMatrixT<double> MatrixDLTermtmp = TMatrixT<double>(fReverseMatrixIndex.size(),fReverseMatrixIndex.size());  // Adding it;
			MatrixDLTermtmp = MatrixDLTermN;
			MatrixDLTermN.Mult(MatrixDLTermtmp, myMatrix );

			MatrixDLTermN *= 1./j;
			MatrixDL += MatrixDLTermN;
			
			NormN = 0;
			for(int m = 0; m < (int)fReverseMatrixIndex.size(); m++)
				for(int n = 0; n < (int)fReverseMatrixIndex.size(); n++)
					NormN += MatrixDLTermN[m][n]*MatrixDLTermN[m][n];
			j++;
			cout << j << " " << NormN << endl;
			
		} while ( NormN != 0 );
	}
	
	return MatrixDL;
}


IsotopicVector IM_Matrix_Decay::GetDecay(IsotopicVector Mother_IV, double time)
{
	
	IsotopicVector DecayIV;
	
	TMatrixT<double>  NMatrix_0 =  TMatrixT<double>(fReverseMatrixIndex.size(),1);
	TMatrixT<double>  NMatrix_t =  TMatrixT<double>(fReverseMatrixIndex.size(),1);
	for(int i = 0; i < (int)fReverseMatrixIndex.size(); i++)
	{
		NMatrix_0[i] = 0;
		NMatrix_t[i] = 0;
	}
	
	
	
	{	// Filling the t=0 State;
		map<ZAI, double > isotopicquantity = Mother_IV.GetIsotopicQuantity();
		
		map<ZAI, double >::iterator it ;
	
		
		for(it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
		{
			/// Need TO change with FP managment
			map<ZAI, int >::iterator it2;
			
			if( (*it).first.Z() < fZAIThreshold )
				it2 = fMatrixIndex.find( ZAI(-2,-2,-2) );
			else it2 = fMatrixIndex.find( (*it).first );
			
			if(it2 == fMatrixIndex.end() )		//If not in index should be TMP, can't be fast decay for new Fuel !!!
				it2 = fMatrixIndex.find( ZAI(-3,-3,-3) );
			
			
			NMatrix_0[ (*it2).second ][0] = (*it).second ;
			
			
		}
		
		isotopicquantity.clear();
	}
	
	
	TMatrixT<double> exp_T_Matrix = TMatrixT<double>(fReverseMatrixIndex.size(),fReverseMatrixIndex.size());
	for(int j = 0; j < (int)fReverseMatrixIndex.size(); j++)
		for(int k = 0; k < (int)fReverseMatrixIndex.size(); k++)
		{
			if(k == j)	exp_T_Matrix[j][k] = exp(time/8E-23);
			else 		exp_T_Matrix[j][k] = 0;
		}

	NMatrix_t = fExponentialDecayMatrix * exp_T_Matrix * NMatrix_0;
	
	for(int i = 0; i < (int)fReverseMatrixIndex.size(); i++)
	{
		DecayIV += fReverseMatrixIndex[i]* NMatrix_t[i][0];
		
	}
	return DecayIV;
	
}


//________________________________________________________________________
void IM_Matrix_Decay::LoadDecay()
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
				string decay_name = IrradiationModel::GetDecay(DecayModes, BR, Iso, start);
				
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
				fDecayConstante += ParentZAI*log(2)/HalfLife*8E-23;
				cout << log(2)/HalfLife*8E-23 << endl;
			}
		}
		
	} while (!infile.eof());
	DBGL

	CleanDecay();
	
	//Build the Matrix index :
	fReverseMatrixIndex = fNormalDecay.GetZAIList();
	for(int i = 0; i< (int)fReverseMatrixIndex.size(); i++)
		fMatrixIndex.insert(pair<ZAI, int> (fReverseMatrixIndex[i], i) );
	

	
	DBGL
	
}




//________________________________________________________________________
void IM_Matrix_Decay::NuclearDataInitialization()
{
	DBGL
	LoadDecay();
	BuildDecayMatrix();
	
	DBGL
}

void IM_Matrix_Decay::CleanDecay()
{
	DBGL
	map<ZAI, IsotopicVector> InitialNucleiFiliation = fNormalDecay.GetNucleiFIliation();		//! Map of all the pathway
	map<ZAI, IsotopicVector> CleanNucleiFiliation;
	map<ZAI, IsotopicVector> CleanNucleiFiliation_bis;
	
	vector<ZAI> ZAIList = fIVList.GetZAIList();
	
	map<ZAI, IsotopicVector>::iterator it;

	for (int i=0; i< (int)ZAIList.size(); i++)
	{
		it = InitialNucleiFiliation.find(ZAIList[i]);
		CleanNucleiFiliation.insert(*it);
	}
	
	bool insertEnded = false;
	
	while(!insertEnded)
	{
		CleanNucleiFiliation_bis = CleanNucleiFiliation;
		for(it = CleanNucleiFiliation_bis.begin(); it!=CleanNucleiFiliation_bis.end(); it++)
		{
			ZAIList = (*it).second.GetZAIList();
			map<ZAI, IsotopicVector>::iterator it2;
			for (int i=0; i< (int)ZAIList.size(); i++)
			{
				it2 = InitialNucleiFiliation.find(ZAIList[i]);
				pair<map<ZAI, IsotopicVector>::iterator, bool> IResult = CleanNucleiFiliation.insert(*it2);
				insertEnded = insertEnded || IResult.second;
			}
		}
		cout << "cleaned size "<< CleanNucleiFiliation.size() << endl;

	}
	
	fNormalDecay.SetNucleiFIliation(CleanNucleiFiliation);
	cout << "out" << endl;
	DBGL
}








