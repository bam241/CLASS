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


#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>



using namespace std;


//________________________________________________________________________
IM_Matrix_Decay::IM_Matrix_Decay():IrradiationModel(new CLASSLogger("IM_Matrix_Decay.log"))
{
	fShorstestHalflife = 0;
	fZAIThreshold = 0;
	NuclearDataInitialization();
	fExponentialDecayMatrix = ExponentialCalculation(fDecayMatrix);

}


IM_Matrix_Decay::IM_Matrix_Decay(CLASSLogger* log):IrradiationModel(log)
{
	fShorstestHalflife = 0;
	fZAIThreshold = 0;
	NuclearDataInitialization();
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
			
		} while ( NormN != 0 );
	}
	
	return MatrixDL;
}


IsotopicVector IM_Matrix_Decay::GetDecay(IsotopicVector Mother_IV, time)
{
	
	
	TMatrixT<double>  NMatrix =  TMatrixT<double>(decayindex.size(),1))
	for(int i = 0; i < (int)fReverseMatrixIndex.size(); i++)
		NMatrix[i] = 0;
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
			
			
			NMatrix[ (*it2).second ][0] = (*it).second ;
			
			
		}
		
		isotopicquantity.clear();
	}
	
	
}
















