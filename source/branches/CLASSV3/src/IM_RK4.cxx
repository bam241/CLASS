//
//  IM_RK4.cxx
//  CLASSSource
//
//  Created by BaM on 04/05/2014.
//  Copyright (c) 2014 BaM. All rights reserved.
//

#include "IM_RK4.hxx"

#include "IsotopicVector.hxx"
#include "CLASSConstante.hxx"
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



using namespace std;


//________________________________________________________________________
IM_RK4::IM_RK4():IrradiationModel(), DynamicalSystem()
{
	fTheNucleiVector = 0;
	fTheMatrix = 0;

	SetForbidNegativeValue();

}


IM_RK4::IM_RK4(LogFile* log):IrradiationModel(log), DynamicalSystem()
{
	fTheNucleiVector = 0;
	fTheMatrix = 0;

	SetForbidNegativeValue();

}









//________________________________________________________________________
/*			Evolution Calculation			*/
//________________________________________________________________________
EvolutionData IM_RK4::GenerateEvolutionData(IsotopicVector isotopicvector, EvolutionData XSSet, double Power, double cycletime)
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
	ReactorType = XSSet.GetReactorType();

	double M_ref = XSSet.GetHeavyMetalMass();
	double M = 0;
	double Power_ref =  XSSet.GetPower();

	// Get the mass of the fuel to irradiate in order to perform the evolution at fixed burnup (the burnup step of the calculation will match the burnup step of the XSSet
	{
		map<ZAI, double >::iterator it ;


		IsotopicVector IVtmp = isotopicvector.GetActinidesComposition();
		map<ZAI, double >isotopicquantity = IVtmp.GetIsotopicQuantity();

		for( it = isotopicquantity.begin(); it != isotopicquantity.end(); it++ )
			M += isotopicvector.GetActinidesComposition().GetZAIIsotopicQuantity( (*it).first )*cZAIMass.fZAIMass.find( (*it).first )->second/AVOGADRO*1e-6;
		isotopicquantity.clear();

	}

	int NStep = XSSet.GetFissionXS().begin()->second->GetN();
	double* DBTimeStep = XSSet.GetFissionXS().begin()->second->GetX();

	int InsideStep = 10;

	double timevector[NStep];
	timevector[0] = 0;

	double  Flux[NStep];

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

	vector< TMatrixT<double> > FissionXSMatrix;	// Store The Fisison XS Matrix
	vector< TMatrixT<double> > CaptureXSMatrix;	// Store The Capture XS Matrix
	vector< TMatrixT<double> > n2nXSMatrix;		// Store The n2N XS Matrix

	for(int i = 0; i < NStep-1; i++)
	{
		double TStepMax = ( (DBTimeStep[i+1]-DBTimeStep[i] ) ) * Power_ref/M_ref / Power*M ;	// Get the next Time step


		TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
		TMatrixT<double> BatemanReactionMatrix = TMatrixT<double>(findex.size(),findex.size());

		TMatrixT<double> NEvolutionMatrix = TMatrixT<double>(findex.size(),1);
		NEvolutionMatrix = NMatrix.back();



		FissionXSMatrix.push_back(GetFissionXsMatrix(XSSet, DBTimeStep[i]));	//Feel the fission reaction Matrix
		CaptureXSMatrix.push_back(GetCaptureXsMatrix(XSSet, DBTimeStep[i]));	//Feel the capture reaction Matrix
		n2nXSMatrix.push_back(Getn2nXsMatrix(XSSet, DBTimeStep[i]));		//Feel the (n,2n)  reaction Matrix

		// ----------------   Evolution

		BatemanReactionMatrix = FissionXSMatrix[i];
		BatemanReactionMatrix += CaptureXSMatrix[i];
		BatemanReactionMatrix += n2nXSMatrix[i];

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

		timevector[i+1] = timevector[i] + TStepMax;

		BatemanMatrix.Clear();
		BatemanReactionMatrix.Clear();
		NEvolutionMatrix.Clear();


	}
	FissionXSMatrix.push_back(GetFissionXsMatrix(XSSet, DBTimeStep[NStep-1])); //Feel the reaction Matrix
	CaptureXSMatrix.push_back(GetCaptureXsMatrix(XSSet, DBTimeStep[NStep-1])); //Feel the reaction Matrix
	n2nXSMatrix.push_back(Getn2nXsMatrix(XSSet, DBTimeStep[NStep-1])); //Feel the reaction Matrix


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
//	GeneratedDB.SetFuelType(fFuelType );
	GeneratedDB.SetReactorType(ReactorType );
	GeneratedDB.SetCycleTime(cycletime);

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

	return GeneratedDB;

}


void IM_RK4::ResetTheMatrix()
{

	if(fTheMatrix)
	{
		for(int i= 0; i<fNVar; i++)
			delete [] fTheMatrix[i];
		delete [] fTheMatrix;
	}
	fTheMatrix = 0;
}

void IM_RK4::SetTheMatrixToZero()
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
void IM_RK4::ResetTheNucleiVector()
{
	if(fTheNucleiVector)
		delete [] fTheNucleiVector;
	fTheNucleiVector = 0;
}

//________________________________________________________________________
void IM_RK4::SetTheNucleiVectorToZero()
{
	ResetTheNucleiVector();
	fTheNucleiVector = new double[fNVar];

#pragma omp parallel for
	for(int i = 0; i < fNVar; i++)
		fTheNucleiVector[i]=0.0;

}

//________________________________________________________________________
void IM_RK4::BuildEqns(double t, double *N, double *dNdt)
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
void IM_RK4::SetTheMatrix(TMatrixT<double> BatemanMatrix)
{
	for (int k = 0; k < (int)fNVar; k++)
		for (int l = 0; l < (int)findex_inver.size(); l++)
			fTheMatrix[l][k] = BatemanMatrix[l][k];
}

//________________________________________________________________________
TMatrixT<double> IM_RK4::GetTheMatrix()
{
	TMatrixT<double> BatemanMatrix = TMatrixT<double>(findex.size(),findex.size());
	for (int k = 0; k < (int)fNVar; k++)
		for (int l = 0; l < (int)findex_inver.size(); l++)
			BatemanMatrix[l][k] = fTheMatrix[l][k];

	return BatemanMatrix;
}

//________________________________________________________________________
void IM_RK4::SetTheNucleiVector(TMatrixT<double> NEvolutionMatrix)
{
	for (int k = 0; k < (int)fNVar; k++)
		fTheNucleiVector[k] = NEvolutionMatrix[k][0];
}

//________________________________________________________________________
TMatrixT<double> IM_RK4::GetTheNucleiVector()
{
	TMatrixT<double> NEvolutionMatrix = TMatrixT<double>(findex.size(),1);
	for (int k = 0; k < (int)fNVar; k++)
		NEvolutionMatrix[k][0] = fTheNucleiVector[k];
	
	return NEvolutionMatrix;
}
















