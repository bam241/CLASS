#include "XSModel.hxx"
#include "XSM_CLOSEST.hxx"
#include "LogFile.hxx"
#include "StringLine.hxx"


#include <TGraph.h>
#include <TString.h>
#include "TSystem.h"
#include "TROOT.h"

#include <sstream>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <cmath>




//________________________________________________________________________
//
//		XSM_CLOSEST
//
//
//
//
//________________________________________________________________________
XSM_CLOSEST::XSM_CLOSEST(LogFile* Log,string DB_index_file, bool oldreadmethod )
{
	SetLog(Log);
	fOldReadMethod = oldreadmethod;
	fDataBaseIndex = DB_index_file;
	fDistanceType = 0;
	fWeightedDistance = false;
	fEvolutionDataInterpolation = false;
	ReadDataBase();
	if(IsLog())
	{
		// Warning
		cout	<< "!!INFO!! !!!XSM_CLOSEST!!! A EvolutionData has been define :" << endl;
		cout	<< "\t His index is : \"" << DB_index_file << "\"" << endl;
		cout	<< "\t " << fFuelDataBank.size() << " EvolutionData have been read."<< endl << endl;
		
		GetLog()->fLog 	<< "!!INFO!! !!!XSM_CLOSEST!!! A EvolutionData has been define :" << endl;
		GetLog()->fLog	<< "\t His index is : \"" << DB_index_file << "\"" << endl;
		GetLog()->fLog	<< "\t " << fFuelDataBank.size() << " EvolutionData have been read."<< endl << endl;
	}

}
//________________________________________________________________________
XSM_CLOSEST::~XSM_CLOSEST()
{
	map<IsotopicVector ,EvolutionData >::iterator it_del;
	for( it_del = fFuelDataBank.begin(); it_del != fFuelDataBank.end(); it_del++)
		(*it_del).second.DeleteEvolutionData();
	fFuelDataBank.clear();


	for( it_del = fFuelDataBankCalculated.begin(); it_del != fFuelDataBankCalculated.end(); it_del++)
		(*it_del).second.DeleteEvolutionData();
	fFuelDataBankCalculated.clear();
}
//________________________________________________________________________
void XSM_CLOSEST::ReadDataBase()
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
/*			Distance Calculation			*/
//________________________________________________________________________
//________________________________________________________________________
//________________________________________________________________________
map<double, EvolutionData> XSM_CLOSEST::GetDistancesTo(IsotopicVector isotopicvector, double t) const
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
EvolutionData XSM_CLOSEST::GetCrossSections(IsotopicVector isotopicvector, double t) 
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
void XSM_CLOSEST::CalculateDistanceParameter()
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
void XSM_CLOSEST::SetDistanceParameter(IsotopicVector DistanceParameter)
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
void XSM_CLOSEST::SetDistanceType(int DistanceType)
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