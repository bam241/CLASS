#include "FabricationPlant.hxx"
#include "MyFabricationPlant.hxx"
#include "Storage.hxx"
#include "Reactor.hxx"
#include "EvolutionData.hxx"
#include "DataBank.hxx"
#include "IsotopicVector.hxx"
#include "CLASS.hxx"
#include "CLASSHeaders.hxx"
#include "LogFile.hxx"




#include "TMatrixT.h"

#include <sstream>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

	//________________________________________________________________________
	//________________________________________________________________________
	//________________________________________________________________________
	//
	//		MyFabricationPlant
	//	a FP specific to one type of Reactor and one type of Fuel
	//
	//
	//
	//________________________________________________________________________
	//________________________________________________________________________
template <class T>  T random(T a, T b) //peak random numebr between a and b
{
	double range = pow(2., 31);
	srand(time(NULL)); //initialize the srand
	return (T)a + (T)(b-a)*rand()/range;
}
MyFabricationPlant::MyFabricationPlant()
{
	fDecayDataBase = 0;
	fStorage = 0;
	fReUsable = 0;
}

MyFabricationPlant::MyFabricationPlant(LogFile* log)
{
	
	SetLog(log);
	fChronologicalTimePriority = false;
	SetCycleTime(-1);
	fUpdateReferenceDBatEachStep = false;
	fSubstitutionFuel = false;
	fDecayDataBase = 0;
	fStorage = 0;
	fReUsable = 0;

	cout	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	cout	<< "\t Chronological Stock Priority set! "<< endl << endl;
	cout	<< "!!WARNING!! !!!FabricationPlant!!! You need to set the different stock manually as well as the Fabrication Time Manualy !! " << endl;
	GetLog()->fLog	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	GetLog()->fLog	<< "\t Chronological Stock Priority set! "<< endl << endl;
	GetLog()->fLog	<< "!!WARNING!! !!!FabricationPlant!!! You need to set the different stock manually as well as the Fabrication Time Manualy !! " << endl;
	
	

}

MyFabricationPlant::MyFabricationPlant(LogFile* log, Storage* storage, Storage* reusable, double fabircationtime)
{
	
	SetLog(log);
	
	fChronologicalTimePriority = false;
	fUpdateReferenceDBatEachStep = false;
	fSubstitutionFuel = false;
	fDecayDataBase = 0;


	SetCycleTime((cSecond)fabircationtime );
	fStorage = storage;
	fReUsable = reusable;
	
	
	cout	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	cout	<< "\t Chronological Stock Priority has been set! "<< endl;
	cout	<< "\t Fabrication time set to \t " << (double)(GetCycleTime()/3600/24/365.25) << " year" << endl << endl;
	
	GetLog()->fLog	<< "!!INFO!! !!!FabricationPlant!!! A FabricationPlant has been define :" << endl;
	GetLog()->fLog	<< "\t Chronological Stock Priority has been set! "<< endl;
	GetLog()->fLog	<< "\t Fabrication time set to \t " << (double)(GetCycleTime()/3600/24/365.25) << " year" << endl << endl;
	
	


}



	//________________________________________________________________________
MyFabricationPlant::~MyFabricationPlant()
{
	
	
}





void MyFabricationPlant::BuildFuelForReactor(int ReactorId)
{
	//cout<<"INFO : This is a specific FabricationPlant"<<endl<<"GOOD JOB!!!"<<endl;
	DataBank<IsotopicVector>* FuelType = GetParc()->GetReactor()[ReactorId]->GetFuelType();
	string ReactorType ="PWR";	
	if(FuelType->GetFuelType() != "MyFuel" || ReactorType !="PWR")//Check if the reactor is the right type and use the right type of fuel
	{
		cout << "!!Bad Trouble!! !!!FabricationPlant!!! Try to do MOX with a not MOXed DB "<< endl;
		GetLog()->fLog << "!!Bad Trouble!! !!!FabricationPlant!!! Try to do MOX with a not MOXed DB" << endl;
		exit (1);
	}	
	double Na = 6.02214129e23;	//N Avogadro

	double HMmass = GetParc()->GetReactor()[ReactorId]->GetHeavyMetalMass();
	double BU = GetParc()->GetReactor()[ReactorId]->GetBurnUp();
	IsotopicVector FullUsedStock;
	IsotopicVector stock;
	
	bool FuelBuild = false;
	while(!FuelBuild)
	{
		FuelBuild = true;
		cout << "!!Bad Trouble!! !!!FabricationPlant!!! This FabricationPlant is not working. You have to complete the BuildFuel function before using it!! "<< endl;
		IsotopicVector IVBeginCycle;
		EvolutionData evolutiondb = BuildEvolutiveDB(ReactorId, IVBeginCycle);		
		
		AddCumulativeIVIn(IVBeginCycle);
		fInsideIV = IVBeginCycle;
		
	}

}
