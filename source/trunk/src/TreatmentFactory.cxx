#include "TreatmentFactory.hxx"
#include "EvolutionDataBase.hxx"
#include "Defines.hxx"


#include <sstream>
#include <string>
#include <iostream>
#include <cmath>
#include <algorithm>

//________________________________________________________________________
//
//		TreatmentFactory
//
//
//
//
//________________________________________________________________________
ClassImp(TreatmentFactory)

TreatmentFactory::TreatmentFactory()
{
	DBGL;
	DBGL;
}

//________________________________________________________________________
TreatmentFactory::TreatmentFactory(long int abstime,
				   long int coolingtime ,
				   long int separationtime, 
				   EvolutionDataBase* evolutivedb)
{
	DBGL;
	fStockManagement = true;
	
	fCoolingTime = coolingtime;
	fInternalTime = 0;
	fDecayDataBase = evolutivedb;
	fSeparationTime = separationtime;
	fCreationTime = abstime;
	IsStarted = false;
	fCoolingLastIndex = 0;
	fSeparatingLastIndex = 0;
	DBGL;
}

//________________________________________________________________________
TreatmentFactory::~TreatmentFactory()
{
	DBGL;
	DBGL;
}



//________________________________________________________________________
void	TreatmentFactory::AddValorisableIV(ZAI zai, double factor)
{
	DBGL;
	pair<map<ZAI, double>::iterator, bool> IResult;
	if(factor > 1) factor = 1;
	
	if(factor > 0)
	{	
		IResult = fValorisableIV.insert( pair<ZAI ,double>(zai, factor));
		if(IResult.second == false)
			IResult.first->second = factor;
	}
	DBGL;
}


//________________________________________________________________________
//	Get Decay
//________________________________________________________________________
IsotopicVector TreatmentFactory::GetDecayProduct(IsotopicVector isotopicvector, long int t)
{
	DBGL;
	IsotopicVector IV;
	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	for(map<ZAI ,double >::iterator it = isotopicquantity.begin(); 
			it != isotopicquantity.end(); it++)
	{
		if((*it).second > 0)
		{
 			IsotopicVector ivtmp = fDecayDataBase->DecayProduction(it->first, t) * (*it).second ;
			IV += ivtmp;
		}
	}
	DBGL;
	return IV;
	DBGL;
}


//________________________________________________________________________
//	Add Temporary IV : 
//		Cooling
//		Separation
//________________________________________________________________________
void TreatmentFactory::AddIVCooling(IsotopicVector IV)
{ 
	DBGL;
	fIVCooling.push_back(IV);
	fCoolingStartingTime.push_back(fInternalTime);
	fCoolingLastIndex++;
	fCoolingIndex.push_back(fCoolingLastIndex);
	fParc->AddTotalCooling(IV);
	DBGL;
}

//________________________________________________________________________
void TreatmentFactory::AddIVSeparating(IsotopicVector IV)
{ 
	DBGL;
	fIVSeparating.push_back(IV); 
	fSeparatingStartingTime.push_back(fInternalTime);
	fSeparatingLastIndex++;
	fSeparatingIndex.push_back(fSeparatingLastIndex);
	fParc->AddTotalSeparating(IV);
	DBGL;
}

//________________________________________________________________________
void TreatmentFactory::AddIVSeparating(IsotopicVector IV, long int absolutadditiontime)
{ 
	DBGL;
	fIVSeparating.push_back(IV); 
	fSeparatingStartingTime.push_back(absolutadditiontime);
	fSeparatingLastIndex++;
	fSeparatingIndex.push_back(fSeparatingLastIndex);
	DBGL;
}


//________________________________________________________________________
void TreatmentFactory::RemoveIVSeparation(int i)	//remove a Treated IsotopicVector
{ 
	DBGL;
	fIVSeparating.erase(fIVSeparating.begin()+i);
	fSeparatingStartingTime.erase(fSeparatingStartingTime.begin()+i);
	fSeparatingIndex.erase(fSeparatingIndex.begin()+i);
	DBGL;
}


//________________________________________________________________________
void TreatmentFactory::RemoveIVCooling(int i)		//!< Remove a Cooling IsotopicVector
{
	DBGL;
	fIVCooling.erase(fIVCooling.begin()+i);
	fCoolingStartingTime.erase(fCoolingStartingTime.begin()+i);
	fCoolingIndex.erase(fCoolingIndex.begin()+i); 
	DBGL;
}



//________________________________________________________________________
//	Time Action with the reste of the Universe : 
//		In/Out stock
//		Evolution : 
//			Waste, Stock, Separating, Cooling
//		Separation
//________________________________________________________________________
void TreatmentFactory::TakeFromStock(IsotopicVector isotopicvector, int index)
{
	DBGL;
	if(fStockManagement == true )
		fIVStock.at(index) -= isotopicvector;
	
	fIVFullStock.Remove(isotopicvector);
	DBGL;
}

//________________________________________________________________________
pair<IsotopicVector, IsotopicVector> TreatmentFactory::Separation(IsotopicVector isotopicvector)
{
	DBGL;
	//[0] = stock ; [1] = waste
	pair<IsotopicVector, IsotopicVector>	IVTmp;
	
	map<ZAI ,double> isotopicquantity = isotopicvector.GetIsotopicQuantity();
	for(map<ZAI ,double >::iterator it = isotopicquantity.begin(); it != isotopicquantity.end(); it++)
	{
		map<ZAI ,double>::iterator it2;
		it2 = fValorisableIV.find((*it).first);
		if ( it2 != fValorisableIV.end() )
		{
			IVTmp.first.Add(	(*it).first, (*it).second * (*it2).second );		//stock
			IVTmp.second.Add(	(*it).first, (*it).second * (1-(*it2).second) );	//waste
		}
		else IVTmp.second.Add(	(*it).first, (*it).second );	//waste
	}
	DBGL;
	return IVTmp;
}


//________________________________________________________________________
void TreatmentFactory::WasteDecay(long int t)
{
	DBGL;
	// Check if the TF has been created ...
	if(t<fCreationTime) return;

	long int EvolutionTime = t - fInternalTime;
	fIVWaste = GetDecayProduct(fIVWaste , EvolutionTime);
	#pragma omp critical(UpdateTotalWasta)
		{fParc->AddTotalWaste(fIVWaste);}
	DBGL;
}

//________________________________________________________________________
void TreatmentFactory::StockDecay(long int t)
{
	DBGL;
	// Check if the TF has been created ...
	if(t<fCreationTime) return;

	int EvolutionTime = t - fInternalTime;
	fIVFullStock = GetDecayProduct(fIVFullStock , EvolutionTime);

	if(fStockManagement == true)
	{
		#pragma omp parallel for
		for (int i=0; i <(int) fIVStock.size() ; i++)
		{
			fIVStock.at(i) = GetDecayProduct(fIVStock.at(i) , EvolutionTime);
			#pragma omp critical(UpdateTotalStock)
				{fParc->AddTotalStock(fIVStock.at(i));}
		}
	
	}
	DBGL;
}


//________________________________________________________________________
void TreatmentFactory::SeparatingEvolution(long int t)
{
	DBGL;
	// Check if the TF has been created ...
	if(t<fCreationTime) return;

	
	long int EvolutionTime = t - fInternalTime;
#pragma omp parallel for
	for (int i = 0 ; i < (int)fIVSeparating.size() ; i++)
	{
		if (fInternalTime - fSeparatingStartingTime.at(i) + EvolutionTime >= fSeparationTime) // ">" should not append, only "=" is normal...
		{
			if (t - fSeparatingStartingTime.at(i) > fSeparationTime) // Warning & Quit
			{
				cout << "!!Warning!! !!!TreamtmentFactory!!! Separation Step : "<< t/365.4/3600/24 << " :"
				     << " An evolution Step is probably missing ! " << endl;
				exit (1);
			}
			fIVSeparating.at(i) = GetDecayProduct( fIVSeparating.at(i) , EvolutionTime);
		#pragma omp critical(DeleteSeprationIVPB)
			{fSeparationEndOfCycle.push_back(i);}
		}
		else 
		{
			fIVSeparating.at(i) = GetDecayProduct( fIVSeparating.at(i) , EvolutionTime);
			#pragma omp critical(UpdateSeparatingStock)
				{fParc->AddTotalSeparating(fIVSeparating.at(i));}
		}
	}
	sort (fSeparationEndOfCycle.begin(), fSeparationEndOfCycle.end());
	DBGL;
}

//________________________________________________________________________
void TreatmentFactory::CoolingEvolution(long int t)
{
	DBGL;
	// Check if the TF has been created ...
	if(t<fCreationTime) return;

	int i;
	int RemainingCoolingTime;
	long int EvolutionTime = t - fInternalTime;
#pragma omp parallel for
	for ( i = 0 ; i < (int)fIVCooling.size() ; i++)
	{
		if (t - fCoolingStartingTime.at(i) >= fCoolingTime) // ">" should not append, only "=" is normal...
		{
			if (t - fCoolingStartingTime.at(i) > fCoolingTime) // Warning & Quit
			{
				cout << "!!Warning!! !!!TreamtmentFactory!!! Cooling Step : " << t/365.4/3600/24 << " :"
				     << " An evolution Step is probably missing ! " << " " << endl;
				exit (1);
			}

			RemainingCoolingTime = fCoolingTime - (fInternalTime - fCoolingStartingTime.at(i));
			//Cooling Decay
			fIVCooling.at(i) = GetDecayProduct( fIVCooling.at(i), RemainingCoolingTime);
		#pragma omp critical(DeleteCoolingIVPB)
			{fCoolingEndOfCycle.push_back(i);}

		}
		else 
		{
			fIVCooling.at(i) = GetDecayProduct( fIVCooling.at(i) , EvolutionTime);
			#pragma omp critical(UpdateCoolingStock)
				{fParc->AddTotalCooling(fIVCooling.at(i));}
		}
	}

	sort (fCoolingEndOfCycle.begin(), fCoolingEndOfCycle.end());
	DBGL;
}


//________________________________________________________________________
void TreatmentFactory::Evolution(long int t)
{
	DBGL;
	// Check if the TF has been created ...
	if(t<fCreationTime) return;
	if(t == fInternalTime) return;
	
	if(fInternalTime ==0 && IsStarted == false)
	{
		fInternalTime = fCreationTime;
		IsStarted = true;
	}
	// Make the evolution for the Waste ...
	WasteDecay(t);
	
	// ... the Stock ...
	StockDecay(t);
	

	// ... the SeparatingIV ...
	SeparatingEvolution(t);
	
	// ... Then Deal the Cooling IV ...
	CoolingEvolution(t);
	
	// ... And Finaly update the AbsoluteInternalTime
	fInternalTime = t;
	

	DBGL;
}


//________________________________________________________________________
void TreatmentFactory::Dump()
{
//------ Cooling ------//
	for(int i = (int)fCoolingEndOfCycle.size()-1; i >=0 ; i--)	// IV End Of Cooling
	{

		int idx = fCoolingEndOfCycle.at(i);			// Get Index number
		pair<IsotopicVector, IsotopicVector> SeparatedIV = Separation(fIVCooling.at(idx));// Separation
		fIVWaste += SeparatedIV.second;				// Add to waste
		AddIVSeparating(SeparatedIV.first, fInternalTime );	// Add to speration
		
		fCoolingEndOfCycle.erase(fCoolingEndOfCycle.begin()+i);	// Remove index entry
		RemoveIVCooling(idx);					// Remove IVcooling
	}
	if((int)fCoolingEndOfCycle.size() != 0 )// Control
	{
		cout << "Problem while Dumping Cooling"<< endl;
		exit (1);
	}



//------ Separation ------//
	for(int i = (int)fSeparationEndOfCycle.size()-1; i >=0 ; i--)	// IV End Of Cooling
	{
		int idx = fSeparationEndOfCycle.at(i);			// Get Index number
		AddIVStock(fIVSeparating.at(idx));
		
		RemoveIVSeparation(idx);				// Remove Separated IV
		fSeparationEndOfCycle.erase(fSeparationEndOfCycle.begin()+i);	// Remove index entry

	}

	if((int)fSeparationEndOfCycle.size() != 0 ) // Control
	{
		cout << "Problem while Dumping Separtion"<< endl;
		exit (1);
	}

}


//________________________________________________________________________
void TreatmentFactory::Write(string TFfilename)
{
	DBGL;
	
	for (int i = 0; i < (int)fIVStock.size(); i++)
	{
		ostringstream ostr;
		ostr << i ;
		string TFStockFilename = TFfilename + "_ST" + ostr.str();
		fIVStock.at(i).Write(TFStockFilename, fInternalTime);
	}

	

	for (int i = 0; i < (int)fIVSeparating.size(); i++)
	{
		ostringstream ostr;
		ostr << fSeparatingIndex.at(i) ;
		string TFSeparationFilename = TFfilename + "_TR" + ostr.str();
		fIVSeparating.at(i).Write(TFSeparationFilename, fInternalTime);
		
	}


	for (int i = 0; i < (int)fIVCooling.size(); i++)
	{
		ostringstream ostr;
		ostr << fCoolingIndex.at(i) ;
		string TFCoolingFilename = TFfilename + "_CO" + ostr.str();
		fIVCooling.at(i).Write(TFCoolingFilename, fInternalTime);
	}
	

	string TFFullStockFilename = TFfilename + "_ST" ;
	
	fIVFullStock.Write(TFFullStockFilename, fInternalTime);


	string TFWasteFilename = TFfilename + "_WASTE";
	fIVWaste.Write(TFWasteFilename, fInternalTime);
		

	
	
	DBGL;
}




