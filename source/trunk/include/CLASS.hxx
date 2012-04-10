#ifndef _CLASS_HXX_
#define _CLASS_HXX_

/*!
 \file 
 \brief Header file for CLASS class. 
*/

#include "IsotopicVector.hxx"
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <string>
#include <map>


using namespace std;

class EvolutionDataBase;
class Reactor;
class TreatmentFactory;


class CLASS 
{
public :
	//! Normal Constructor.
 	CLASS();
 	CLASS(long int abstime);
  	CLASS(string name);
  	CLASS(string name, long int abstime);
 	
 	//! Normal Destructor.
 	~CLASS();
	

//********* Get Method *********//
	long int			GetAbsoluteTime()	{return fAbsoluteTime;}		//!< Return the Absolute Clock
 	map<long int, int>		GetTimeStep()		{return fTimeStep;}		//!< Return the Time Step Vector
 	vector<TreatmentFactory*>	GetTreatmentFactory()	{return fTreatmentFactory;}	//!< Return the TF Vector
 	vector<Reactor*>		GetReactor()		{return fReactor;}		//!< Return the Reactor Vector
	long int			GetPrintSet()		{return fPrintStep;}		//!< Return the Print Step Periodicity
	bool				GetStockManagement()	{return fStockManagement;}
	string				GetOutputName()		{return fOutputName;}

//********* Set Method *********//
	void	SetTimeStep(int timestep) 	{fPrintStep = timestep;}	//!< Set the Printing Step periodicity
	void	SetStockManagement(bool val)	{fStockManagement = val;}
	void	SetOutpurName(string name)	{fOutputName =name;}

//********* Modification Method *********//
	IsotopicVector	BuildIsotopicVector(IsotopicVector isotopicvector );	//!< Build The needed Isotopic Vector from the stock
 	void	AddTreatmentFactory(TreatmentFactory* treatmentfactory);	//!< Add A TF to the Park
 	void	AddReactor(Reactor* reactor);					//!< Add a Reactor to the Park 
	
//********* Evolution Method *********//
	void	BuildTimeVector(long int t);				//!< Build the Time Evolution Vector
	void	Evolution(long int t);					//!< Do the Evolution
	void	TreatmentEvolution();					//!< Do TF Evolution
	void	ReactorEvolution();					//!< Do the Reactor Evolution
	void	RemoveReactor();					//!<  

	void	AddTotalWaste(IsotopicVector IV)		{fTotalWaste +=IV;
								 fIVTotal +=IV;}

	void	AddTotalStock(IsotopicVector IV)		{fTotalStock +=IV;
								 fIVInCycleTotal +=IV;
								 fIVTotal +=IV;}
	void	RemoveTotalStock(IsotopicVector IV)		{fTotalStock -=IV;
								 fIVInCycleTotal -=IV;
								 fIVTotal -=IV;}
	void	AddTotalCooling(IsotopicVector IV)		{fTotalCooling +=IV;
								 fIVInCycleTotal +=IV;
								 fIVTotal +=IV;}
	void	AddTotalSeparating(IsotopicVector IV)		{fTotalSeparating +=IV;
								 fIVInCycleTotal +=IV;
								 fIVTotal +=IV;}
	void	AddTotalInReactor(IsotopicVector IV)		{fTotalInReactor +=IV;
								 fIVInCycleTotal +=IV;
								 fIVTotal +=IV;}

	void	AddIVInCycleTotal(IsotopicVector IV)		{fIVInCycleTotal +=IV;}
	void	AddIVTotal(IsotopicVector IV)			{fIVTotal +=IV;}
	void	AddTotalGodIncome(IsotopicVector IV)		{fTotalGodIncome +=IV;}

//********* Other Method *********//
	void	Print();
	void	Write();
	
	void	OpenOutputTree();
	void	CloseOutputTree();
	void	OutAttach();
	void	ResetQuantity();
	
protected :
 	long int			fPrintStep;
	long int			fAbsoluteTime;		//!< Absolute Clock
 	map<long int, int>		fTimeStep;		//!< Time Step Vector for the evolution : 
 								//!< 1 to print, 
 								//!< 2 reactor destruction, 4 end of reactor cycle,
 								//!< 8 end of Cooling, 16 end of Separation. 
 	vector<TreatmentFactory*>	fTreatmentFactory;		//!< Vector of Treament Factory
 	vector<int>			fTreatmentFactoryIndex;		//!< Vector of TF Index
 	int				fTreatmentFactoryLastIndex;	//!< Quantity of Factory Added
 	vector<Reactor*>		fReactor;			//!< Vector of Reactor
 	vector<int>			fReactorIndex;			//!< Vecotr of Reactor Index
	int				fReactorLastIndex;		//!< Quantity of Reactor Added

	bool				fStockManagement;		//!< ...
	
	string				fOutputName;
	TFile				*fOutTree;
	TTree				*fOutT;

	IsotopicVector			fTotalWaste;
	IsotopicVector			fTotalStock;
	IsotopicVector			fTotalGodIncome;
	IsotopicVector			fTotalCooling;
	IsotopicVector			fTotalSeparating;
	
	IsotopicVector			fTotalInReactor;
	
	IsotopicVector			fIVInCycleTotal;
	IsotopicVector			fIVTotal;


};


#endif
