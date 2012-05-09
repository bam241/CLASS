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
#include <iostream>


using namespace std;

class EvolutionDataBase;
class Reactor;
class TreatmentFactory;
class LogFile;


class CLASS 
{
public :
	///< Normal Constructor.
 	CLASS();
 	CLASS(long int abstime);
  	CLASS(string name);
  	CLASS(string name, long int abstime);
 	
 	///< Normal Destructor.
 	~CLASS();
	

//********* Get Method *********//
	long int			GetAbsoluteTime()	{return fAbsoluteTime;}		///< Return the Absolute Clock
 	map<long int, int>		GetTimeStep()		{return fTimeStep;}		///< Return the Time Step Vector
 	vector<TreatmentFactory*>	GetTreatmentFactory()	{return fTreatmentFactory;}	///< Return the TF Vector
 	vector<Reactor*>		GetReactor()		{return fReactor;}		///< Return the Reactor Vector
	long int			GetPrintSet()		{return fPrintStep;}		///< Return the Print Step Periodicity
	bool				GetStockManagement()	{return fStockManagement;}
	string				GetOutputName()		{return fOutputName;}
	LogFile*			GetLog()		{return fLog;}
//********* Set Method *********//
	void	SetTimeStep(int timestep) 	{fPrintStep = timestep;}	///< Set the Printing Step periodicity
	void	SetStockManagement(bool val)	{fStockManagement = val;}
	void	SetBuildingMethod(int val)	{fStockManagement = true; fBuildingMethod = val;}
	void	SetOutpurName(string name)	{fOutputName =name;}

//********* Add Method *********//
	void	AddTreatmentFactory(TreatmentFactory* treatmentfactory);	///< Add A TF to the Park
 	void	AddReactor(Reactor* reactor);					///< Add a Reactor to the Park 

	
//********* IV Creation Method *********//
	IsotopicVector	BuildIsotopicVector(IsotopicVector isotopicvector );	///< Build The needed Isotopic Vector from the stock
	IsotopicVector	MonteCarloBuild(IsotopicVector isotopicvector );	///< Build the needed IV with the MonteCarlo Method
	IsotopicVector	MinimizationBuild(IsotopicVector isotopicvector );	///< Build the needed IV with the Minimization Method

	void	UpdateParcStock();
	void	DumpParcStock();
	
	
	
	
	
 	
//********* Evolution Method *********//
	void	BuildTimeVector(long int t);				///< Build the Time Evolution Vector
	void	Evolution(long int t);					///< Do the Evolution
	void	TreatmentEvolution();					///< Do TF Evolution
	void	ReactorEvolution();					///< Do the Reactor Evolution
	void	RemoveReactor();					///<  


	void	AddTotalWaste(IsotopicVector IV)		{fIVTotal +=IV; fTotalWaste +=IV;}
	void	AddTotalStock(IsotopicVector IV)		{fIVTotal +=IV; fIVInCycleTotal +=IV; fTotalStock +=IV;}
	void	RemoveTotalStock(IsotopicVector IV)		{fIVTotal -=IV; fIVInCycleTotal -=IV; fTotalStock -=IV;}
	void	AddTotalCooling(IsotopicVector IV)		{fIVTotal +=IV; fIVInCycleTotal +=IV; fTotalCooling +=IV;}
	void	RemoveTotalCooling(IsotopicVector IV)		{fIVTotal -=IV; fIVInCycleTotal -=IV; fTotalCooling -=IV;}
	void	AddTotalSeparating(IsotopicVector IV)		{fIVTotal +=IV; fIVInCycleTotal +=IV; fTotalSeparating +=IV;}
	void	RemoveTotalSeparating(IsotopicVector IV)	{fIVTotal -=IV; fIVInCycleTotal -=IV; fTotalSeparating -=IV;}
	void	AddTotalInReactor(IsotopicVector IV)		{fIVTotal +=IV; fIVInCycleTotal +=IV; fTotalInReactor +=IV;}
	void	RemoveTotalInReactor(IsotopicVector IV)		{fIVTotal -=IV; fIVInCycleTotal -=IV; fTotalInReactor -=IV;}
	void	AddIVInCycleTotal(IsotopicVector IV)		{fIVTotal +=IV; fIVInCycleTotal +=IV;}
	void	AddIVTotal(IsotopicVector IV)			{fIVTotal +=IV;}
	void	AddTotalGodIncome(IsotopicVector IV)		{fTotalGodIncome +=IV;}

//********* In/Out related Method *********//
	void	ProgressPrintout(int starttime, long int t);
	
	void	Print();
	void	Write();
	
	void	OpenOutputTree();
	void	CloseOutputTree();
	void	OutAttach();
	void	ResetQuantity();
	
	
	
	
protected :
 	long int			fPrintStep;
	long int			fAbsoluteTime;		///< Absolute Clock
 	map<long int, int>		fTimeStep;		///< Time Step Vector for the evolution : 
 								///< 1 to print, 
 								///< 2 reactor destruction, 4 end of reactor cycle,
 								///< 8 end of Cooling, 16 end of Separation. 
 	vector<TreatmentFactory*>	fTreatmentFactory;		///< Vector of Treament Factory
 	vector<int>			fTreatmentFactoryIndex;		///< Vector of TF Index
 	int				fTreatmentFactoryLastIndex;	///< Quantity of Factory Added
 	vector<Reactor*>		fReactor;			///< Vector of Reactor
 	vector<int>			fReactorIndex;			///< Vecotr of Reactor Index
	int				fReactorLastIndex;		///< Quantity of Reactor Added

	bool				fStockManagement;		///< ...
	int				fBuildingMethod;
	
	
	LogFile*			fLog;
	
	string				fOutputName;
	TFile*				fOutTree;
	TTree*				fOutT;

	IsotopicVector			fTotalWaste;
	IsotopicVector			fTotalStock;
	IsotopicVector			fTotalGodIncome;
	IsotopicVector			fTotalCooling;
	IsotopicVector			fTotalSeparating;
	
	IsotopicVector			fTotalInReactor;
	
	IsotopicVector			fIVInCycleTotal;
	IsotopicVector			fIVTotal;
	
	map< pair<int,int>, IsotopicVector >	fParcStock;


};


#endif
