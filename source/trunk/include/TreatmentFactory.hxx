#ifndef __TreatmentFactory_HXX__
#define __TreatmentFactory_HXX__

/*!
 \file 
 \brief Header file for TreatmentFactory class.
*/


#include "EvolutionDataBase.hxx"
#include "CLASS.hxx"
#include "IsotopicVector.hxx"
#include "LogFile.hxx"

#include <vector>

using namespace std;



class TreatmentFactory : public TObject
{
public :
 	TreatmentFactory();						///< Normal constructor

	TreatmentFactory(long int abstime,
			 long int coolingtime = 5,
			 long int separationtime = 2,
			 EvolutionDataBase* evolutivedb = NULL);		///< Advanced Constructor
	
 	~TreatmentFactory();						///< Normal Destructor.




//********* Set Method *********//
	void SetParc(CLASS* parc)			{fParc = parc;}				//!< Set the Pointer to the Parc
	void SetLog(LogFile* Log)			{fLog = Log;}				//!< Set the Pointer to the Log

	void SetDecayDataBase(EvolutionDataBase* ddb)	{fDecayDataBase = ddb;}			//!< Set the pointer to the Evolution DataBase
	void SetIVWaste(IsotopicVector IV) 		{fIVWaste = IV;}			//!< Set the initial composition of the Ultimate Waste IV
	void SetCoolingTime(long int time) 		{fCoolingTime = time;}			//!< Set Cooling Time
	void SetSeparationTime(long int time) 		{fSeparationTime = time;}		//!< Set Separation Time
	void SetStockManagement(bool bmanagement)	{fStockManagement = bmanagement;}	//!< Set the way to manage the stock
												//!< True  : to treat each IV stock individualy
												//!< False : to treat the sotck as a Soup 

//********* Get Method *********//
	LogFile*		GetLog()		{return fLog;}			//!< Return the Pointer to the Log
	CLASS*			GetParc()		{return fParc;}			//!< return the Pointer to the Parc
	
	long int GetInternalTime() const		{return fInternalTime;}		//!< Return Creation Time
	long int GetCreationTime() const		{return fCreationTime;}		//!< Return Internal Time
	long int GetCoolingTime() const			{return fCoolingTime;}		//!< Return the Cooling Time
	long int GetSeparationTime() const		{return fSeparationTime;}	//!< Return the Separation Time
	
	EvolutionDataBase* 	GeDecayDataBase() const	{return fDecayDataBase;}	
								//!< Return the pointer to the Evolution DataBase
	
	map<ZAI ,double>	GetValorisableIV() const {return fValorisableIV;} 	///< Return the Valorisable Table
	void			AddValorisableIV(ZAI zai, double factor);		///< Add Valorisable Element



//********* IsotopicVector Method *********//

//--------- Cooling ---------//
	vector<long int>	GetCoolingStartingTime() const	{return fCoolingStartingTime;}
											//!< Return the vector of Cooling Sstarting Time
	vector<IsotopicVector>	GetIVCooling() const		{return fIVCooling;}	//!< Return the vector of Cooling IsotopicVector
	void			AddIVCooling(IsotopicVector IV);			//!< Add Cooling IsotopicVector
	void			RemoveIVCooling(int i);					//!< Remove a Cooling IsotopicVector


//---------- Waste ----------//
	IsotopicVector		GetIVWaste() const			{return fIVWaste;}		//!< Return the Waste IsotopicVector
	void			AddIVWaste(ZAI zai, double quantity)	{fIVWaste += zai*quantity; }	//!< Add a ZAI*quantity to the waste 
	void			AddIVWaste(IsotopicVector IV)		{fIVWaste += IV; }		//!< Add a IsotopicVector to the waste


//------- Separation --------//
	vector<long int>	GetSeparatingStartingTime() const	{return fSeparatingStartingTime;}	//!< Return the vector of Separation Starting Time
	vector<IsotopicVector>	GetIVSeparating() const			{return fIVSeparating;}			//!< Return the vector of Separation IsotopicVector
	
	void			AddIVSeparating(IsotopicVector IV);						//!< Add Separation IsotopicVector
	void			AddIVSeparating(IsotopicVector IV, long int absolutadditiontime);		//!< Add Separation IsotopicVector at an absolute time
	void			RemoveIVSeparation(int i);							//!< Remove a Treated IsotopicVector

//---------- Stock ----------//
	vector<IsotopicVector>	GetIVStock() const		{return fIVStock;}			//!< Return the Stock IsotopicVector
	void SetIVStock(vector<IsotopicVector> IVsStock)	{fIVStock =IVsStock;}			//!< Set The Stock isotopicVector
	void ClearIVStock();
	
	void AddIVStock(ZAI zai, double quantity)		{AddIVStock(zai*quantity);}		//!< Add a ZAI*quantity to the stock
	void AddIVStock(IsotopicVector isotopicvector)		{fIVStock.push_back(isotopicvector); fIVFullStock+= isotopicvector; fParc->AddTotalStock(isotopicvector);}
													//!< Add an Isotopicvector to the stock
	IsotopicVector		GetIVFullStock() const		{return fIVFullStock;}			//!< Return the Full Stock
	void AddIVFullStock(ZAI zai, double quantity)		{fIVFullStock.Add(zai*quantity); }	//!< Add a ZAI*quantity to the stock
	void AddIVFullStock(IsotopicVector isotopicvector)	{fIVFullStock.Add(isotopicvector); }	//!< Add a IsotopicVector to the stock

	void TakeFromStock(IsotopicVector isotopicvector, int index); 		//!< Take isotopicvector from the (index)st vector of the stock 	
	
//-------- GodIncome --------//

	IsotopicVector		GetIVGodIncome() const		{return fIVGodIncome;}	//!< Return the God Providings IsotopicVector
	void AddIVGodIncome(ZAI zai, double quantity)		{AddIVGodIncome(zai*quantity);}	//!< Add a ZAI*quantity to GodIncome
	void AddIVGodIncome(IsotopicVector isotopicvector)	{fIVGodIncome.Add(isotopicvector); fParc->AddTotalGodIncome(isotopicvector);}	//!< Add a isotopicVector to GodIncome



//********* Other Method *********//
	void Evolution(long int t);		//!< Performe the evolution until the Time t
	void Dump();				//!< Write Modification (exchange between Cooling, Separation and Stock)
	
	
protected :
	
	long int 		fInternalTime;		///< Internal Clock
	bool 			IsStarted;		///< True if Running, False Otherwise
	
//********* Internal Parameter *********//
	CLASS* 			fParc;			//!< Pointer to the Parc
	LogFile*		fLog;			//!< Pointer to the Log

	map<ZAI ,double>	fValorisableIV;		///< The Valorisable Table
	EvolutionDataBase*	fDecayDataBase;		//!< Pointer to the Evolution DataBase

	long int 		fCreationTime;		///< Date of Creation of the Factory
	long int 		fCoolingTime;		///< Cooling Duration Time
	long int 		fSeparationTime;	///< Separation Duration Time


//********* Isotopic Quantity *********//
//--------- Cooling ---------//
	vector<IsotopicVector>	fIVCooling;		///< Vector of the Cooling Isotopic Vector
	vector<long int>	fCoolingStartingTime;	///< Vector of the Cooling Starting Time
	vector<int>		fCoolingIndex;		///< Vector of the Cooling Index
	int			fCoolingLastIndex;	//!< Number of Cooling IV Treated
	vector<int>		fCoolingEndOfCycle;	//!< Index of the Cooling IV reaching the End of a Cooling Cycle

//------- Separation --------//
	vector<IsotopicVector>	fIVSeparating;		///< Vector of the Treated Isotopic Vector
	vector<long int>	fSeparatingStartingTime;///< Vector of the Treated Starting Time
	vector<int>		fSeparatingIndex;	///< Vector of the Treated Index
	int			fSeparatingLastIndex;	//!< Number of Separation IV Treated
	vector<int>		fSeparationEndOfCycle;	//!< Index of the Separation IV reaching the End of a Cooling Cycle
	
//---------- Stock ----------//
	vector<IsotopicVector>	fIVStock;		///< The Stock IsotopicVector
	IsotopicVector		fIVGodIncome;		///< Stock quantity providing by a godlike entitie
	IsotopicVector		fIVFullStock;		///< Full stock conglomerat
	bool 			fStockManagement;	///< if true use the stock as the full conglomerat, not individualy

	IsotopicVector		fIVWaste;		///< The Ultimate Waste IsotopicVector 

//********* Private Method *********//
	IsotopicVector GetDecayProduct(IsotopicVector isotopicvector, long int t);	//!< Get IsotopicVector Decay at the t time
	void	CoolingEvolution(long int t);						//!< Deal the cooling and then send it to Separation
	void	SeparatingEvolution(long int t);					//!< Deal the Separating IV Decay Evolution and then send it to stock 
	void	StockDecay(long int t);							//!< Deal the Stock Decay Evolution


	void	WasteDecay(long int t);							//!< Deal the Waste Decay Evolution
	pair<IsotopicVector, IsotopicVector> Separation(IsotopicVector isotopicvector);	//!< Make the Separation 
											//!< return IV[0] -> To Stock / IV[1] -> To Waste



	ClassDef(TreatmentFactory,3);
};

#endif
