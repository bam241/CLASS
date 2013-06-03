#ifndef __TreatmentFactory_HXX__
#define __TreatmentFactory_HXX__




#include "TObject.h"
#include <string>
#include <map>
#include "IsotopicVector.hxx"

using namespace std;
typedef long long int cSecond;

class Storage;
class CLASS;
class LogFile;
template <class T> 
class EvolutionDataBase;


/*!
 \file
 \brief Header file for TreatmentFactory class.
 
 The aim of the Class is to manage evolution of all out reactor fuel. from Cooling to Waste
 
 
 @author BaM
 @version 0.
 */



class TreatmentFactory : public TObject
{
public :
	///< Normal constructor
 	TreatmentFactory();
	///< Advanced Constructor
	TreatmentFactory(double abstime,
			 double coolingtime = 5*3600.*24.*365.25); //!<
	
	TreatmentFactory(Storage* Storage,
			 double abstime = 0,
			 double coolingtime = 5*3600.*24.*365.25); //!<

	///< Normal Destructor.
	~TreatmentFactory();

  


//********* Set Method *********//
	void SetId(int id)			{ fId = id; }			//!< Set The TF Parc'Id
	void SetParc(CLASS* parc)		{ fParc = parc; }		//!< Set the Pointer to the Parc
	void SetLog(LogFile* Log)		{ fLog = Log; }			//!< Set the Pointer to the Log
	void SetStorage(Storage* storage)	{ fStorage = storage; fPutToWaste = true; }		//!< Set the Pointer to the Storage
	void SetPutToWaste(bool val)		{ fPutToWaste = val; }		//!< Set True if IV goes to waste after cooling false instead

	void SetDecayDataBase(EvolutionDataBase<ZAI>* ddb)	{ fDecayDataBase = ddb; }		//!< Set the pointer to the Decay DataBase

	void SetCoolingTime(double time) 		{ fCoolingTime = (cSecond)time; }			//!< Set Cooling Time

//********* Get Method *********//
	int 		GetId()		const	{ return fId; }			//!< Return the TF Parc'Is
	LogFile*	GetLog()	const	{ return fLog; }		//!< Return the Pointer to the Log
	CLASS*		GetParc()	const	{ return fParc; }		//!< Return the Pointer to the Parc
	Storage*	GetStorage()	const	{ return fStorage; }		//!< Return the Pointer to the Storage
	bool		GetPutToWaste()	const	{ return fPutToWaste; }		//!< Return True if IV goes to waste after cooling false instead

	
	
	cSecond GetInternalTime() const			{ return fInternalTime; }		//!< Return Creation Time
	cSecond GetCreationTime() const			{ return fCreationTime; }		//!< Return Internal Time
	cSecond GetCoolingTime() const			{ return fCoolingTime; }		//!< Return the Cooling Time

	
	EvolutionDataBase<ZAI>* 	GeDecayDataBase() const	{ return fDecayDataBase; }	//!< Return the pointer to the Decay DataBase



//********* IsotopicVector Method *********//

//--------- Cooling ---------//
	vector<cSecond>	GetCoolingStartingTime() const	{ return fCoolingStartingTime; }
											//!< Return the vector of Cooling Sstarting Time
	vector<IsotopicVector>	GetIVCooling() const		{ return fIVCooling; }	//!< Return the vector of Cooling IsotopicVector
	void			AddIVCooling(IsotopicVector IV);			//!< Add Cooling IsotopicVector
	void			RemoveIVCooling(int i);					//!< Remove a Cooling IsotopicVector
	IsotopicVector		GetFullCooling();


//********* Other Method *********//
	void Evolution(cSecond t);		//!< Performe the evolution until the Time t
	void Dump();				//!< Write Modification (exchange between Cooling, Separation and Storage)
	
	
protected :
	int		fId;			//!< Identity of the Reactor inside the Parc
	cSecond		fInternalTime;		///< Internal Clock
	bool 		IsStarted;		///< True if Running, False Otherwise
	
	
	
//********* Internal Parameter *********//
	CLASS* 			fParc;			//!< Pointer to the Parc
	Storage*		fStorage;		//!< Pointer to the Stock
	bool			fPutToWaste;		//!< True if IV goes to waste after cooling false instead
	LogFile*		fLog;			//!< Pointer to the Log


	EvolutionDataBase<ZAI>*	fDecayDataBase;		//!< Pointer to the Decay DataBase

	cSecond 		fCreationTime;		///< Date of Creation of the Factory
	cSecond 		fCoolingTime;		///< Cooling Duration Time


//********* Isotopic Quantity *********//
//--------- Cooling ---------//
	vector<IsotopicVector>	fIVCooling;		///< Vector of the Cooling Isotopic Vector
	vector<cSecond>		fCoolingStartingTime;	///< Vector of the Cooling Starting Time
	vector<int>		fCoolingIndex;		///< Vector of the Cooling Index
	int			fCoolingLastIndex;	//!< Number of Cooling IV Treated
	vector<int>		fCoolingEndOfCycle;	//!< Index of the Cooling IV reaching the End of a Cooling Cycle


//********* Private Method *********//
	IsotopicVector GetDecay(IsotopicVector isotopicvector, cSecond t);	//!< Get IsotopicVector Decay at the t time
	void	CoolingEvolution(cSecond t);					//!< Deal the cooling and then send it to Separation




	ClassDef(TreatmentFactory,1);
};

#endif
