#ifndef __Storage_HXX__
#define __Storage_HXX__

/*!
 \file 
 \brief Header file for Storage class.

 The aim of the Class is to manage evolution of Storage

 
 @author BaM
 @version 0.
 */

#include "IsotopicVector.hxx"


#include <vector>

using namespace std;
typedef long long int cSecond;

class CLASS;
class LogFile;
template <class T> 
class EvolutionDataBase;


class Storage : public TObject
{
public :
	///< Normal Constructor.
 	Storage();
 	///< Advanced Constructor
	Storage(EvolutionDataBase<ZAI>* evolutivedb);						//!<

 	///< Normal Destructor.
 	~Storage();



//********* Set Method *********//
	void SetId(int id)				{ fId = id; }				//!< Set The Storage Parc'Id
	void SetParc(CLASS* parc)			{ fParc = parc; }			//!< Set the Pointer to the Parc
	void SetLog(LogFile* Log)			{ fLog = Log; }				//!< Set the Pointer to the Log

	void SetDecayDataBase(EvolutionDataBase<ZAI>* ddb)	{ fDecayDataBase = ddb; }	//!< Set the pointer to the Decay DataBase
	void SetStock(vector<IsotopicVector> IVsstock)	{ fIVStock = IVsstock; }		//!< Set The Storage isotopicVector


//********* Get Method *********//
	int 			GetId()		const		{ return fId; }			//!< Return the Storage Parc'Is
	LogFile*		GetLog()	const		{ return fLog; }		//!< Return the Pointer to the Log
	CLASS*			GetParc()	const		{ return fParc; }		//!< return the Pointer to the Parc

	double GetInternalTime() const				{ return fInternalTime; }	//!< Return Creation Time
	
	//!<
	EvolutionDataBase<ZAI>* GeDecayDataBase() const		{ return fDecayDataBase; }	//!< Return the pointer to the Decay DataBase
	vector<IsotopicVector> GetStock()	const		{ return fIVStock; }		//!< Return the Storage IsotopicVector
	IsotopicVector GetFullStock()		const		{ return fIVFullStock; }	//!< Return the Full Storage

//********* IsotopicVector Method *********//

//---------- Storage ----------//


	void ClearStock();
	
	void AddToStock(ZAI zai, double quantity)		{ AddToStock(zai*quantity); }		//!< Add a ZAI*quantity to the Storage
	void AddToStock(IsotopicVector isotopicvector);							//!< Add an Isotopicvector to the Storage
	void AddToFullStock(ZAI zai, double quantity)		{ fIVFullStock += zai*quantity; }	//!< Add a ZAI*quantity to the Storage
	void AddToFullStock(IsotopicVector isotopicvector)	{ fIVFullStock += isotopicvector; }	//!< Add a IsotopicVector to the Storage

	void TakeFractionFromStock(int IVId,double fraction);						//!< Take a part from an IV in sotck;
	void TakeFromStock(IsotopicVector isotopicvector);						//!<

//********* Other Method *********//
	void Evolution(cSecond t);		//!< Performe the evolution until the Time t
	
	
protected :
	int		fId;			//!< Identity of the Reactor inside the Parc
	cSecond		fInternalTime;		///< Internal Clock
	
//********* Internal Parameter *********//
	CLASS* 		fParc;			//!< Pointer to the Parc
	LogFile*	fLog;			//!< Pointer to the Log

	EvolutionDataBase<ZAI>*	fDecayDataBase;	//!< Pointer to the Decay DataBase



//********* Isotopic Quantity *********//

//---------- Storage ----------//
	vector<IsotopicVector>	fIVStock;	///< The Storage IsotopicVector
	IsotopicVector		fIVFullStock;	///< Full Storage conglomerat



//********* Private Method *********//
	IsotopicVector GetDecay(IsotopicVector isotopicvector, cSecond t);	//!< Get IsotopicVector Decay at the t time
	void	StorageEvolution(cSecond t);					//!< Deal the Storage Decay Evolution




	ClassDef(Storage,1);
};

#endif
