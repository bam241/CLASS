#ifndef __Pool_HXX__
#define __Pool_HXX__


/*!
 \file
 \brief Header file for Pool class.
 
 The aim of the Class is to manage evolution of all out reactor fuel. from Cooling to Waste or storage
 
 
 @author BaM
 @version 2.0
 */




#include <string>
#include <map>

#include "CLSSFacility.hxx"
#include "IsotopicVector.hxx"

using namespace std;
typedef long long int cSecond;

class Storage;
class CLASS;
class LogFile;
template <class T> 
class DataBank;



class Pool : public CLSSFacility
{
public :
	///< Normal constructor
 	Pool();
 	Pool(LogFile* log);
	///< Advanced Constructor
	Pool(LogFile* log, double abstime,
			 double coolingtime = 5*3600.*24.*365.25); //!<
	
	Pool(LogFile* log, Storage* Storage,
			 double abstime = 0,
			 double coolingtime = 5*3600.*24.*365.25); //!<

	///< Normal Destructor.
	~Pool();

  


//********* Set Method *********//
	void SetStorage(Storage* storage)	{ fStorage = storage; fPutToWaste = true; }		//!< Set the Pointer to the Storage
	void SetPutToWaste(bool val)		{ fPutToWaste = val; }		//!< Set True if IV goes to waste after cooling false instead

	void SetDecayDataBase(DataBank<ZAI>* ddb)	{ fDecayDataBase = ddb; }		//!< Set the pointer to the Decay DataBase

	void SetCoolingTime(double time) 		{ SetCycleTime((cSecond)time); }			//!< Set Cooling Time

//********* Get Method *********//
	Storage*	GetStorage()	const	{ return fStorage; }		//!< Return the Pointer to the Storage
	bool		GetPutToWaste()	const	{ return fPutToWaste; }		//!< Return True if IV goes to waste after cooling false instead

	
	
	cSecond GetCoolingTime() const		{ return GetCycleTime(); }		//!< Return the Cooling Time

	
	DataBank<ZAI>* 	GeDecayDataBase() const	{ return fDecayDataBase; }	//!< Return the pointer to the Decay DataBase



//********* IsotopicVector Method *********//

//--------- Cooling ---------//
	vector<cSecond>	GetCoolingStartingTime() const		{ return fCoolingStartingTime; }
											//!< Return the vector of Cooling Sstarting Time
	vector<IsotopicVector>	GetIVCooling() const		{ return fIVCooling; }	//!< Return the vector of Cooling IsotopicVector
	void			AddIVCooling(IsotopicVector IV);			//!< Add Cooling IsotopicVector
	void			RemoveIVCooling(int i);					//!< Remove a Cooling IsotopicVector
	IsotopicVector		GetFullCooling()		{return GetInsideIV(); }


//********* Other Method *********//
	void Evolution(cSecond t);		//!< Performe the evolution until the Time t
	void Dump();				//!< Write Modification (exchange between Cooling, Separation and Storage)
	
	
protected :
	
	
	
//********* Internal Parameter *********//
	Storage*		fStorage;		//!< Pointer to the Stock
	bool			fPutToWaste;		//!< True if IV goes to waste after cooling false instead
							//	LogFile*		fLog;			//!< Pointer to the Log


	DataBank<ZAI>*	fDecayDataBase;		//!< Pointer to the Decay DataBase


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




	ClassDef(Pool,2);
};

#endif
