#ifndef __Pool_HXX__
#define __Pool_HXX__
/*!
 \file
 \brief Header file for Pool class.
 */

#include <string>
#include <map>

#include "CLASSConstante.hxx"
#include "CLASSBackEnd.hxx"
#include "IsotopicVector.hxx"

using namespace std;
typedef long long int cSecond;

class CLASSBackEnd;
class CLASSLogger;
class DecayDataBank;

//-----------------------------------------------------------------------------//
//! Defines the spent fuel pool

/*!
 This class deal with the management of the spent fuel pool


 @author BaM
 @version 2.0
 */
//________________________________________________________________________



class Pool : public CLASSBackEnd
{
public :


//********* Constructor/Destructor Method *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	
	Pool();				///< Normal Constructor.

	//{
	/// Special Constructor.
	/*!
	 Make a new Pool
	 \param log : used for the log.
	 \param coolingtime : duration of the cooling (default : 5 years) in [s].
	 */
	Pool(CLASSLogger* Log, cSecond coolingtime = 5*cYear); //!<
	//}


	//{
	/// Special  Constructor.
	/*!
	 Make a new Pool
	 \param log : used for the log.
	 \param backend : CLASSBackend which get the fuel after the cooling
	 \param coolingtime : duration of the cooling (default : 5 years) in [s].
	 */
	Pool(CLASSLogger* log, CLASSBackEnd* backend,
			 cSecond coolingtime = 5*cYear); //!<
	//}


	~Pool();	///< Normal Destructor.
	//@}




//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{

	void SetOutBackEndFacility(CLASSBackEnd* befacility)
						{  fOutBackEndFacility = befacility;
						   SetIsStorageType(false);
						   fPutToWaste = false; }		//!< Set the pointer to facility at the back end of the pool

	void SetPutToWaste(bool val)		{ fPutToWaste = val; }		//!< Set true if IV goes to waste after cooling false instead

	void SetIVArray(vector<IsotopicVector> ivarray);			//! not use there (does nothing)
	void SetIVArray(vector<IsotopicVector> ivarray, vector<cSecond> timearray); //!< Set the IsotopicVector Array at the corresponding time


	using CLASSBackEnd::SetName;

	//@}




//********* Get Method *********//

	/*!
	 \name Get Method
	 */
	//@{

	bool		GetPutToWaste()	const	{ return fPutToWaste; }		//!< Return true if IV goes to waste after cooling, false instead

	//@}




//********* IsotopicVector Managment Method *********//

	/*!
	 \name IsotopicVector Managment Method
	 */
	//@{

	vector<cSecond>	GetCoolingStartingTime() const
						{ return GetIVArrayArrivalTime(); }	//!< Return vector of the arrival time of each IV in the Pool
	void	RemoveIVCooling(int i);					//!< Remove a IsotopicVector from cooling

	void	AddIV(IsotopicVector isotopicvector);			//!< Add an Isotopicvector to the IVArray
	//@}




//********* Other Method *********//

	//@}
	/*!
	 \name Other Method
	 */
	//@{

	void Evolution(cSecond t);		//!< Perform the evolution until time t
	void Dump();				//!< Write modification (exchange between Cooling, Separation and Storage)
	
	//@}

protected :
	
	
	
//********* Internal Parameter *********//
	bool			fPutToWaste;	//!< True if IV goes to waste after cooling false instead


//********* Isotopic Quantity *********//
//--------- Cooling ---------//
	vector<int>		fCoolingIndex;		///< Vector of the cooling index
	int			fCoolingLastIndex;	//!< Number of cooling IV handle
	vector<int>		fCoolingEndOfCycle;	//!< Index of the cooling IV reaching the end of a cooling cycle


//********* Private Method *********//
	void	CoolingEvolution(cSecond t);					//!< Deal the cooling evolution




	ClassDef(Pool,3);
};

#endif
