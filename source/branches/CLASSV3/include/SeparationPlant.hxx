#ifndef __SeparationPlant_HXX__
#define __SeparationPlant_HXX__
/*!
 \file
 \brief Header file for SeparationPlant class.
 */

#include <string>
#include <map>

#include "CLASSBackEnd.hxx"
#include "Storage.hxx"
#include "IsotopicVector.hxx"

using namespace std;
typedef long long int cSecond;

class CLASSBackEnd;
class CLASSLogger;
class DecayDataBank;

//-----------------------------------------------------------------------------//
/*!
 Define a SeparationPlant.
 The aim of the Class is to separate an IV into several IV (MA, Pu, PF, etc...) and to send it to corresponding storage 

 @author NT
 @version 1.0
 */
//________________________________________________________________________



class SeparationPlant : public CLASSBackEnd
{
public :


//********* Constructor/Destructor Method *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	
	SeparationPlant();				///< Normal Constructor.

	//{
	/// Special Constructor.
	/*!
	 Make a new SeparationPlant
	 \param Log CLASSLogger used for the log...
	 \param separationtime duration of the SeparationPlant
	 */
	SeparationPlant(CLASSLogger* Log, cSecond separationtime = 0. *3600.*24.*365.25); //!<
	//}



	~SeparationPlant();	///< Normal Destructor.
	//@}




//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{


	void SetStorageDestination(CLASSBackEnd*  storagedestination, IsotopicVector isotopicvector, cSecond destinationstartingtime);

	void AddIV(IsotopicVector IV);

	void SetPutToWaste(bool val)		{ fPutToWaste = val; }		//!< Set True if IV goes to waste after cooling false instead


	using CLASSBackEnd::SetName;

	//@}




//********* Get Method *********//

	/*!
	 \name Get Method
	 */
	//@{

	bool		GetPutToWaste()	const	{ return fPutToWaste; }		//!< Return True if IV goes to waste after cooling false instead

	//@}




//********* IsotopicVector Managment Method *********//

	/*!
	 \name IsotopicVector Managment Method
	 */
	//@{

	vector<cSecond>	GetCoolingStartingTime() const
						{ return GetIVArrayArrivalTime(); }	//!< Return the vector of Cooling Sstarting Time
	//@}




//********* Other Method *********//

	//@}
	/*!
	 \name Other Method
	 */
	//@{

	
	//@}

protected :
	
	
	
//********* Internal Parameter *********//
	bool			fPutToWaste;	//!< True if IV goes to waste after cooling false instead
	vector<CLASSBackEnd* > 	fDestinationStorage;	//!< Vector containing destination storage of the IV in the Separation Plant
	vector<IsotopicVector >	fDestinationStorageIV;	//!< Vector containing destination storage of the IV in the Separation Plant
	vector<cSecond>			fDestinationStorageStartingTime; 	//!< Vector containing destination storage starting time of the IV in the Separation Plant

//********* Isotopic Quantity *********//
//--------- Cooling ---------//
	vector<int>		fCoolingIndex;		///< Vector of the Cooling Index
	int			fCoolingLastIndex;	//!< Number of Cooling IV Treated
	vector<int>		fCoolingEndOfCycle;	//!< Index of the Cooling IV reaching the End of a Cooling Cycle


//********* Private Method *********//




	ClassDef(SeparationPlant,3);
};

#endif
