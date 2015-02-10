#ifndef __SeparationPlant_HXX__
#define __SeparationPlant_HXX__
/*!
 \file
 \brief Header file for SeparationPlant class.
 */

#include <string>
#include <map>

#include "CLASSConstante.hxx"
#include "CLASSBackEnd.hxx"
#include "Storage.hxx"
#include "IsotopicVector.hxx"

using namespace std;
typedef long long int cSecond;

class CLASSBackEnd;
class CLASSLogger;
class DecayDataBank;

//-----------------------------------------------------------------------------//
//!  Defines a SeparationPlant.

/*!
 The aim of this class is to separate an IV into several IV (MA, Pu, PF, etc...) and
 to send it to one or several Storage

 @author NT
 @author BaM
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
	 \param log : used for the log.
	 */
	SeparationPlant(CLASSLogger* Log);
	//}

	~SeparationPlant();	///< Normal Destructor.
	//@}


//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{


	void SetBackEndDestination(CLASSBackEnd*  storagedestination, IsotopicVector isotopicvector, cSecond destinationstartingtime); //!< Tell Separation plant to begin separation  at time destinationstartingtime according efficiency (between [0-1]) defined in isotopicvector and to send the separated materials to storagedestination

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


	map<cSecond,int> GetTheBackEndTimePath();


//********* IsotopicVector Managment Method *********//

	/*!
	 \name IsotopicVector Managment Method
	 */
	//@{

	vector<cSecond>	GetCoolingStartingTime() const
						{ return GetIVArrayArrivalTime(); }	//!< Return the vector of Cooling starting Time
	//@}

protected :
	
	
	
//********* Internal Parameter *********//
	bool			fPutToWaste;	//!< True if IV goes to waste after cooling false instead
	vector<CLASSBackEnd* > 	fDestinationStorage;	//!< Vector containing destination storage of the IV in the Separation Plant
	vector<IsotopicVector >	fDestinationStorageIV;	//!< Vector containing destination storage of the IV in the Separation Plant
	vector<cSecond>		fDestinationStorageStartingTime; 	//!< Vector containing destination storage starting time of the IV in the Separation Plant

//********* Private Method *********//




	ClassDef(SeparationPlant,3);
};

#endif
