#ifndef __Pool_HXX__
#define __Pool_HXX__
/*!
 \file
 \brief Header file for Pool class.
 */

#include <string>
#include <map>

#include "CLASSBackEnd.hxx"
#include "IsotopicVector.hxx"

using namespace std;
typedef long long int cSecond;

class CLASSBackEnd;
class CLASS;
class LogFile;
class DecayDataBank;

//-----------------------------------------------------------------------------//
/*!
 Define a Pool.
 The aim of the Class is to manage evolution of all out reactor fuel. from Cooling to Waste or storage


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
	/// LogFile Constructor.
	/*!
	 Use create an empty Pool loading a LogFile
	 \param LogFile LogFile used for the log...
	 */
 	Pool(LogFile* log);
	//}


	//{
	/// Special Constructor.
	/*!
	 Make a new EvolutionData
	 \param Log LogFile used for the log...
	 \param abstime time to start the Pool
	 \param coolingtime duration of the cooling.
	 */
	Pool(LogFile* Log, double abstime,
			 double coolingtime = 5*3600.*24.*365.25); //!<
	//}


	//{
	/// Special Special Constructor.
	/*!
	 Make a new EvolutionData
	 \param Log LogFile used for the log...
	 \param Storage storage which get the fuel after the cooling
	 \param abstime time to start the Pool
	 \param coolingtime duration of the cooling.
	 */
	Pool(LogFile* log, CLASSBackEnd* Storage,
			 double abstime = 0,
			 double coolingtime = 5*3600.*24.*365.25); //!<
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
						   SetIsStorageType();
						   fPutToWaste = true; }		//!< Set the Pointer to the Storage

	void SetPutToWaste(bool val)		{ fPutToWaste = val; }		//!< Set True if IV goes to waste after cooling false instead

	void SetIVArray(vector<IsotopicVector> ivarray);			//! not use there (Does nothing!!!)
	void SetIVArray(vector<IsotopicVector> ivarray, vector<cSecond> timearray); //!< Set The isotopicVector Array at the corresponding time

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
						{ return fCoolingStartingTime; }	//!< Return the vector of Cooling Sstarting Time
	void	RemoveIVCooling(int i);					//!< Remove a Cooling IsotopicVector

	void	AddIV(IsotopicVector isotopicvector);			//!< Add an Isotopicvector to the IVArray
	//@}




//********* Other Method *********//

	//@}
	/*!
	 \name Other Method
	 */
	//@{

	void Evolution(cSecond t);		//!< Performe the evolution until the Time t
	void Dump();				//!< Write Modification (exchange between Cooling, Separation and Storage)
	
	//@}

protected :
	
	
	
//********* Internal Parameter *********//
	bool			fPutToWaste;	//!< True if IV goes to waste after cooling false instead


//********* Isotopic Quantity *********//
//--------- Cooling ---------//
	vector<cSecond>		fCoolingStartingTime;	///< Vector of the Cooling Starting Time
	vector<int>		fCoolingIndex;		///< Vector of the Cooling Index
	int			fCoolingLastIndex;	//!< Number of Cooling IV Treated
	vector<int>		fCoolingEndOfCycle;	//!< Index of the Cooling IV reaching the End of a Cooling Cycle


//********* Private Method *********//
	void	CoolingEvolution(cSecond t);					//!< Deal the cooling and then send it to Separation




	ClassDef(Pool,3);
};

#endif
