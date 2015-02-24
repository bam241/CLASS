
#ifndef _CLASSBACKEND_
#define _CLASSBACKEND_


/*!
 \file
 \brief Header file for CLASSFacility class.
 
 */

#include <string>
#include <fstream>
#include <map>


#include "CLASSFacility.hxx"
#include "IsotopicVector.hxx"
#include "DecayDataBank.hxx"

#include "TNamed.h"

using namespace std;

//-----------------------------------------------------------------------------//
//! Class defining the common properties of all back end fuel cycle facilities

/*!
 Define a CLASS Facility.
 The aim of these class is to gather all the commom properties of the
 facilities which are involve in the BackEnd Fuel cycle.
 
 @author BaM
 @version 2.0
 */
//________________________________________________________________________


class CLASSBackEnd : public CLASSFacility
{
	public :
	//{
	/// Default Constructor.
	/*!Create an empty CLASSBackEnd
	 \param type
	 \li -2 :SeparationPlant
	 \li -1 :Storage
	 \li 8 :Pool
	 */
	CLASSBackEnd(int type = 0);
	//}
	
	//{
	/// CLASSLogger Constructor.
	/*!
	 Create an empty CLASSBackEnd loading a CLASSLogger
	 \param log : used for the log.
	 \param type
	 \li -2 :SeparationPlant
	 \li -1 :Storage
	 \li 8 :Pool     */
	CLASSBackEnd(CLASSLogger* log,int type = 0);
	//}
	
	//{
	/// Cycle time Constructor.
	/*!
	 Create an empty CLASSBackEnd loading a CLASSLogger
	 \param log : used for the log.
	 \param cycletime Cycle time of the CLASSBackend (e.g. Cooling time for the pool) in [s],
	 \param type -2 :SeparationPlant -1 : Storage ; 8 :Pool
	 */
	CLASSBackEnd(CLASSLogger* log, cSecond cycletime, int type = 0);
	//}
	//@}

	//********* Get Method *********//
	/*!
	 \name Get Function
	 */
	//@{
	
	
	vector<IsotopicVector> GetIVArray()	const	{ return fIVArray; }		//!< Return the IsotopicVector Array
	vector<cSecond>	GetIVArrayArrivalTime()	const	{ return fIVArrayArrivalTime;}	//!<Vector of arrival time of each IV in the CLASSBackEnd
	
	int		GetIVNumber()		const	{ return fIVArray.size();} //!< Return the number of Isotopic Vector present in the CLASSBackEnd object
	bool		GetStorageType()	const	{ return fIsStorageType;}	//!< Return the storageType : True if it is a Storage
	IsotopicVector  GetIV(int i)		const	{ if(i < (int)fIVArray.size()) return fIVArray[i];
									else return IsotopicVector(); }
#ifndef __CINT__
	DecayDataBank*	GetDecayDataBank()		{ return fDecayDataBase;}	//!< Return the pointer to the decay DataBank
	CLASSBackEnd*	GetOutBackEndFacility()	const	{ return fOutBackEndFacility;}	//!<Return the pointer to the OUtBackEndFacility
	virtual map<cSecond,int> GetTheBackEndTimePath();	//!< Get the full path 
	
#endif
	
	//@}
	
	//********* Set Method *********//
	/*!
	 \name Set Function
	 */
	//@{
	void		SetIsStorageType(bool val = true)		{ fIsStorageType = val;}	//! Set the fIsStorage bool
	virtual	void	SetIVArray(vector<IsotopicVector> ivarray)	{ fIVArray = ivarray; }		//!< Set The isotopicVector Array
#ifndef __CINT__
	void		SetDecayDataBank(DecayDataBank* decayDB)	{ fDecayDataBase = decayDB;}	//! Set the Decay DataBank
	virtual void	SetOutBackEndFacility(CLASSBackEnd* befacility)	{ fOutBackEndFacility = befacility;
		fIsStorageType = false; } //! Set an out Facility
	
#endif
	
	using CLASSFacility::SetName;
	
	//@}
	
	
	/*!
	 \name BackEndFacility specific Method
	 */
	//@{
	virtual void	AddIV(IsotopicVector isotopicvector);	//!< Add an Isotopicvector to the IVArray
	void		ClearIVArray();					//!< Empty the IVArray removing all fuel stored
	
	//@}
	virtual void Evolution(cSecond t)	{}	//!< Performs the Evolution until time t
	void UpdateInsideIV();
	
	
	protected :
	IsotopicVector		GetDecay(IsotopicVector isotopicvector, cSecond t);	//!< Get IsotopicVector decay at time t [s]
	vector<IsotopicVector>	fIVArray;						///< Vector containning all the fuel stored.
	vector<cSecond>		fIVArrayArrivalTime;					///< Vector containning the arrival time of each fuel in [s]
	
#ifndef __CINT__
	CLASSBackEnd*		fOutBackEndFacility;					//!< Facility getting the fuel at the end of the cycle
#endif
	
	//********* Internal Parameter *********//
	private :
	bool		fIsStorageType;		//!< True if there is no out CLASSBackEnd facility (like a Storage...)
	
#ifndef __CINT__
	DecayDataBank*	fDecayDataBase;		//!< Pointer to the DecayDataBank
#endif
	
	ClassDef(CLASSBackEnd,2);
};

#endif

