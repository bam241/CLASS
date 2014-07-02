
#ifndef _CLASSBACKEND_HXX
#define _CLASSBACKEND_HXX


/*!
 \file
 \brief Header file for CLASSFacility class.

 */

#include <string>
#include <fstream>

#include "CLASSFacility.hxx"
#include "IsotopicVector.hxx"
#include "DecayDataBank.hxx"

#include "TNamed.h"

using namespace std;
typedef long long int cSecond;

//-----------------------------------------------------------------------------//
/*!
 Define a CLASS Facility.
 The aim of these class is synthetyse all the commum properties of the nuclear facilities which are involve in the BackEnd Fuel cycle.


 @author BaM
 @version 2.0
 */
//________________________________________________________________________




class CLASSBackEnd : public CLASSFacility
{
	public :
	///< Normal Constructor.
	CLASSBackEnd();
	CLASSBackEnd(LogFile* log);
	CLASSBackEnd(LogFile* log, cSecond cycletime);

	//********* Get Method *********//
	/*!
	 \name Get Function
	 */
	//@{

	DecayDataBank*	GetDecayDataBank()		{ return fDecayDataBase;}	//!< Return the pointer to the decay DataBank

	vector<IsotopicVector> GetIVArray()	const	{ return fIVArray; }		//!< Return the IsotopicVector Array
	vector<cSecond>	GetIVArrayArrivalTime()	const	{ return fIVArrayArrivalTime;}	//!<Return the pointer to the OUtBackEndFacility

	int		GetIVNumber()		const	{ return fIVArray.size();}
	bool		GetStorageType()	const	{ return fIsStorageType;}	//!< Return the storageType
	CLASSBackEnd*	GetOutBackEndFacility()	const	{ return fOutBackEndFacility;}	//!<Return the pointer to the OUtBackEndFacility
	IsotopicVector  GetIV(int i)		const	{ if(i < (int)fIVArray.size()) return fIVArray[i];
								else return IsotopicVector(); }


	//@}

	//********* Set Method *********//
	/*!
	 \name Set Function
	 */
	//@{
	void		SetIsStorageType(bool val = true)		{ fIsStorageType = val;}	//! Set the fIsStorage bool
	void		SetDecayDataBank(DecayDataBank* decayDB)	{ fDecayDataBase = decayDB;}	//! Set the Decay DataBank
	virtual	void	SetIVArray(vector<IsotopicVector> ivarray)	{ fIVArray = ivarray; }		//!< Set The isotopicVector Array
	virtual void	SetOutBackEndFacility(CLASSBackEnd* befacility)	{ fOutBackEndFacility = befacility;
									  fIsStorageType = true; } //! Set a Out Facility for the fuel
	//@}


	/*!
	 \name BackEndFacility specific Method
	 */
	//@{
	virtual void	AddIV(IsotopicVector isotopicvector);	//!< Add an Isotopicvector to the IVArray
	void		ClearIVArray();					//!< Empty the IVArray removing all fuel stored

	//@}
	virtual void Evolution(cSecond t)	{}	//!< Performe the Evolution to the Time t


	protected :
	IsotopicVector		GetDecay(IsotopicVector isotopicvector, cSecond t);	//!< Get IsotopicVector Decay at the t time
	vector<IsotopicVector>	fIVArray;					///< Vector containning all the fuel stored.
	vector<cSecond>		fIVArrayArrivalTime;					///< Vector containning all the fuel stored.
	CLASSBackEnd*		fOutBackEndFacility;					//!< Facility getting the fuel at the end of the cycle

	//********* Internal Parameter *********//
	private :
	bool		fIsStorageType;		//!< True if there is not OutBAckEndFacility (like a storage...)
	DecayDataBank*	fDecayDataBase;		//!< Pointer to the Decay DataBase

	ClassDef(CLASSBackEnd,2);
};

#endif

