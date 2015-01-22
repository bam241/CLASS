
#ifndef _CLASSBACKEND_HXX
#define _CLASSBACKEND_HXX


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
	CLASSBackEnd(int type = 0);
	CLASSBackEnd(CLASSLogger* log,int type = 0);
	CLASSBackEnd(CLASSLogger* log, cSecond cycletime, int type = 0);

	//********* Get Method *********//
	/*!
	 \name Get Function
	 */
	//@{


	vector<IsotopicVector> GetIVArray()	const	{ return fIVArray; }		//!< Return the IsotopicVector Array
	vector<cSecond>	GetIVArrayArrivalTime()	const	{ return fIVArrayArrivalTime;}	//!<Return the pointer to the OUtBackEndFacility

	int		GetIVNumber()		const	{ return fIVArray.size();}
	bool		GetStorageType()	const	{ return fIsStorageType;}	//!< Return the storageType
	IsotopicVector  GetIV(int i)		const	{ if(i < (int)fIVArray.size()) return fIVArray[i];
								else return IsotopicVector(); }
#ifndef __CINT__
	DecayDataBank*	GetDecayDataBank()		{ return fDecayDataBase;}	//!< Return the pointer to the decay DataBank
	CLASSBackEnd*	GetOutBackEndFacility()	const	{ return fOutBackEndFacility;}	//!<Return the pointer to the OUtBackEndFacility
	virtual map<cSecond,int> GetTheBackEndTimePath();

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
										fIsStorageType = false; } //! Set a Out Facility for the fuel

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
	virtual void Evolution(cSecond t)	{}	//!< Performe the Evolution to the Time t
	void UpdateInsideIV();


	protected :
	IsotopicVector		GetDecay(IsotopicVector isotopicvector, cSecond t);	//!< Get IsotopicVector Decay at the t time
	vector<IsotopicVector>	fIVArray;					///< Vector containning all the fuel stored.
	vector<cSecond>		fIVArrayArrivalTime;					///< Vector containning all the fuel stored.

#ifndef __CINT__
	CLASSBackEnd*		fOutBackEndFacility;					//!< Facility getting the fuel at the end of the cycle
#endif

	//********* Internal Parameter *********//
	private :
	bool		fIsStorageType;		//!< True if there is not OutBAckEndFacility (like a storage...)

#ifndef __CINT__
	DecayDataBank*	fDecayDataBase;		//!< Pointer to the Decay DataBase
#endif

	ClassDef(CLASSBackEnd,2);
};

#endif

