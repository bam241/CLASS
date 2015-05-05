#ifndef _Storage_
#define _Storage_

/*!
 \file 
 \brief Header file for Storage class.
 */


#include <vector>

#include "CLASSBackEnd.hxx"
#include "IsotopicVector.hxx"


using namespace std;
typedef long long int cSecond;

class CLASSLogger;
class DecayDataBank;

//-----------------------------------------------------------------------------//
//! Defines a Storage object

/*!
 A Storage is a CLASSBackEnd facility. It is almost the same as a Pool with a
 infinite cooling time.
 A CLASSFacility can take IsotopicVector(s) contained in a Storage but a Storage 
 cannot send its content in other CLASSFacility (its a kind of passive facility)
 
 @author BaM
 @version 2.0
 */
//________________________________________________________________________



class Storage : public CLASSBackEnd
{
public :


//********* Constructor/Destructor Method *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{

 	Storage();		///< Normal Constructor.

	//{
	/// CLASSLogger Constructor.
	/*!
	 Use to create an empty Storage with a CLASSLogger
	 \param log : used for the log.
	 */
 	Storage(CLASSLogger* log);
	//}


	//{
	/// Special Constructor.
	/*!
	 Make a new Storage
	 \param log : used for the log.
	 \param evolutivedb : DecayDataBank for decay management
	 */
	Storage(CLASSLogger* log, DecayDataBank* evolutivedb);
	//}


 	~Storage(); 	///< Normal Destructor.

	//@}




//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{

	using CLASSBackEnd::SetName;
	using CLASSBackEnd::SetIsStorageType;

	//@}

//********* Storage specific Method *********//

	/*!
	 \name Storage specific methods
	 */
	//@{
		
	void TakeFractionFromStock(int IVId,double fraction);		//!< Take a part from an IV in sotck;
	void TakeFromStock(IsotopicVector isotopicvector);		//!< Take an entire IV from stock


	void AddIV(IsotopicVector isotopicvector);					//!< Add an Isotopicvector to the IVArray
	void AddToStock(IsotopicVector isotopicvector) {AddIV(isotopicvector);}		//!< Add an Isotopicvector to the IVArray
	void RemoveEmptyStocks(); //!< delete the empty Isotopicvector(s) contained in IVArray

	//@}




//********* Evolution Method *********//

	/*!
	 \name Evolution Method
	 */
	//@{
	
	void Evolution(cSecond t);		//!< Perform the evolution until time t

	//@}

	//********* In/Out Method *********//

	/*!
	 \name In/Out Method
	 */
	//@{

	//{
	/// Write the Isotope composition of all IsotopicVector stored.
	/*!
	 Make a new reactor
	 \param filename : CLASSLogger used for the log.
	 \param date : only use to write a date in the file
	 */
	void Write(string filename,cSecond date = -1);
	//}

	//@}

protected :
	
//********* Isotopic Quantity *********//



//********* Private Method *********//
	void	StorageEvolution(cSecond t);					//!< Deal the Storage Decay Evolution




	ClassDef(Storage,3);
};

#endif
