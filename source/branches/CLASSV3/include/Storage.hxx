#ifndef __Storage_HXX__
#define __Storage_HXX__

/*!
 \file 
 \brief Header file for Storage class.
 */


#include <vector>

#include "CLASSBackEnd.hxx"
#include "IsotopicVector.hxx"


using namespace std;
typedef long long int cSecond;

class CLASS;
class LogFile;
class DecayDataBank;

//-----------------------------------------------------------------------------//
/*!
 Define a Storage.
 The aim of this class is to deal the store used fuel after the cooling dealing the evolution of all radiaoactive nuclei.

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
	/// LogFile Constructor.
	/*!
	 Use create an empty Stotarage loading a LogFile
	 \param LogFile LogFile used for the log...
	 */
 	Storage(LogFile* log);
	//}


	//{
	/// Special Constructor.
	/*!
	 Make a new reactor
	 \param LogFile LogFile used for the log...
	 \param evolutivedb DataBank for decay management
	 */
	Storage(LogFile* log, DecayDataBank* evolutivedb);
	//}


 	~Storage(); 	///< Normal Destructor.

	//@}




//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{
	//@}


//********* Get Method *********//

	/*!
	 \name Get Method
	 */
	//@{
	//@}




//********* Storage specific Method *********//

	/*!
	 \name Storage specific Method
	 */
	//@{
		
	void AddToFullStock(IsotopicVector isotopicvector)	{ fInsideIV += isotopicvector; }	//!< Add a IsotopicVector to the Storage

	void TakeFractionFromStock(int IVId,double fraction);						//!< Take a part from an IV in sotck;
	void TakeFromStock(IsotopicVector isotopicvector);						//!<


	void AddIV(IsotopicVector isotopicvector);			//!< Add an Isotopicvector to the IVArray

	//@}




//********* Evolution Method *********//

	/*!
	 \name Evolution Method
	 */
	//@{
	
	void Evolution(cSecond t);		//!< Performe the evolution until the Time t

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
	 \param filenam LogFile used for the log...
	 \param data only use to srite a date in the file, theyr is not treatment of the date in this method....
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
