#ifndef __Storage_HXX__
#define __Storage_HXX__

/*!
 \file 
 \brief Header file for Storage class.
 */


#include <vector>

#include "CLSSFacility.hxx"
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



class Storage : public  CLSSFacility
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
	
	void SetDecayDataBase(DecayDataBank* ddb)	{ fDecayDataBase = ddb; }	//!< Set the pointer to the Decay DataBase
	void SetStock(vector<IsotopicVector> IVsstock)	{ fIVStock = IVsstock; }		//!< Set The Storage isotopicVector

	//@}




//********* Get Method *********//

	/*!
	 \name Get Method
	 */
	//@{
	
	DecayDataBank* GeDecayDataBase()	const		{ return fDecayDataBase; }	//!< Return the pointer to the Decay DataBase
	vector<IsotopicVector> GetStock()	const		{ return fIVStock; }		//!< Return the Storage IsotopicVector
	IsotopicVector GetFullStock()		const		{ return GetInsideIV(); }	//!< Return the Full Storage

	//@}




//********* Storage specific Method *********//

	/*!
	 \name Storage specific Method
	 */
	//@{
	
	void ClearStock();										//!< Empty the stock removing all fuel stored
	
	void AddToStock(ZAI zai, double quantity)		{ AddToStock(zai*quantity); }		//!< Add a ZAI*quantity to the Storage
	void AddToStock(IsotopicVector isotopicvector);							//!< Add an Isotopicvector to the Storage
	void AddToFullStock(ZAI zai, double quantity)		{ fInsideIV += zai*quantity; }		//!< Add a ZAI*quantity to the Storage
	void AddToFullStock(IsotopicVector isotopicvector)	{ fInsideIV += isotopicvector; }	//!< Add a IsotopicVector to the Storage

	void TakeFractionFromStock(int IVId,double fraction);						//!< Take a part from an IV in sotck;
	void TakeFromStock(IsotopicVector isotopicvector);						//!<

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
	
//********* Internal Parameter *********//

	DecayDataBank*	fDecayDataBase;	//!< Pointer to the Decay DataBase


//********* Isotopic Quantity *********//

//---------- Storage ----------//
	vector<IsotopicVector>	fIVStock;		///< Vector containning all the fuel stored.



//********* Private Method *********//
	IsotopicVector GetDecay(IsotopicVector isotopicvector, cSecond t);	//!< Get IsotopicVector Decay at the t time
	void	StorageEvolution(cSecond t);					//!< Deal the Storage Decay Evolution




	ClassDef(Storage,2);
};

#endif
