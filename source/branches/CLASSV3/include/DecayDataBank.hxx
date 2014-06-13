#ifndef __DecayDataBank_HXX__
#define __DecayDataBank_HXX__

/*!
 \file
 \brief Header file for DecayDataBank class.
 @version 2.0
 */

#include "CLASSObject.hxx"
#include "TMatrix.h"
#include "EvolutionData.hxx"
#include "IsotopicVector.hxx"
#include "DynamicalSystem.hxx"

#include <map>
#include <vector>


using namespace std;
typedef long long int cSecond;

class ZAI;
class LogFile;

double ReactionRateWeightedDistance(IsotopicVector IV1, EvolutionData DB );
double ReactionRateWeightedDistance(EvolutionData DB, IsotopicVector IV1  );

//-----------------------------------------------------------------------------//
/*!
 Define a DecayDataBank.
 The aim of these class is describe the evolution of "all" evoluting system in CLASS.
 
 For the Decay Matrix the DecayDataBank iq mainly contain a map of <ZAI,EvolutionData>.This map do the correspondance between a ZAI and its decay evolution (containing all the daughter nuclei comming from the decay a the original ZAI).
 
 @author BaM
 @author Marc
 @author PTO for a part the Decay management -- steal from MURE (Even if he does not kown it!! :))
 @version 2.0
 */
//________________________________________________________________________



class DecayDataBank : public CLASSObject
{
	
	public :
	
	
	//********* Constructor/Desctructor *********//
	
	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	/// Normal Constructor.
	DecayDataBank();
	
	//{
	/// Special Constructor.
	/*!
	 Use to load a LogFile
	 \param LogFile LogFile used for the log...
	 \param DB_index_file path to the index file
	 \param setlog if the log are stored in the LogFile
	 \param olfreadmethod true if the old format of EvolutionData are used (ie without the key word such as Inv, XSFiss...)
	 */
	DecayDataBank(LogFile* Log, string DB_index_file, bool setlog = true, bool olfreadmethod = true );
	//}
	
	//{
	/// Normal Destructor.
	/*!
	 Delete de DecayDataBank and all associated EvolutionData...
	 */
	~DecayDataBank();
	//}
	
	//{
	/// Reset the DecayDataBank.
	/*!
	 Use to reset the DecayDataBank to its default values  whihout deleting the EvolutionData (which contain pointer... ).
	 it does just clear the different maps
	 */
	void Clear();
	//}
	//@}
	
	
	
	
	//********* Get Method *********//
	/*!
	 \name Get Method
	 */
	//@{
	map<ZAI ,EvolutionData >	GetDecayDataBank()	const	{ return fDecayDataBank; }	//!< Return the DecayDataBank
	string 				GetDataBaseIndex()	const	{ return fDataBaseIndex; }	//!< Return the index Name
	bool 				IsDefine(const ZAI& zai)	const;					//!< True the key is define, false unstead
	
	string	GetDataFileName()	const { return fDataFileName; }
	string	GetDataDirectoryName()  const { return fDataDirectoryName; }
	
	double  GetShorstestHalflife()	const { return fShorstestHalflife; }
	
	//@}
	
	
	
	
	//********* Set Method *********//
	
	/*!
	 \name Set Method
	 */
	//@{
	
	void SetDecayDataBank(map<ZAI ,EvolutionData > mymap)	{ fDecayDataBank = mymap; }	//!< Set the DecayDataBank map
	
	void SetDataBaseIndex(string database) { fDataBaseIndex = database;; ReadDataBase(); }	//!< Set the Name of the database index
	
	void SetOldReadMethod(bool val)			{ fOldReadMethod = val; ReadDataBase();}			///< use the old reading method
	
	//}
	
	
	
	
	//********* Evolution Method *********//
	
	//@}
	/*!
	 \name Evolution Method
	 */
	//@{
	
	
	IsotopicVector	Evolution(const ZAI& zai, double dt);	///< Return the Product IsotopicVector evolution from zai during a dt time
	
	//@}
	
	
	
	
	//********* Other Method *********//
	/*!
	 \name Other Method
	 */
	//@{
	void	ReadDataBase();				///< read the index file and fill the evolutionData map
	
	void Print() const;
	
	//@}
	
	
	
	
	
	protected :
	
	double  fShorstestHalflife;
	int	fZAIThreshold;	//!< Highest Mass deal bye the evolution (default 90)
	
	string			fDataFileName;		///< Name of the decay list
	string			fDataDirectoryName;	///< Path to the decay list file
	
	map<ZAI, EvolutionData>	fDecayDataBank;		///< DataBanck map
	
 	string			fDataBaseIndex;			///< Name of the index
	
	bool			fOldReadMethod;		///< use old DB format
	
};



#endif
