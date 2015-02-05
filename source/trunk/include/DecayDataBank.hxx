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
class CLASSLogger;

double ReactionRateWeightedDistance(IsotopicVector IV1, EvolutionData DB );
double ReactionRateWeightedDistance(EvolutionData DB, IsotopicVector IV1  );

//-----------------------------------------------------------------------------//
//! Describes outcore radioactive decays 

/*!
 Define a DecayDataBank.
 The aim of these class is to describe the evolution of "all" evoluting systems in CLASS.
 
 For the Decay Matrix the DecayDataBank  mainly contains a map of <ZAI,EvolutionData>.This map do the correspondance between a ZAI and its decay evolution (containing all the daughter nuclei comming from the decay of the original ZAI and quantities evolutions).
 
 @author BaM
 @author Marc
 @author PTO for a part of the Decay management -- steal from MURE (Even if he does not kown it!! :))
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
	 Use to load a DecayDataBank
	 \param DB_index_file : path to the index file
	 \param olfreadmethod : true if the old format of EvolutionData are used (deprecated) (ie without the key word such as Inv, XSFiss...)
	 */
	DecayDataBank(string DB_index_file, bool olfreadmethod = false );
	//}
	//{
	/// Special Constructor.
	/*!
	 Use to load a DecayDataBank
	 \param Log : CLASSLogger used for the log.
	 \param DB_index_file : path to the index file
	 \param olfreadmethod : true if the old format of EvolutionData are used (ie without the key word such as Inv, XSFiss...)
	 */
	DecayDataBank(CLASSLogger* Log, string DB_index_file, bool olfreadmethod = false );
	//}
	
	//{
	/// Normal Destructor.
	/*!
	 Delete the DecayDataBank and all associated EvolutionData(s)...
	 */
	~DecayDataBank();
	//}
	
	//{
	/// Reset the DecayDataBank.
	/*!
	 Use to reset the DecayDataBank to its default values whihout deleting the EvolutionData (which contain pointer... ).
	 it just clears the different maps
	 */
	void Clear();
	//}
	//@}
	
	
	
	
	//********* Get Method *********//
	/*!
	 \name Get Method
	 */
	//@{
	map<ZAI ,EvolutionData > GetDecayDataBank()	const	{ return fDecayDataBank; }	//!< Return the DecayDataBank
	bool 			IsDefine(const ZAI& zai)	const;				//!< True if the key is define, false unstead

	string 			GetDataBaseIndex()	const	{ return fDataBaseIndex; }	//!< Return the index name

	IsotopicVector		GetDecay(IsotopicVector isotopicvector, cSecond t);	//!< Get IsotopicVector decay at time t

	//@}
	
	
	
	
	//********* Set Method *********//
	
	/*!
	 \name Set Method
	 */
	//@{
	
	void SetDecayDataBank(map<ZAI ,EvolutionData > mymap)
						{ fDecayDataBank = mymap; }	//!< Set the DecayDataBank map
	
	void SetDataBaseIndex(string database)	{ fDataBaseIndex = database;; ReadDataBase(); }	//!< Set the name of the database index
	
	void SetOldReadMethod(bool val)		{ fOldReadMethod = val; ReadDataBase();}			///< use the old reading method
	
	//}
	
	
	
	
	//********* Evolution Method *********//
	
	//@}
	/*!
	 \name Evolution Method
	 */
	//@{
	
	IsotopicVector	Evolution(const ZAI& zai, double dt);	///< Return the IsotopicVector from the decay of zai during a dt period
	
	//@}
	
	
	
	
	//********* Other Method *********//
	/*!
	 \name Other Method
	 */
	//@{
	void	ReadDataBase();				///< Read the index file and fill the EvolutionData map
	
	void Print() const;
	
	//@}
	
	
	
	
	
	protected :
	
	map<ZAI, EvolutionData>	fDecayDataBank;		///< DataBank map
 	string			fDataBaseIndex;		///< Name of the index
	bool			fOldReadMethod;		///< use old DB format
	
};



#endif
