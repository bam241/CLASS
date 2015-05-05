#ifndef _CLASSNucleiFiliation_
#define _CLASSNucleiFiliation_

/*!
 \file
 \brief Header file for CLASSNucleiFiliation classes.
 */

#include "CLASSObject.hxx"
#include "IsotopicVector.hxx"

using namespace std;

//-----------------------------------------------------------------------------//
//! Handles connection between nuclei (decay/reaction)

/*!
 Define a nuclei as : Z A I.
 The aim of this class is to discribe each CLASSNucleiFiliation. It connects nuclei to their daughter nuclei through a nuclear process using the correct branching ratio.
 Note that this class can be used to connect nulei to their daughter through any reaction (Fission, capture...) or natural decay process...
 
 In the case of decay process,  it is possible to put (by multiplying) the decay constant to the branching ratio (ie : BR*log(2)/HalfLife, see IrradiationModel)....
 
 
 @author BaM
 @version 2.0
 */
//________________________________________________________________________



class CLASSNucleiFiliation : public CLASSObject
{
public:
	
	
	//********* Constructor/Destructor Method *********//
	
	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	
	CLASSNucleiFiliation();	///< Default constructor
	
	
	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	//{
	/*!
	 Use to load a CLASSNucleiFiliation
	 \param log : used for the log.
	 */
	CLASSNucleiFiliation(CLASSLogger* log);	///< Default constructor
	//}
	
	CLASSNucleiFiliation(const CLASSNucleiFiliation& CNF); ///< Copy Constructor
	
	
	~CLASSNucleiFiliation();	///< Normal Destructor.
	
	//@}
	
	
	//********* Get Method *********//
	
	/*!
	 \name Set/Get Method
	 */
	//@{
	map<ZAI, IsotopicVector> GetNucleiFIliation() const {return fNucleiFiliation;}	//!< Return the full filiation list
	vector<ZAI> GetZAIList() const;			//!< Return the list of mother ZAI present in the filiation list
	int size() const{return (int)fNucleiFiliation.size();}	//!< Return the number of mother ZAI (then filiation path)
	
	IsotopicVector GetFiliation(ZAI Mother);	//!< Return the filiation isotopic vector of the ZAI mother
	
	ZAI GetArtificialDecay(ZAI Mother);			//!< Make an artificial and instantaneus decay of the ZAI, (desexcitation, or Beta decay)
	
	
	void SetNucleiFIliation(map<ZAI, IsotopicVector> Fiiliation) { fNucleiFiliation = Fiiliation;}	//!< Return the full filiation list

	//}
	
	//********* Add Method *********//
	
	/*!
	 \name Adding Method
	 */
	//@{
	void Add(ZAI Mother, IsotopicVector Daughter );		//!< Add A ZAI and its IsotopicVector of daughter(s) to the filiation
	//}
	
	
	//********* Modification Method *********//
	
	/*!
	 \name Modification Method
	 */
	//@{
	
	
	void FiliationCleanUp(map<ZAI, int> GoodNuclei, CLASSNucleiFiliation CuttedNuclei);	//!< Cutting all pathway until each path ends on a nuclei in the GoodList following the CuttedNuclei. If nuclei are neither in the GoodNuclei list or in CuttedNuclei, then artificial decay are performed
	
	void SelfFiliationCleanUp(map<ZAI, int> GoodNuclei);	//!< Cutting all the pathway ending on a nuclei not present as a motehr nuclei
	
	void NormalizeBranchingRatio(double Value = 1);		//!< Normalization of all the branching ratio to 1
	
	void NormalizeBranchingRatio(ZAI Mother, double Value);	//!< Normalize the branching ratio pathway of the Mother ZAI to the set value
	
	//}
	
	
	protected :
	
	map<ZAI, IsotopicVector> fNucleiFiliation;		//! Map of all the pathway
	
	
};


#endif
