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
/*!
 Define a nuclei as : Z A I.
 The aim of this class is to discribe each CLASSNucleiFiliation connecting nuclei to their daughter nuclei through a nuclear process using the correct branching ratio.
 Note That this class can be used to connect nulei to their daughter through any reaction (Fission, capture...) or natural decay process...
 
 In the case of Decay process,  it is possible to add the decay constant to the branching ratio (ie : log(2)/HalfLife, see IrradiationModel.cxx )....
 
 
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
	 \name Get Method
	 */
	//@{
	map<ZAI, IsotopicVector> GetNucleiFIliation() const {return fNucleiFiliation;}	//!< Return the full filiation list
	vector<ZAI> GetZAIList() const;			//!< Return the list of mother ZAI present in the Filiation list
	int size() const{return (int)fNucleiFiliation.size();}	//!< Return the Number of Mother ZAI (then filiation path)
	
	IsotopicVector GetFiliation(ZAI Mother);	//!< Return the filiation isotopic Vector of the ZAI Mother
	
	ZAI GetArtificialDecay(ZAI Mother);			//!< Make an artificial and instateneus decay on the ZAI, (desexcitation, or Beta decay)
	
	//}
	
	//********* Add Method *********//
	
	/*!
	 \name Adding Method
	 */
	//@{
	void Add(ZAI Mother, IsotopicVector Daughter );		//!< Add A ZAI and its IsotopicVector of Daughter to the Filiation
	//}
	
	
	//********* Modification Method *********//
	
	/*!
	 \name Modification Method
	 */
	//@{
	
	
	void FiliationCleanUp(map<ZAI, int> GoodNuclei, CLASSNucleiFiliation CuttedNuclei);	//!< Cutting all pathway to end each path on a nuclei in the GoodList following the CuttedNuclei, if nuclei is neither in the GoodNuclei list or in CuttedNuclei, then Artificial Decay are performed
	
	void SelfFiliationCleanUp(map<ZAI, int> GoodNuclei);	//!< Cutting All the pathway ending on a nuclei not present as a Motehr Nuclei
	
	void NormalizeBranchingRatio(double Value = 1);		//!< Normalization of all the branching ratio to 1
	
	void NormalizeBranchingRatio(ZAI Mother, double Value);	//!< Normalize the branching ratio pathway of the Mother ZAI to the set value
	
	//}
	
	
	protected :
	
	map<ZAI, IsotopicVector> fNucleiFiliation;		//! Map of all the pathway
	
	
};


#endif
