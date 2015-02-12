#ifndef _IRRADIATIONMODEL_
#define _IRRADIATIONMODEL_


/*!
 \file
 \brief Header file for IrradiationModel class.
 
 
 @author BaM
 @author BLG
 @version 2.0
 */


#include "CLASSObject.hxx"

#include "IsotopicVector.hxx"
#include "CLASSNucleiFiliation.hxx"
#include "EvolutionData.hxx"

#include "TMatrix.h"


#include <map>
#include <vector>


using namespace std;
typedef long long int cSecond;

class ZAI;
class CLASSLogger;
//-----------------------------------------------------------------------------//
//! The Bateman equation solver

/*!
 Define an IrradiationModel.
 An IrradiationModel is a Bateman equation solving method.
 This is the mother class.
 see derivated classes :
 \li @see IM_Matrix
 \li @see IM_RK4
 
 The aim of these class is to gather all the commom properties of all the
 derivated Irradiation Model.
 
 @author BaM
 @author BLG
 @version 3.0
 */
//________________________________________________________________________

class IrradiationModel : public CLASSObject
{
	
	public :
	/*!
	 \name Constructors
	 */
	//@{
	IrradiationModel(); //!< Default constructor

	IrradiationModel(CLASSLogger* log); //!< Logger constructor
	
	//@}
	
	//{
	/// virtual method called to perform the irradiation calculation using a set of cross section.
	/*!
	 Perform the Irradiation Calcultion using the XSSet data
	 \param IV : isotopic vector to irradiate
	 \param XSSet : set of mean cross section to use in order to perform the depletion calculation
	 \param Power : constant power to use for irradation [W]
	 \param irradiationtime : irradiation time [s]
	 */
	virtual	 EvolutionData GenerateEvolutionData(IsotopicVector IV, EvolutionData XSSet, double Power, double cycletime) = 0 ;
	//}

	
	
	
	//********* Get Method *********//
	/*!
	 \name Get Method
	 */
	//@{
	string	GetDataFileName()	const { return fDataFileName; }		// File name where decay constants and branching ratios are written.
	string	GetDataDirectoryName()  const { return fDataDirectoryName; }	//!< Path to fDataFileName

	double  GetShorstestHalflife()	const { return fShorstestHalflife; }	//!< Nuclei with HL below fShorstestHalflife are cut (replaced by their daughter(s))
	

	void GetNuclearProcessMatrix(TMatrixT<double> &myMatrix, ZAI Mother, IsotopicVector ProductedIV, double XSValue = 1);
	
	void BuildReactionFiliation();

	string GetSpectrumType(){return fSpectrumType;}				//!< Get the type of neutron spectrum (thermal or fast)

	//@}
	
	
	
	
	//********* Set Method *********//
	
	/*!
	 \name Set Method
	 */
	//@{
	

	//{
	/// set Fission Energy using a file
	/*!
	 // This method fill the Fission Energy [eV] map using a file
	 // \param FissionEnergyFile: filename containing the Fission Energy of some nuclei 
	 (format : Z A I Energy[eV])
	 */
	void SetFissionEnergy(string FissionEnergyFile);
	//}
	
	//{
	/// set Fission Energy for a ZAI using ZAI(Z,A,I)
	/*!
	 // This method fill the Fission Energy map of a set ZAI
	 // \param zai : the ZAI
	 // \param E : Energy released by fission for nuclei zai [eV]
	 */
	void SetFissionEnergy(ZAI zai, double E);
	//}
	
	//{
	/// set Fission Energy for a ZAI using the Z, A, I
	/*!
	 // This method fill the Fission Energy map of a set ZAI
	 // \param Z : Z of the ZAI
	 // \param A : A of the ZAI
	 // \param I : I of the ZAI
	 // \param E : Fission energy of the ZAI [eV]
	 */
	void SetFissionEnergy(int Z, int A, int I, double E )   { SetFissionEnergy(ZAI(Z,A,I), E);}
	//}
	
	void SetShortestHalfLife(double halflife)	{ fShorstestHalflife = halflife;}	//!< Set the Half Life cut
	void LoadFPYield(string SponfaneusYield, string ReactionYield);				//!< Build Fision Yields maps
	
	void SetSpectrumType(string type);					//!< Set the type of neutron spectrum (thermal or fast)
	
	//@}
	
	
	
	//********* Evolution Method *********//
	/*!
	 \name Evolution Method
	 */
	//@{
	

	void	BuildDecayMatrix();			//!< Build the Decay Matrix for the futur time step
	void    LoadDecay();				//!< Load the decay properties (HL,BR)

	void	NuclearDataInitialization();		//!< Build Decay matrices & read FpYields if any
	//@}
	

	//********* Other Method *********//
	/*!
	 \name Other Method
	 */
	//@{
	void Print() const;
	
	int  GetZAIThreshold(){return fZAIThreshold;} //!< Gives the threshold (in charge number Z). The nuclei below this threshold are not managed
	//@}
	
	
	
	protected :
	
	double  fShorstestHalflife;	//!< Limit on the half life of nuclei to take it into account
	int	fZAIThreshold;		//!< Lowest Mass deal by the evolution (default 90)
	string	fDataFileName;		//!< Name of the decay list
	string	fDataDirectoryName;	//!< Path to the decay list file

	map<ZAI, double >	fFissionEnergy;	//!< Store the Energy per fission use for the flux normalisation.
	
	map<ZAI, int> fMatrixIndex;		//!< correspondance matrix from ZAI to the column (or line) of the different Reaction/Decay matrix
	vector<ZAI> fReverseMatrixIndex;	//!< correspondance matrix from the column (or line) of the different Reaction/Decay matrix to the ZAI
	
	TMatrixT<double>	fDecayMatrix;	//!< Matrix with half life for each nuclei
	
	CLASSNucleiFiliation	fFastDecay;	//!< Store the nuclei being cut (HL threshold)
	CLASSNucleiFiliation	fNormalDecay;	//!< Store the uncut nuclei
	IsotopicVector			fDecayConstante; //!< List of decay constants
	
	CLASSNucleiFiliation	fSpontaneusYield;	//!< Store the spontaneus fission yield
	CLASSNucleiFiliation	fReactionYield;		//!< Store the reaction fission yield

	CLASSNucleiFiliation	fCaptureReaction;	//!< Store the reaction capture Filiation
	CLASSNucleiFiliation	fn2nReaction;		//!< Store the reaction n,2n Filiation
	
	string	fSpontaneusYieldFile;	//!< Store the name of the spontaneus fission yield file
	string	fReactionYieldFile;		//!< Store the name of the reaction fission yield file

	string	fSpectrumType;			//!< Type of the spectrum : thermal or fast. (needed for Isomeric branching ratios)
	
	//{
	/// Return the Fission XS Matrix at the time TStep
	/*!
	 // This method extract the fission cross section of an EvolutionData at the set time
	 // \param EvolutionDataStep : EvolutionData
	 // \param TStep : time
	 */
	TMatrixT<double> GetFissionXsMatrix(EvolutionData EvolutionDataStep,double TStep);
	//}
	
	//{
	/// Return the capture cross section matrix at the time TStep
	/*!
	 // This Method extract the capture cross section of an EvolutionData at the set time
	 // \param EvolutionDataStep : EvolutionData
	 // \param TStep :  time
	 */
	TMatrixT<double> GetCaptureXsMatrix(EvolutionData EvolutionDataStep,double TStep);
	//}
	
	//{
	/// Return the n2n XS matrix at the time TStep
	/*!
	 // This Method extract the (n,2n) Cross section of an EvolutionData at the set time
	 // \param EvolutionDataStep : EvolutionData
	 // \param TStep :  time
	 */
	TMatrixT<double> Getn2nXsMatrix(EvolutionData EvolutionDataStep,double TStep);
	//}
	
	
	//{
	//! Returns a particular decay mode.
	/*!
	 \param DecayModes : a list of decay modes with their branching ratios and isomeric state of the Daughters.
	 \param BR : branching ratio of the current decay mode
	 \param Iso : isomeric state of the Daughter of the current decay mode.
	 \param StartPos : the current decay mode to extract.
	 */
	string GetDecay(string DecayModes, double &BR,int &Iso, int &StartPos);
	//}
	
	CLASSNucleiFiliation ReadFPYield(string Yield);	///< Read a CLASSYield file and return the correpsponding map
	
	private :
 	
};

#endif

