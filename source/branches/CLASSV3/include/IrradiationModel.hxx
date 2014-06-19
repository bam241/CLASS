#ifndef _IRRADIATIONMODEL_HXX
#define _IRRADIATIONMODEL_HXX


/*!
 \file
 \brief Header file for IrradiationModel class.
 
 
 @author BaM
 @version 2.0
 */


#include "CLASSObject.hxx"
#include "TMatrix.h"
#include "IsotopicVector.hxx"
#include "EvolutionData.hxx"

#include <map>
#include <vector>


using namespace std;
typedef long long int cSecond;

class ZAI;
class LogFile;
//-----------------------------------------------------------------------------//
/*!
 Define a IrradiationModel.
 The aim of these class is synthetyse all the commum properties to all Irradiation Model.
 
 
 @author BaM
 @version 3.0
 */
//________________________________________________________________________

class EvolutionData;

class IrradiationModel : public CLASSObject
{
	
	public :
	
	IrradiationModel();

	IrradiationModel(LogFile* log);

	/// virtueal method called to perform the irradiation calculation using a set of cross section.
	/*!
	 Perform the Irradiation Calcultion using the XSSet data
	 \param IsotopicVector IV isotopic vector to irradiate
	 \param EvolutionData XSSet set of corss section to use to perform the evolution calculation
	 \param double Power, constant power to use for irradation
	 \param double irradiationtime, time of the irradiation
	 */
	virtual	 EvolutionData GenerateEvolutionData(IsotopicVector IV, EvolutionData XSSet, double Power, double cycletime) { return 0;}
	//}

	
	
	
	//********* Get Method *********//
	/*!
	 \name Get Method
	 */
	//@{
	string	GetDataFileName()	const { return fDataFileName; }
	string	GetDataDirectoryName()  const { return fDataDirectoryName; }

	double  GetShorstestHalflife()	const { return fShorstestHalflife; }
	
	//@}
	
	
	
	
	//********* Set Method *********//
	
	/*!
	 \name Set Method
	 */
	//@{
	

	//{
	/// set Fission Energy using a file
	/*!
	 // This method fill the Fission Energy map using a file
	 // \param FissionEnergyFile: filename containing the Fission Energy of some nuclei (form : Z A I Energy)
	 */
	void SetFissionEnergy(string FissionEnergyFile);
	//}
	
	//{
	/// set Fission Energy for a ZAI using ZAI(Z,A,I)
	/*!
	 // This method fill the Fission Energy map of a set ZAI
	 // \param zai ZAI
	 // \param E Fission energy of the ZAI
	 */
	void SetFissionEnergy(ZAI zai, double E);
	//}
	
	//{
	/// set Fission Energy for a ZAI using the Z, A, I
	/*!
	 // This method fill the Fission Energy map of a set ZAI
	 // \param Z Z of the ZAI
	 // \param A A of the ZAI
	 // \param I I of the ZAI
	 // \param E Fission energy of the ZAI
	 */
	void SetFissionEnergy(int Z, int A, int I, double E )   { SetFissionEnergy(ZAI(Z,A,I), E);}
	//}
	
	void SetShortestHalfLife(double halflife)	{ fShorstestHalflife = halflife;}	///< Set the Half Life cut
	void LoadFPYield(string SponfaneusYield, string ReactionYield);				///< Build Fision Yields maps;
	
	
	
	
	
	
	//********* Evolution Method *********//
	
	//@}
	/*!
	 \name Evolution Method
	 */
	//@{
	

	void	BuildDecayMatrix();			///w Build the Decay Matrix for the futur evolution...
	
	//@}
	

	//********* Other Method *********//
	/*!
	 \name Other Method
	 */
	//@{
	void Print() const;
	
	//@}
	
	
	
	
	
	protected :
	
	double  fShorstestHalflife;
	int	 fZAIThreshold;	//!< Lowest Mass deal by the evolution (default 90)
	string			fDataFileName;		///< Name of the decay list
	string			fDataDirectoryName;	///< Path to the decay list file

	
	TMatrixT<double>		fDecayMatrix;	///< Matrix with half life of each nuclei
	map<ZAI, double >		fFissionEnergy;	///< Store the Energy per fission use for the flux normalisation.
	map<ZAI, map<ZAI, double> >	fFastDecay;	///< Store the cut decay
	map<ZAI, IsotopicVector>	fSpontaneusYield;	///< Store the Spontaneus fission yield
	map<ZAI, IsotopicVector>	fReactionYield;		///< Store the reaction fission yield
	

	
	map<ZAI, int> findex_inver;	///< correspondance matrix from ZAI to the column (or line) of the different Reaction/Decay matrix
	map<int, ZAI> findex;		///< correspondance matrix from the column (or line) of the different Reaction/Decay matrix to the ZAI
	
	//{
	/// Return the Fission XS Matrix at the time TStep
	/*!
	 // This Method extract the Fission Cross section of an EvolutionData at the set time
	 // \param EvolutionDataStep: EvolutionData
	 // \param TStep:  time
	 */
	TMatrixT<double> GetFissionXsMatrix(EvolutionData EvolutionDataStep,double TStep);
	//}
	
	//{
	/// Return the Capture XS Matrix at the time TStep
	/*!
	 // This Method extract the capture Cross section of an EvolutionData at the set time
	 // \param EvolutionDataStep: EvolutionData
	 // \param TStep:  time
	 */
	TMatrixT<double> GetCaptureXsMatrix(EvolutionData EvolutionDataStep,double TStep);
	//}
	
	//{
	/// Return the n2n XS Matrix at the time TStep
	/*!
	 // This Method extract the (n,2n) Cross section of an EvolutionData at the set time
	 // \param EvolutionDataStep: EvolutionData
	 // \param TStep:  time
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
	
	map< ZAI,IsotopicVector > ReadFPYield(string Yield);	///< Read a CLASSYield file and return the correpsponding map
	
	private :
 	
};

#endif

