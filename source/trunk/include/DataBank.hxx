#ifndef __DataBank_HXX__
#define __DataBank_HXX__

/*!
 \file
 \brief Header file for DataBank class. 
 @version 2.0
 */

#include "CLASSObject.hxx"
#include "TMatrix.h"
#include "IsotopicVector.hxx"
#include "DynamicalSystem.hxx"

#include <map>
#include <vector>


using namespace std;
typedef long long int cSecond;

class ZAI;
class EvolutionData;
class LogFile;

double ReactionRateWeightedDistance(IsotopicVector IV1, EvolutionData DB );
double ReactionRateWeightedDistance(EvolutionData DB, IsotopicVector IV1  );

//-----------------------------------------------------------------------------//
/*!
 Define a DataBank.
 The aim of these class is describe the evolution of "all" evoluting system in CLASS.
 2 kind of evoluting system can be defined :
	\li the Decay Matrix
	\li the fuel 
 
 For the Decay Matrix the Databank take the form of the ZAI template (ie DataBank<ZAI>) which mainly contain a map of <ZAI,EvolutionData>.This map do the correspondance between a ZAI and its decay evolution (containing all the daughter nuclei comming from the decay a the original ZAI).
 
 For the Fuel the Databank take the form of the IsotopicVector template (ie DataBank<Isotopic>), which mainly contain a map of <IsotopicVector,EvolutionData>. This map do the correspondance between a IsotopicVector and its evolution throw irradiation (containing all the nuclei produced by the reaction on the original IsotopicVector)

 @author BaM
 @author Marc
 @author PTO for a part the Decay management -- steal from MURE (Even if he does not kown it!! :))
 @version 2.0
 */
//________________________________________________________________________



template <class T> 
class DataBank : public CLASSObject, DynamicalSystem
{

public :


//********* Constructor/Desctructor *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	/// Normal Constructor.
	DataBank();

	//{
	/// Special Constructor.
	/*!
	 Use to load a LogFile
	 \param LogFile LogFile used for the log...
	 \param DB_index_file path to the index file
	 \param setlog if the log are stored in the LogFile
	 \param olfreadmethod true if the old format of EvolutionData are used (ie without the key word such as Inv, XSFiss...)
	 */
	DataBank(LogFile* Log, string DB_index_file, bool setlog = true, bool olfreadmethod = true );
	//}

	//{
	/// Normal Destructor.
	/*!
		Delete de DataBank and all associated EvolutionData...
	 */
	~DataBank();
	//}

	//{
	/// Reset the DataBank.
	/*!
	 Use to reset the DataBank to its default values  whihout deleting the EvolutionData (which contain pointer... ).
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
	map<T ,EvolutionData >		GetDataBank()		const	{ return fDataBank; }		//!< Return the DataBank
	string 				GetDataBaseIndex()	const	{ return fDataBaseIndex; }	//!< Return the index Name
	string				GetFuelType()		const	{ return fFuelType; }		//!< Return the fuel type of the DB
	vector<double>			GetFuelParameter()	const	{ return fFuelParameter; }	//!< Return the Fuel parameter of the DB
	pair<double,double>		GetBurnUpRange()	const	{ return fBurnUpRange;}		//!< Return the BurnUp range of the DB
	bool 				IsDefine(const T& key)	const;					//!< True the key is define, false unstead

	map<double, EvolutionData>	GetDistancesTo(IsotopicVector isotopicvector, double t = 0) const;	//! Return a map containing the distance of each EvolutionData in the DataBase to the set IV at the t time
	EvolutionData	GetClosest(IsotopicVector isotopicvector, double t = 0) const;	//! Return the closest EvolutionData from the DataBank.


	string	GetDataFileName()	const { return fDataFileName; }
	string	GetDataDirectoryName()  const { return fDataDirectoryName; }

	double  GetShorstestHalflife()	const { return fShorstestHalflife; }

	//@}




//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{

	void SetDataBank(map<T ,EvolutionData > mymap)	{ fDataBank = mymap; }	//!< Set the Databank map

	void SetDataBaseIndex(string database) { fDataBaseIndex = database;; ReadDataBase(); }	//!< Set the Name of the database index

	EvolutionData GenerateEvolutionData(IsotopicVector isotopicvector, double cycletime, double Power); //!< Generation of a New EvolutionData From the one already present

	void SetOldReadMethod(bool val)			{ fOldReadMethod = val; ReadDataBase();}			///< use the old reading method

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

	void SetDataFileName(string name)		{ fDataFileName = name;}		///< Set the name of the reaction file
	void SetDataDirectoryName(string name)		{ fDataDirectoryName = name;}		///< Set the Path to the reaction file
	void SetShortestHalfLife(double halflife)	{ fShorstestHalflife = halflife;}	///< Set the Half Life cut
	void LoadFPYield(string SponfaneusYield, string ReactionYield);				///< Build Fision Yields maps;

	void SetWeightedDistanceCalculation(bool val = true) { fWeightedDistance = val;}		///< Set weighted Distance calculation
	void SetEvolutionDataInterpolation(bool val = true) { fEvolutionDataInterpolation = val;}		///< Set weighted Distance calculation

	void SetDistanceParameter(IsotopicVector DistanceParameter);		///< Define mannually the weight for each ZAI in the distance calculation


	//{
	/// Define the way to decide if two isotopic vectors are close.
	/*!
	// The different algorythm are:
	// \li 0 is for the standard norme,
	// \li 1 for each ZAI weighted with its XS,
	// \li 2 for each ZAI weighted with coefficient given by the user.
	*/
	 void SetDistanceType(int DistanceType);
	//}


	

//********* Evolution Method *********//

	//@}
	/*!
	 \name Evolution Method
	 */
	//@{


	IsotopicVector	Evolution(const T &key, double dt);	///< Return the Product IsotopicVector evolution from zai during a dt time
	void	CalculateDistanceParameter();		///< Calcul of the weight for each ZAI in the distance calculation from the mean XS of the DataBank
	void	BuildDecayMatrix();			///w Build the Decay Matrix for the futur evolution...


	

//********* RK4 Method *********//

	//@}
	/*!
	 \name RK4 Method
	 */
	//@{

	void	UseRK4EvolutionMethod(bool usemethod = true)	{fUseRK4EvolutionMethod = usemethod;}

	
	using	DynamicalSystem::RungeKutta;
	//!	Pre-treatment Runge-Kutta method.
	/*!
	// This method does initialisation and then call DynamicalSystem::RungeKutta
	// \param t1: initial time
	// \param t2: final time
	*/


   	void BuildEqns(double t, double *N, double *dNdt);
	void SetTheMatrixToZero();			//!< Initialize the evolution Matrix
	void ResetTheMatrix();
	void SetTheMatrix(TMatrixT<double> BatemanMatrix);	//!< Set the Evolution Matrix (Bateman equations)
	TMatrixT<double> GetTheMatrix();		//!< return the Evolution Matrix (Bateman equations)

	void SetTheNucleiVectorToZero();			//!< Initialize the evolution Matrix
	void ResetTheNucleiVector();
	void SetTheNucleiVector(TMatrixT<double> NEvolutionMatrix);	//!< Set the Evolution Matrix (Bateman equations)
	TMatrixT<double> GetTheNucleiVector();		//!< return the Evolution Matrix (Bateman equations)
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

	map<T, EvolutionData>	fDataBank;		///< DataBanck map
	map<T, EvolutionData>	fDataBankCalculated;	///< Map of the already calculated EvolutionData (to avoid recalculation...)
	
 	string			fDataBaseIndex;			///< Name of the index

	bool			fUseRK4EvolutionMethod; ///< if true use RK4 calculation, mtriciel unstead
	bool			fOldReadMethod;		///< use old DB format
	bool			fWeightedDistance;	///< USe XS weighted distance calculation
	bool			fEvolutionDataInterpolation;	///< USe XS weighted distance calculation


 	string 			fFuelType;		///< Type of fuel of the DataBank
 	pair<double,double>	fBurnUpRange;		///< Range of the Burn-up range of the DataBank
 	vector<double>		fFuelParameter;		///< Parameter needed by the equivalence model



	int 	fDistanceType;		///< Set the distance calculation algorytm
					/// \li 0 is for the standard norm (Default = 0),
					/// \li 1 for each ZAI weighted with its XS,
					/// \li 2 for each ZAI weighted with coefficient given by the user.
	
	IsotopicVector		fDistanceParameter;	///< weight for each ZAI in the distance calculation

	TMatrixT<double>		fDecayMatrix;	///< Matrix with half life of each nuclei
	map<ZAI, double >		fFissionEnergy;	///< Store the Energy per fission use for the flux normalisation.
	map<ZAI, map<ZAI, double> >	fFastDecay;	///< Store the cut decay
	map<ZAI, IsotopicVector>	fSpontaneusYield;	///< Store the Spontaneus fission yield
	map<ZAI, IsotopicVector>	fReactionYield;		///< Store the reaction fission yield


	double	*fTheNucleiVector;	//!< The evolving atoms copied from Material proportions.
	double 	**fTheMatrix;  		//!< The evolution Matrix
	int	fNVar;		 //!< The size of the composition vector and /or number of ZAIs involved.

	double	fPrecision;	//!< Precision of the RungeKutta
	double	fHestimate;	//!< RK Step estimation.
	double	fHmin;		//!< RK minimum Step.
	double	fMaxHdid;	//!< store the effective RK max step
	double	fMinHdid;	//!< store the effective RK min step
	bool	fIsNegativeValueAllowed; //!< whether or not negative value are physical.

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


};



#endif
