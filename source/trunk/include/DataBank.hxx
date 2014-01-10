#ifndef __DataBank_HXX__
#define __DataBank_HXX__

/*!
 \file
 \brief Header file for DataBank class. 
 The aim of this Class is to store the evolution Database of the all decay nuclei.

 @author BaM, Marc
 @version 2.0
 */

#include "CLSSObject.hxx"
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


template <class T> 
class DataBank : public CLSSObject, DynamicalSystem
{

public :
//********* Constructor/Destructor Method *********//
	///< Normal Constructor.
	DataBank();
	
	DataBank(LogFile* Log, string DB_index_file, bool setlog = true, bool olfreadmethod = true );
	
	///< Normal Destructor.
	~DataBank();

//********* Get Method *********//
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




//********* Set Method *********//
	
	void SetDataBank(map<T ,EvolutionData > mymap)	{ fDataBank = mymap; }

	void SetDataBaseIndex(string database) { fDataBaseIndex = database; }
	EvolutionData GenerateEvolutionData(IsotopicVector isotopicvector, double cycletime, double Power); //!< Genration of a New EvolutionData From the one already present

	void SetOldReadMethod(bool val)			{ fOldReadMethod = val;}			// use the old reading method
	void SetFissionEnergy(string FissionEnergyFile);						// set Fission Energy using a file
	void SetFissionEnergy(ZAI zai, double E);							// set Fission Energy for a ZAI
	void SetFissionEnergy(int Z, int A, int I, double E )   { SetFissionEnergy(ZAI(Z,A,I), E);}	// set Fission Energy for a ZAI

	void SetDataFileName(string name)		{ fDataFileName = name;}		// Set the name of the reaction file
	void SetDataDirectoryName(string name)		{ fDataDirectoryName = name;}		// Set the Path to the reaction file
	void SetShartestHalfLife(double halflife)	{ fShorstestHalflife = halflife;}	// Set the Half Life cut
	void LoadFPYield(string SponfaneusYield, string ReactionYield);			//Build Fision Yields maps;


//********* Modification Method *********//
	
	IsotopicVector	Evolution(const T &key, double dt);	///< Return the Product IsotopicVector evolution from zai during a dt time
	void	ReadDataBase();				///< ...
	void	CalculateDistanceParameter();		///< Calcul of the weight for each ZAI in the distance calculation from the mean XS of the DataBank
	void	SetDistanceParameter(IsotopicVector DistanceParameter);///< Define mannually the weight for each ZAI in the distance calculation
	void	SetDistanceType(int DistanceType);	///< Define the way to decide if two isotopic vectors are close :
								///< 0 is for the standard norme,
								///< 1 for each ZAI weighted with its XS,
								///< 2 for each ZAI weighted with coefficient given by the user.

	void	BuildDecayMatrix();

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

    
    
//********* Printing Method *********//
	void Print() const;
	
protected :

	double  fShorstestHalflife;
	int	fZAIThreshold;	//!< Highest Mass deal bye the evolution (default 90)

	string			fDataFileName;
	string			fDataDirectoryName;

	map<T, EvolutionData>	fDataBank;
	map<T, EvolutionData>	fDataBankCalculated;
	
 	string			fDataBaseIndex;

	bool			fUseRK4EvolutionMethod;
	bool			fOldReadMethod;

 	string 			fFuelType;
 	pair<double,double>	fBurnUpRange;
 	vector<double>		fFuelParameter;
	int 	fDistanceType;		///< 0 is for the standard norm (Default = 0),
					///< 1 for each ZAI weighted with its XS,
					///< 2 for each ZAI weighted with coefficient given by the user.
	
	IsotopicVector		fDistanceParameter;	///< weight for each ZAI in the distance calculation
	

	TMatrixT<double> GetFissionXsMatrix(EvolutionData EvolutionDataStep,double TStep);
	TMatrixT<double> GetCaptureXsMatrix(EvolutionData EvolutionDataStep,double TStep);
	TMatrixT<double> Getn2nXsMatrix(EvolutionData EvolutionDataStep,double TStep);
	
	TMatrixT<double> ExtractXS(EvolutionData EvolutionDataStep,double TStep);

	string GetDecay(string DecayModes, double &BR,int &Iso, int &StartPos);

	TMatrixT<double>		fDecayMatrix;		///< Matrix with half life of each nuclei
	map<ZAI, double >		fFissionEnergy; ///< Store the Energy per fission use for the flux normalisation.
	map<ZAI, map<ZAI, double> >	fFastDecay;
	map<ZAI, IsotopicVector>	fSpontaneusYield;
	map<ZAI, IsotopicVector>	fReactionYield;


	double	*fTheNucleiVector;	//!< The evolving atoms copied from Material proportions.
	double 	**fTheMatrix;  		//!< The evolution Matrix

	int	fNVar;		 //!< The size of the composition vector and /or number of ZAIs involved.
	double	fPrecision;	//!< Precision of the RungeKutta
	double	fHestimate;	//!< RK Step estimation.
	double	fHmin;		//!< RK minimum Step.
	double	fMaxHdid;	//!< store the effective RK max step
	double	fMinHdid;	//!< store the effective RK min step
	bool	fIsNegativeValueAllowed; //!< whether or not negative value are physical.

	map<ZAI, int> findex_inver;
	map<int, ZAI> findex;

	map< ZAI,IsotopicVector > ReadFPYield(string Yield);


};



#endif
