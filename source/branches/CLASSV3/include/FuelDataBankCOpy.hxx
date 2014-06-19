#ifndef __FuelDataBank_HXX__
#define __FuelDataBank_HXX__

/*!
 \file
 \brief Header file for FuelDataBank class.
 @version 2.0
 */

#include "CLASSObject.hxx"
#include "TMatrix.h"
#include "IsotopicVector.hxx"
#include "DynamicalSystem.hxx"
#include "EvolutionData.hxx"

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
 Define a FuelDataBank.
 The aim of these class is describe the evolution of fuel evolution in CLASS.

 \li the fuel


 For the Fuel the FuelDataBank take the form of the IsotopicVector template which mainly contain a map of <IsotopicVector,EvolutionData>. This map do the correspondance between a IsotopicVector and its evolution throw irradiation (containing all the nuclei produced by the reaction on the original IsotopicVector)

 @author BaM
 @author Marc
 @author PTO for a part the Decay management -- steal from MURE (Even if he does not kown it!! :))
 @version 2.0
 */
//________________________________________________________________________



class FuelDataBank : public CLASSObject, DynamicalSystem
{

	public :


	//********* Constructor/Desctructor *********//

	/*!
	 \name Constructor/Desctructor
	 */
	//@{
	/// Normal Constructor.
	FuelDataBank();

	//{
	/// Special Constructor.
	/*!
	 Use to load a LogFile
	 \param LogFile LogFile used for the log...
	 \param DB_index_file path to the index file
	 \param setlog if the log are stored in the LogFile
	 \param olfreadmethod true if the old format of EvolutionData are used (ie without the key word such as Inv, XSFiss...)
	 */
	FuelDataBank(LogFile* Log, string DB_index_file, bool setlog = true, bool olfreadmethod = true );
	//}

	//{
	/// Normal Destructor.
	/*!
	 Delete de FuelDataBank and all associated EvolutionData...
	 */
	~FuelDataBank();
	//}

	//{
	/// Reset the FuelDataBank.
	/*!
	 Use to reset the FuelDataBank to its default values  whihout deleting the EvolutionData (which contain pointer... ).
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
	map<IsotopicVector ,EvolutionData >	GetFuelDataBank()	const	{ return fFuelDataBank; }	//!< Return the FuelDataBank
	string 				GetDataBaseIndex()	const	{ return fDataBaseIndex; }	//!< Return the index Name
	string				GetFuelType()		const	{ return fFuelType; }		//!< Return the fuel type of the DB
	vector<double>			GetFuelParameter()	const	{ return fFuelParameter; }	//!< Return the Fuel parameter of the DB
	pair<double,double>		GetBurnUpRange()	const	{ return fBurnUpRange;}		//!< Return the BurnUp range of the DB
	bool 				IsDefine(IsotopicVector IV)	const;					//!< True the key is define, false unstead

	map<double, EvolutionData>	GetDistancesTo(IsotopicVector isotopicvector, double t = 0) const;	//! Return a map containing the distance of each EvolutionData in the DataBase to the set IV at the t time
	EvolutionData	GetClosest(IsotopicVector isotopicvector, double t = 0) const;	//! Return the closest EvolutionData from the FuelDataBank.



	//@}




	//********* Set Method *********//

	/*!
	 \name Set Method
	 */
	//@{

	void SetFuelDataBank(map< IsotopicVector ,EvolutionData > mymap)	{ fFuelDataBank = mymap; }	//!< Set the FuelDataBank map

	void SetDataBaseIndex(string database) { fDataBaseIndex = database;; ReadDataBase(); }	//!< Set the Name of the database index

	EvolutionData GenerateEvolutionData(IsotopicVector isotopicvector, double cycletime, double Power); //!< Generation of a New EvolutionData From the one already present

	void SetOldReadMethod(bool val)			{ fOldReadMethod = val; ReadDataBase();}			///< use the old reading method



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


	void	CalculateDistanceParameter();		///< Calcul of the weight for each ZAI in the distance calculation from the mean XS of the FuelDataBank








	//********* Other Method *********//
	/*!
	 \name Other Method
	 */
	//@{
	void	ReadDataBase();				///< read the index file and fill the evolutionData map

	void Print() const;

	//@}





	protected :


	map<IsotopicVector, EvolutionData>	fFuelDataBank;		///< DataBanck map
	map<IsotopicVector, EvolutionData>	fFuelDataBankCalculated;	///< Map of the already calculated EvolutionData (to avoid recalculation...)

 	string			fDataBaseIndex;			///< Name of the index

	bool			fUseRK4EvolutionMethod; ///< if true use RK4 calculation, mtriciel unstead
	bool			fOldReadMethod;		///< use old DB format
	bool			fWeightedDistance;	///< USe XS weighted distance calculation
	bool			fEvolutionDataInterpolation;	///< USe XS weighted distance calculation


 	string 			fFuelType;		///< Type of fuel of the FuelDataBank
 	pair<double,double>	fBurnUpRange;		///< Range of the Burn-up range of the FuelDataBank
 	vector<double>		fFuelParameter;		///< Parameter needed by the equivalence model



	int 	fDistanceType;		///< Set the distance calculation algorytm
	/// \li 0 is for the standard norm (Default = 0),
	/// \li 1 for each ZAI weighted with its XS,
	/// \li 2 for each ZAI weighted with coefficient given by the user.

	IsotopicVector		fDistanceParameter;	///< weight for each ZAI in the distance calculation


};



#endif
