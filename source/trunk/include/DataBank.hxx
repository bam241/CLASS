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
class DataBank : public CLSSObject
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

//********* Set Method *********//
	
	void SetDataBank(map<T ,EvolutionData > mymap)	{ fDataBank = mymap; } 

	void SetDataBaseIndex(string database) { fDataBaseIndex = database; }
	EvolutionData GenerateEvolutionData(IsotopicVector isotopicvector, double cycletime, double Power); //!< Genration of a New EvolutionData From the one already present
	EvolutionData OldGenerateEvolutionData(IsotopicVector isotopicvector, double cycletime, double Power); //!< Genration of a New EvolutionData From the one already present
	void SetUpdateReferenceDBatEachStep(bool val)	{fUpdateReferenceDBatEachStep = val;}

	void SetOldReadMethod(bool val)			{ fOldReadMethod = val;}
	void SetFissionEnergy(string FissionEnergyFile);
	void SetFissionEnergy(ZAI zai, double E);
	void SetFissionEnergy(int Z, int A, int I, double E )   { SetFissionEnergy(ZAI(Z,A,I), E);}

//********* Modification Method *********//
	
	IsotopicVector	Evolution(const T &key, double dt);	///< Return the Product IsotopicVector evolution from zai during a dt time
	void	ReadDataBase();				///< ...
	void	CalculateDistanceParameter();		///< Calcul of the weight for each ZAI in the distance calculation from the mean XS of the DataBank
	void	SetDistanceParameter(IsotopicVector DistanceParameter);///< Define mannually the weight for each ZAI in the distance calculation
	void	SetDistanceType(int DistanceType);	///< Define the way to decide if two isotopic vectors are close :
								///< 0 is for the standard norme,
								///< 1 for each ZAI weighted with its XS,
								///< 2 for each ZAI weighted with coefficient given by the user.

	void UseDBTimeStep(bool oldmethod = true)		{fUseDBTimeStep = oldmethod;}
	
	void UseOldGeneration(bool oldmethod = true)		{fUseOldGeneration = oldmethod;}
//********* Printing Method *********//
	void Print() const;
	
protected :

	map<T, EvolutionData>	fDataBank;
	map<T, EvolutionData>	fDataBankCalculated;
	
 	string			fDataBaseIndex;

	bool			fUpdateReferenceDBatEachStep;
	bool			fOldReadMethod;
	bool			fUseOldGeneration;
	bool			fUseDBTimeStep;
	
 	string 			fFuelType;
 	pair<double,double>	fBurnUpRange;
 	vector<double>		fFuelParameter;
	int 	fDistanceType;		///< 0 is for the standard norm (Default = 0),
					///< 1 for each ZAI weighted with its XS,
					///< 2 for each ZAI weighted with coefficient given by the user.
	
	IsotopicVector		fDistanceParameter;	///< weight for each ZAI in the distance calculation
	
	TMatrixT<double>	fDecayMatrix;		///< Matrix with half life of each nuclei
	void	BuildDecayMatrix();
	TMatrixT<double> GetFissionXsMatrix(EvolutionData EvolutionDataStep,double TStep);
	TMatrixT<double> GetCaptureXsMatrix(EvolutionData EvolutionDataStep,double TStep);
	TMatrixT<double> Getn2nXsMatrix(EvolutionData EvolutionDataStep,double TStep);
	
	
	
	TMatrixT<double> ExtractXS(EvolutionData EvolutionDataStep,double TStep);

	map<ZAI, double >		fFissionEnergy; ///< Store the Energy per fission use for the flux normalisation.
	map<ZAI, map<ZAI, double> >	fFastDecay;
	map<ZAI, int> findex_inver;
	map<int, ZAI> findex;
	//0 TMP
	//1 PF
	//2 232Th
	//3 233U
	//4 234U
	//5 235U
	//6 236U
	//7 238U
	//8 237Np
	//9 238Pu
	//10 239Pu
	//11 240Pu
	//12 241Pu
	//13 242Pu
	//14 241Am
	//15 242Am*
	//16 243Am
	//17 242Cm
	//18 243Cm
	//19 244Cm
	//20 245Cm
	//21 246Cm
	//22 247Cm
	//23 248Cm

};



#endif
