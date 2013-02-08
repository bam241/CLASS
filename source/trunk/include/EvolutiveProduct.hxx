#ifndef __EVOLUTIVEPRODUCT_HXX__
#define __EVOLUTIVEPRODUCT_HXX__

/*!
 \file
 \brief Header file for EvolutiveProduct classes. 
  The aim of this Class is to manage evolution of a system, such as a radioactiv nuclei or a reactor. It store the evolution of daughter nuclei proportion as the function of time of the 
 
 @author BaM
 @version 0.
 */

#include <string>
#include <map>
#include "IsotopicVector.hxx"

#include "ZAI.hxx"
#include "TObject.h"
#include "TMatrix.h"

class TGraph;
class EvolutiveProduct;
class LogFile;

using namespace std;
typedef long long int cSecond;

///< A ZAIIDataBase defined a database which contain the evolution of faction of all product, subproduct (or sub(sub...sub)product) for a nucleus. 
/*!
 The aim of this class is to handle the evolution all information of all Nucleus product, subproduct (or sub(sub...sub)product) as the fonction of the time.
 
 @author BaM
 @version 1.0
 */

EvolutiveProduct operator*(EvolutiveProduct const& evol, double F);
EvolutiveProduct operator*(double F, EvolutiveProduct const& evol);




class EvolutiveProduct : public TObject
{
	
public :

//********* Constructor/Destructor Method *********//
	///< Normal DB Constructor.
	
	EvolutiveProduct();  
	EvolutiveProduct(LogFile* Log); 	///< Make a new Evolutive Product evolution 
	EvolutiveProduct(LogFile* Log, string DB_file, bool stable = false, ZAI zai = ZAI(0,0,0) ); 	///< Make a new Evolutive Product evolution 
	 
	///< Normal Destructor.
	~EvolutiveProduct();

//********* Set Method *********//
	void 	SetReactorType(string reactortype)	{fReactorType = reactortype;}
	void	SetFuelType(string fueltype)		{fFuelType = fueltype;}
	void 	SetPower(double power)			{fPower = power;}
	void 	SetHMMass(double HMMass)		{fHMMass = HMMass;}

	
//********* Get Method *********//
	map<ZAI ,TGraph* >	GetEvolutiveProduct()	const { return fEvolutiveProduct; }	//!<
	map<ZAI ,TGraph* >	GetFissionXS()		const { return fFissionXS; }		//!<
	map<ZAI ,TGraph* >	GetCaptureXS()		const { return fCaptureXS; }		//!<
	map<ZAI ,TGraph* >	Getn2nXS()		const { return fn2nXS; }		//!<
	TGraph*			GetKeff()		const { return fKeff; }
	TGraph*			GetFlux()		const { return fFlux; }

	double			GetPower()		const { return fPower; }			//!<
	string			GetDB_file()		const { return fDB_file; }

	TGraph*	GetEvolutionTGraph(const ZAI& zai); 
								///< Return the A,Z product proportion evolution TGraph
	IsotopicVector	GetIsotopicVectorAt(double t); 		///< Return the Product IsotopicVector evolution TGraph




	bool Insert(pair<ZAI, TGraph*> zaitoinsert);			//!<

//********* Get Method *********//
	EvolutiveProduct GenerateDBFor(IsotopicVector isotopicvector);	///< Build A DB from a close one


protected :
	
	string fDB_file;
	map<ZAI ,TGraph* >	fEvolutiveProduct;	//!< 
	map<ZAI ,TGraph* >	fFissionXS;	//!< 
	map<ZAI ,TGraph* >	fCaptureXS;	//!< 
	map<ZAI ,TGraph* >	fn2nXS;	//!< 
	TGraph*		fKeff;
	TGraph*		fFlux;
	
	cSecond 	fDBendTime;
	bool		fIsCrossSection;
	
	
	
	string 		fReactorType;
	string		fFuelType;
	double 		fPower;
	double 		fHMMass;
	
	void		ReadDB(string DBfile);
	double		Interpolate(double t, TGraph& EvolutionGraph); 
								///< Interpolating the value of EvolutionGraph at the t time
	void		AddAsStable(ZAI zai);
	LogFile*	fLog;

	ClassDef(EvolutiveProduct,0);
};

#endif
