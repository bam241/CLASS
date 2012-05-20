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

class ZAI;
class IsotopicVector;
class TGraphErrors;
class LogFile;

using namespace std;

///< A ZAIIDataBase defined a database which contain the evolution of faction of all product, subproduct (or sub(sub...sub)product) for a nucleus. 
/*!
 The aim of this class is to handle the evolution all information of all Nucleus product, subproduct (or sub(sub...sub)product) as the fonction of the time.
 
 @author BaM
 @version 1.0
 */


class EvolutiveProduct 
{
	
public :

//********* Constructor/Destructor Method *********//
	///< Normal ZAI DB Constructor.
	EvolutiveProduct(LogFile* Log, int A = 0, int Z = 0, int I = 0 , string DBindexfile = "Default_Index.dat"); 	///< Make a new ZAI evolution 
	
	EvolutiveProduct(string DBindexfile = "Default_Reactor.dat");						///< Make a new Reactor evolution
	 
	///< Normal Destructor.
	~EvolutiveProduct();
	
//********* Get Method *********//
	map<ZAI ,TGraphErrors* >	GetEvolutiveProduct() const {return fEvolutiveProduct;}
	TGraphErrors*	GetEvolutionTGraphErrors(const ZAI& zai); 
								///< Return the A,Z product proportion evolution TGraphErrors
	IsotopicVector	GetIsotopicVectorAt(double t); 		///< Return the Product IsotopicVector evolution TGraphErrors
	

protected :
	map<ZAI ,TGraphErrors* >	fEvolutiveProduct;	///< 
	void		ReadDB(string DBfile);
	int 		fDatabaseEndTime;
	double		Interpolate(double t, TGraphErrors& EvolutionGraph); 
								///< Interpolating the value of EvolutionGraph at the t time
	void		AddAsStable(int Z, int A, int I=0);
	LogFile*	fLog;

};

#endif
