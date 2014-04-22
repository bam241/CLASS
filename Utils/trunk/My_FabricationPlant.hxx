#ifndef __MyFabricationPlant_HXX__
#define __MyFabricationPlant_HXX__

/*!
 \file 
 \brief Header file for FabricationPlant class.

 The aim of the Class is to manage evolution of FabricationPlant

 
 @author BaM, Marc
 @version 2.0
 */



#include <vector>
#include <map>

#include "CLSSFacility.hxx"
#include "IsotopicVector.hxx"
#include "EvolutionData.hxx"
#include "CLASS.hxx"
#include "Storage.hxx"
#include "Reactor.hxx"
#include "DataBank.hxx"
#include "LogFile.hxx"
#include "ZAI.hxx"

using namespace std;
typedef long long int cSecond;

class MyFabricationPlant : public FabricationPlant
{
	//on utilise les constructeur de la classe FabricationPlant
	// dans le code :MyFabricationPlant *FP_MOX = new FabricationPlant(gCLASS.GetLog(),Stock, ReUsable);
	//on a besoin de surcherger le constructeur que si
public :
	MyFabricationPlant();
	MyFabricationPlant(LogFile* log);
	
	MyFabricationPlant(LogFile* log, Storage* storage, Storage* reusable, double fabricationtime = 365.25*24*3600*2);
	///< Normal Destructor.
	~MyFabricationPlant();
	
	void	BuildFuelForReactor(int ReactorId);			//obligatoire, c'est un peu la base
protected :
//********* Private Method *********//
// mettre ici les méthodes dont on a besoin dans sa fonction BuildFuelForReactor perso

};

#endif
