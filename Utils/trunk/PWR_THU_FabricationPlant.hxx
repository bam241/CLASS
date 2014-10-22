#ifndef __PWR_THU_FabricationPlant_HXX__
#define __PWR_THU_FabricationPlant_HXX__

/*!
 \file 
 \brief Header file for FabricationPlant class.

 The aim of the Class is to manage evolution of FabricationPlant

 
 @author BaM, Marc
 @version 2.0
 */



#include <vector>
#include <map>

#include "FabricationPlant.hxx"


using namespace std;
typedef long long int cSecond;

class PWR_THU_FabricationPlant : public FabricationPlant
{
	//on utilise les constructeur de la classe FabricationPlant
	// dans le code :MyFabricationPlant *FP_MOX = new FabricationPlant(gCLASS.GetLog(),Stock, ReUsable);
	//on a besoin de surcherger le constructeur que si
public :
	PWR_THU_FabricationPlant();
	PWR_THU_FabricationPlant(CLASSLogger* log);
	
	PWR_THU_FabricationPlant(CLASSLogger* log, Storage* storage, Storage* reusable, double fabricationtime = 365.25*24*3600*2);
	///< Normal Destructor.
	~PWR_THU_FabricationPlant();
	
	void	BuildFuelForReactor(int ReactorId);			//obligatoire, c'est un peu la base
protected :
//********* Private Method *********//
// mettre ici les méthodes dont on a besoin dans sa fonction BuildFuelForReactor perso

};

#endif