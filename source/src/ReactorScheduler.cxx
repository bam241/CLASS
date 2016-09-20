#include "ReactorScheduler.hxx"

#include "CLASSLogger.hxx"

using namespace std;


//________________________________________________________________________
//
//		ReactorModel
//
//________________________________________________________________________

ReactorModel::ReactorModel(EvolutionData* evo)
{
	fEvolutionData = evo;
	fPhysicsModels = 0;
}

//________________________________________________________________________
ReactorModel::ReactorModel(PhysicsModels* evo)
{
	fEvolutionData = 0;
	fPhysicsModels = evo;
}

//________________________________________________________________________
//
//		ScheduleEntry
//
//________________________________________________________________________

ScheduleEntry::ScheduleEntry(ReactorModel* ED_Or_PhysMod, double BurnUp, double Power, double HMMass )
{
	DBGL

    fReactorModel   = ED_Or_PhysMod ;
    fBurnUp         = BurnUp ;
    fPower          = Power ;
    fHeavyMetalMass = HMMass ;

	DBGL
}

//________________________________________________________________________
//
//		ReactorScheduler
//
//________________________________________________________________________

//________________________________________________________________________
ReactorScheduler::ReactorScheduler():CLASSObject()
{
	DBGL
}

//________________________________________________________________________
ReactorScheduler::ReactorScheduler(CLASSLogger* log):CLASSObject(log)
{
	DBGL
}

//________________________________________________________________________
void ReactorScheduler::AddEntry(cSecond time,  ReactorModel* Model, double BurnUp, double Power, double HMMass)
{

	DBGV("time : "<< time<<"  BU : " << BurnUp<<" Power : "<<Power <<" HMMass " << HMMass <<endl);
	fReactorSchedulerMap.insert( pair<cSecond,ScheduleEntry*>( time, new ScheduleEntry( Model, BurnUp, Power, HMMass )));
	DBGL
}

//________________________________________________________________________
ScheduleEntry* ReactorScheduler::GetEntryAt(cSecond t)
{
/*	map<cSecond , ScheduleEntry*>::iterator it;
	for(it=fReactorSchedulerMap.begin();it!=fReactorSchedulerMap.end();it++){
		DBGV("time : "<<it->first <<"  BU : "<<it->second->GetBurnUp()
									<<" Power : "<<it->second->GetPower()
									<<" HMMass : "<<it->second->GetHeavyMetalMass()<<endl
			);
	}
*/
	DBGV("Getting ScheduleEntry at time : "<<t<<" s"<<endl);
	map<cSecond , ScheduleEntry*>::iterator itlow;
	itlow = fReactorSchedulerMap.lower_bound(t);
	DBGL

    if(itlow->first == t){
    	DBGL
        return itlow->second;
    }
    else if(itlow != fReactorSchedulerMap.begin()){
    	DBGL
        return (--itlow)->second;
    }
    else {
		DBGL
        return fReactorSchedulerMap.begin()->second;

    }

}








