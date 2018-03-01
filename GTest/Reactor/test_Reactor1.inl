#include "CLASSHeaders.hxx"
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>
#include "Reactor.hxx"
#include "XSM_MLP.hxx"
#include "IM_RK4.hxx"
#include "EQ_OneParameter.hxx"

TEST ( Reactor, GetBurnUp ) {

    double HMMass     = 1.;
    double Power      = 1.;
    double LoadFactor = 1.;
    double CycleTime  = 1.;

    Reactor* R1 = new Reactor();
    R1->SetHMMass(HMMass); // Mass should be defined first
    R1->SetPower(LoadFactor * Power);
    R1->SetCycleTime(CycleTime);

    double ExpectedBurnUp = R1->GetBurnUp();
    double RealBurnUp     = Power * LoadFactor * CycleTime / 3600. / 24. / HMMass;

    ASSERT_NEAR(RealBurnUp, ExpectedBurnUp, RealBurnUp*1e-15);
}

TEST ( Reactor, GetCycleTime ) {

    double HMMass     = 1.;
    double Power      = 1.;
    double LoadFactor = 1.;
    double BurnUp     = 1.;

    Reactor* R1 = new Reactor();
    R1->SetHMMass(HMMass); // Mass should be defined first
    R1->SetBurnUp(BurnUp);
    R1->SetPower(LoadFactor * Power);

    cSecond ExpectedCycleTime = R1->GetCycleTime();
    cSecond RealCycleTime     = (cSecond) (BurnUp*1e9 * HMMass  *3600*24 / (Power * LoadFactor));

    ASSERT_NEAR(RealCycleTime, ExpectedCycleTime, 0.5);
}
