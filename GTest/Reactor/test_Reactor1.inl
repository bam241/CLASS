#include "CLASSHeaders.hxx"
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>
#include "Reactor.hxx"
#include "XSM_MLP.hxx"
#include "IM_RK4.hxx"
#include "EQ_OneParameter.hxx"

TEST ( TestReactor, getBurnUp ) {

    double HMMass     = 1;
    double Power      = 1;
    double LoadFactor = 1;
    double CycleTime  = 1;

    Reactor* R1 = new Reactor();
    R1->SetCycleTime(1);
    double ExpectedBurnUp = R1->GetBurnUp();
    double RealBurnUp     = Power * CycleTime / 3600. / 24. / HMMass;

    ASSERT_NEAR(RealBurnUp, ExpectedBurnUp, RealBurnUp*1e-15);
}

TEST ( TestReactor, getCycleTime ) {

    double HMMass     = 1;
    double Power      = 1;
    double LoadFactor = 1;
    double BurnUp     = 1;

    Reactor* R1 = new Reactor();
    R1->SetBurnUp(BurnUp);
    double ExpectedCycleTime = R1->GetCycleTime();
    double RealCycleTime     = (cSecond) (BurnUp*1e9 / (Power) * HMMass  *3600*24);

    ASSERT_NEAR(RealCycleTime, ExpectedCycleTime, RealCycleTime*1e-15);
}
