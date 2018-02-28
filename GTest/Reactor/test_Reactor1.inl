#include "CLASSHeaders.hxx"
#include <sstream>
#include <iomanip>
#include <math.h>
#include <string>
#include "Reactor.hxx"
#include "XSM_MLP.hxx"
#include "IM_RK4.hxx"
#include "EQ_OneParameter.hxx"

TEST ( TestReactor1, getSize ) {

    XSM_MLP*            XSMUOX      = new XSM_MLP(gCLASS->GetLog(), PATH_TO_DATA + "PWR/UOX/XSModel/30Wg_FullUOX");
    EQ_OneParameter*    EQMPWRUOX   = new EQ_OneParameter(gCLASS->GetLog(), PATH_TO_DATA + "PWR/UOX/EQModel/XML/PWR_UOX.xml", PATH_TO_DATA + "PWR/UOX/EQModel/NFO/PWR_UOX.nfo");
    EQMPWRUOX->SetModelParameter("kThreshold", 1.034);
    EQMPWRUOX->SetModelParameter("NumberOfBatch", 3);
    PhysicsModels*      PMUOX       = new PhysicsModels(XSMUOX, EQMPWRUOX, IMRK4);

    double HMMass        = 1;
    double Power         = 1;
    double BurnUp        = 1;
    double LoadFactor    = 1;

    cSecond StartingTime =  0;
    cSecond LifeTime     =  100;

    Reactor* PWR_UOX = new Reactor(gCLASS->GetLog(), PMUOX, FP_UOX, PoolUOx, StartingTime, LifeTime, Power, HMMass, BurnUp, LoadFactor);

    PWR_UOX->SetName("PWR_UOx");

    EXPECT_EQ(  ,  );
}
