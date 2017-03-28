// =====================================================================================

#include "Scenar_t.hxx"

ClassImp(Scenar_t) // Needed by ROOT

//    
// Const
// =====================================================================================
Scenar_t::Scenar_t():TObject()
{
    NTimeStep = 0;
    S_IsOk = false ;
    //
    UOx_CT = 0.0;
    MOx_CT = 0.0;
    // Cumulated Number of Load
    UOx_NLOAD = 0;
    MOx_NLOAD = 0;
    MOx_MLOAD = 0;
}
Scenar_t::Scenar_t(const unsigned short vlength):TObject()
{

        NTimeStep = vlength;
        Time.resize(vlength); // [Year]
        Power.resize(vlength);
        //
        S_IsOk = false ;
        //
        UOx_CT = 0.0;
        MOx_CT = 0.0;
        // Cumulated Number of Load
        UOx_NLOAD = 0;
        MOx_NLOAD = 0;
        MOx_MLOAD = 0;

        // Facilities IV of Interest
        //
        v_PWR_UOx_BOC.resize(vlength);
        v_PWR_MOx_BOC.resize(vlength);
               
        v_PWR_UOx_EOC.resize(vlength);
        v_PWR_MOx_EOC.resize(vlength);

        v_StockUOx_CIN.resize(vlength);
        v_StockUOx_COU.resize(vlength);
        v_StockUOx.resize(vlength);

        v_StockMOx_CIN.resize(vlength);
        v_StockMOx_COU.resize(vlength);
        v_StockMOx.resize(vlength);
        
        // Isotopic vectors, CYCLE
        //
        v_IVTOTAL.resize(vlength);
        v_IVINCYCLE.resize(vlength);
        v_IVWASTE.resize(vlength);
}
//
// Dest
// =====================================================================================
Scenar_t::~Scenar_t()
{}

//
// Methods
// =====================================================================================

void Scenar_t::Print()
{

    cout << " =============================================================== " << endl;
    cout << " Scenario Info " << endl;
    cout << " BU_UOx   : " << this->BU_UOx    << endl;
    cout << " TC_UOx   : " << this->TC_UOx    << endl;
    cout << " BU_MOx   : " << this->BU_MOx    << endl;
    cout << " TC_MOx   : " << this->TC_MOx    << endl;
    cout << " TF_MOx   : " << this->TF_MOx    << endl;
    cout << " IsMOxAm  : " << this->IsMOxAm   << endl;
    cout << " Fr_MOx   : " << this->Fr_MOx    << endl;
    cout << " NB_Fuel  : " << this->NB_Fuel   << endl;
    cout << " Ks_Fuel  : " << this->Ks_Fuel   << endl;
    cout << " LF_Fuel  : " << this->LF_Fuel   << endl;
    cout << " StMMOx   : " << this->StMMOx    << endl;
    cout << " NTimeStep: " << this->NTimeStep << endl;
    cout << " TimeStep : " << this->TimeStep  << endl;
    cout << " S_IsOk   : " << this->S_IsOk    << endl;
    cout << " =============================================================== " << endl;

}


void Scenar_t::Clear()
{

        Time.clear(); // [Year]
        Power.clear();

        // Facilities IV of Interest
        //
        v_PWR_UOx_BOC.clear();
        v_PWR_MOx_BOC.clear();
               
        v_PWR_UOx_EOC.clear();
        v_PWR_MOx_EOC.clear();

        v_StockUOx_CIN.clear();
        v_StockUOx_COU.clear();
        v_StockUOx.clear();

        v_StockMOx_CIN.clear();
        v_StockMOx_COU.clear();
        v_StockMOx.clear();
        
        // Isotopic vectors, CYCLE
        //
        v_IVTOTAL.clear();
        v_IVINCYCLE.clear();
        v_IVWASTE.clear();
}


// =====================================================================================

vector<double> Scenar_t::GetPWR_UOx_BOC(unsigned short Z)
{//
    vector<double> v_tBOC;
    v_tBOC.reserve(this->NTimeStep);
    for (unsigned short i = 0; i < this->NTimeStep; i++){
        v_tBOC.push_back(this->v_PWR_UOx_BOC[i]->GetSpeciesComposition(Z).GetTotalMass());
    }
    return v_tBOC;
}

vector<double> Scenar_t::GetPWR_UOx_BOC(unsigned short Z, unsigned short A, unsigned short I)
{//
    IsotopicVector tIV; tIV.Add(Z,A,I,1);
    vector<double> v_tBOC;
    for (unsigned short i = 0; i < this->NTimeStep; i++){
        v_tBOC.push_back(this->v_PWR_UOx_BOC[i]->GetThisComposition(tIV).GetTotalMass());
    }
    return v_tBOC;
}
///

vector<double> Scenar_t::GetPWR_MOx_BOC(unsigned short Z)
{//
    vector<double> v_tBOC;
    v_tBOC.reserve(this->NTimeStep);
    for (unsigned short i = 0; i < this->NTimeStep; i++){
        v_tBOC.push_back(this->v_PWR_MOx_BOC[i]->GetSpeciesComposition(Z).GetTotalMass());
    }
    return v_tBOC;
}

vector<double> Scenar_t::GetPWR_MOx_BOC(unsigned short Z, unsigned short A, unsigned short I)
{//
    IsotopicVector tIV; tIV.Add(Z,A,I,1);
    vector<double> v_tBOC;
    v_tBOC.reserve(this->NTimeStep);
    for (unsigned short i = 0; i < this->NTimeStep; i++){
        v_tBOC.push_back(this->v_PWR_MOx_BOC[i]->GetThisComposition(tIV).GetTotalMass());
    }
    return v_tBOC;
}

////

vector<double> Scenar_t::GetPWR_UOx_EOC(unsigned short Z)
{//
    vector<double> v_tEOC;
    v_tEOC.reserve(this->NTimeStep);
    for (unsigned short i = 0; i < this->NTimeStep; i++){
        v_tEOC.push_back(this->v_PWR_UOx_EOC[i]->GetSpeciesComposition(Z).GetTotalMass());
    }
    return v_tEOC;
}

vector<double> Scenar_t::GetPWR_UOx_EOC(unsigned short Z, unsigned short A, unsigned short I)
{//
    IsotopicVector tIV; tIV.Add(Z,A,I,1);
    vector<double> v_tEOC;
    v_tEOC.reserve(this->NTimeStep);
    for (unsigned short i = 0; i < this->NTimeStep; i++){
        v_tEOC.push_back(this->v_PWR_UOx_EOC[i]->GetThisComposition(tIV).GetTotalMass());
    }
    return v_tEOC;
}
///

vector<double> Scenar_t::GetPWR_MOx_EOC(unsigned short Z)
{//
    vector<double> v_tEOC;
    v_tEOC.reserve(this->NTimeStep);
    for (unsigned short i = 0; i < this->NTimeStep; i++){
        v_tEOC.push_back(this->v_PWR_MOx_EOC[i]->GetSpeciesComposition(Z).GetTotalMass());
    }
    return v_tEOC;
}

vector<double> Scenar_t::GetPWR_MOx_EOC(unsigned short Z, unsigned short A, unsigned short I)
{//
    IsotopicVector tIV; tIV.Add(Z,A,I,1);
    vector<double> v_tEOC;
    v_tEOC.reserve(this->NTimeStep);
    for (unsigned short i = 0; i < this->NTimeStep; i++){
        v_tEOC.push_back(this->v_PWR_MOx_EOC[i]->GetThisComposition(tIV).GetTotalMass());
    }
    return v_tEOC;
}


// =====================================================================================

vector<double> Scenar_t::GetTotalInvOfIsotope(unsigned short Z)
{
    vector<double> v_tTotalInv;
    v_tTotalInv.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){//
        v_tTotalInv.push_back((this->v_IVTOTAL)[i]->GetSpeciesComposition(Z).GetTotalMass());
    }
    return v_tTotalInv;
}

vector<double> Scenar_t::GetTotalInvOfIsotope(unsigned short Z, unsigned short A, unsigned short I)
{
    IsotopicVector tIV; tIV.Add(Z,A,I,1);
    vector<double> v_tTotalInv;
    v_tTotalInv.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){//
        v_tTotalInv.push_back((this->v_IVTOTAL)[i]->GetThisComposition(tIV).GetTotalMass());
    }
    return v_tTotalInv;
}

double Scenar_t::GetTotalInvOfIsotopeAtTime(unsigned short Z, double time)
{
    if (this->v_IVTOTAL.empty()) return 0.0;
    else{
        unsigned short i = (unsigned short) (time/(this->TimeStep));
        return  (this->v_IVTOTAL)[i]->GetSpeciesComposition(Z).GetTotalMass();
    }
}

double Scenar_t::GetTotalInvOfIsotopeAtTime(unsigned short Z, unsigned short A, unsigned short I, double time)
{
    if (this->v_IVTOTAL.empty()) return 0.0;
    else{
        unsigned short i = (unsigned short) (time/(this->TimeStep));
        IsotopicVector tIV; tIV.Add(Z,A,I,1);
        //
        return  (this->v_IVTOTAL)[i]->GetThisComposition(tIV).GetTotalMass();
    }
}

// =====================================================================================

vector<double> Scenar_t::GetIncycleInvOfIsotope(unsigned short Z)
{
    vector<double> v_tIncycleInv;
    v_tIncycleInv.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){//
        v_tIncycleInv.push_back((this->v_IVINCYCLE)[i]->GetSpeciesComposition(Z).GetTotalMass());
    }
    return v_tIncycleInv;
}

vector<double> Scenar_t::GetIncycleInvOfIsotope(unsigned short Z, unsigned short A, unsigned short I)
{
    IsotopicVector tIV; tIV.Add(Z,A,I,1);
    vector<double> v_tIncycleInv;
    v_tIncycleInv.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){//
        v_tIncycleInv.push_back((this->v_IVINCYCLE)[i]->GetThisComposition(tIV).GetTotalMass());
    }
    return v_tIncycleInv;
}

double Scenar_t::GetIncycleInvOfIsotopeAtTime(unsigned short Z, double time)
{
    if (this->v_IVINCYCLE.empty()) return 0.0;
    else{
        unsigned short i = (unsigned short) (time/(this->TimeStep));
        return  (this->v_IVINCYCLE)[i]->GetSpeciesComposition(Z).GetTotalMass();
    }
}

double Scenar_t::GetIncycleInvOfIsotopeAtTime(unsigned short Z, unsigned short A, unsigned short I, double time)
{
    if (this->v_IVINCYCLE.empty()) return 0.0;
    else{
        unsigned short i = (unsigned short) (time/(this->TimeStep));
        IsotopicVector tIV; tIV.Add(Z,A,I,1);
        //
        return  (this->v_IVINCYCLE)[i]->GetThisComposition(tIV).GetTotalMass();
    }
}


// =====================================================================================

vector<double> Scenar_t::GetWasteInvOfIsotope(unsigned short Z)
{
    vector<double> v_tWasteInv;
    v_tWasteInv.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){//
        v_tWasteInv.push_back((this->v_IVWASTE)[i]->GetSpeciesComposition(Z).GetTotalMass());
    }
    return v_tWasteInv;
}

vector<double> Scenar_t::GetWasteInvOfIsotope(unsigned short Z, unsigned short A, unsigned short I)
{
    IsotopicVector tIV; tIV.Add(Z,A,I,1);
    vector<double> v_tWasteInv;
    v_tWasteInv.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){//
        v_tWasteInv.push_back((this->v_IVWASTE)[i]->GetThisComposition(tIV).GetTotalMass());
    }
    return v_tWasteInv;
}

double Scenar_t::GetWasteInvOfIsotopeAtTime(unsigned short Z, double time)
{
    if (this->v_IVWASTE.empty()) return 0.0;
    else{
        unsigned short i = (unsigned short) (time/(this->TimeStep));
        return  (this->v_IVWASTE)[i]->GetSpeciesComposition(Z).GetTotalMass();
    }
}


double Scenar_t::GetWasteInvOfIsotopeAtTime(unsigned short Z, unsigned short A, unsigned short I, double time)
{
    if (this->v_IVWASTE.empty()) return 0.0;
    else{
        unsigned short i = (unsigned short) (time/(this->TimeStep));
        IsotopicVector tIV; tIV.Add(Z,A,I,1);
        //
        return  (this->v_IVWASTE)[i]->GetThisComposition(tIV).GetTotalMass();
    }
}


// =====================================================================================

vector<double> Scenar_t::GetStockUOxFlux(unsigned short Z)
{
    vector<double> v_tStockUOxFlux;
    v_tStockUOxFlux.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){
        double cumin = (this->v_StockUOx_CIN[i])->GetSpeciesComposition(Z).GetTotalMass();
        double cumout = (this->v_StockUOx_COU[i])->GetSpeciesComposition(Z).GetTotalMass();

        v_tStockUOxFlux.push_back(cumin - cumout);
    }
    return v_tStockUOxFlux;
}

vector<double> Scenar_t::GetStockUOxFlux(unsigned short Z, unsigned short A, unsigned short I)
{
    IsotopicVector tIV; tIV.Add(Z,A,I,1);
    vector<double> v_tStockUOxFlux;
    v_tStockUOxFlux.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){
        double cumin = (this->v_StockUOx_CIN[i])->GetThisComposition(tIV).GetTotalMass();
        double cumout = (this->v_StockUOx_COU[i])->GetThisComposition(tIV).GetTotalMass();
        v_tStockUOxFlux.push_back(cumin - cumout);
    }
    return v_tStockUOxFlux;
}

vector<double> Scenar_t::GetStockUOxInvOfIsotope(unsigned short Z)
{
    vector<double> v_tStockInv;
    v_tStockInv.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){//
        v_tStockInv.push_back((this->v_StockUOx)[i]->GetSpeciesComposition(Z).GetTotalMass());
    }
    return v_tStockInv;
}


vector<double> Scenar_t::GetStockUOxInvOfIsotope(unsigned short Z, unsigned short A, unsigned short I)
{
    IsotopicVector tIV; tIV.Add(Z,A,I,1);
    vector<double> v_tStockInv;
    v_tStockInv.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){//

        v_tStockInv.push_back((this->v_StockUOx)[i]->GetThisComposition(tIV).GetTotalMass());

    }
    return v_tStockInv;
}


double Scenar_t::GetStockUOxInvOfIsotopeAtTime(unsigned short Z, double time)
{
    if (this->v_StockUOx.empty()) return 0.0;
    else{
        unsigned short i = (unsigned short) (time/(this->TimeStep));
        //
        return (this->v_StockUOx)[i]->GetSpeciesComposition(Z).GetTotalMass();
    }
}


double Scenar_t::GetStockUOxInvOfIsotopeAtTime(unsigned short Z, unsigned short A, unsigned short I, double time)
{
    if (this->v_StockUOx.empty()) return 0.0;
    else{
        unsigned short i = (unsigned short) (time/(this->TimeStep));
        IsotopicVector tIV; tIV.Add(Z,A,I,1);
        //
        return (this->v_StockUOx)[i]->GetThisComposition(tIV).GetTotalMass();
    }
}

////

vector<double> Scenar_t::GetStockMOxFlux(unsigned short Z)
{
    vector<double> v_tStockMOxFlux;
    v_tStockMOxFlux.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){
        double cumin = (this->v_StockMOx_CIN[i])->GetSpeciesComposition(Z).GetTotalMass();
        double cumout = (this->v_StockMOx_COU[i])->GetSpeciesComposition(Z).GetTotalMass();

        v_tStockMOxFlux.push_back(cumin - cumout);
    }
    return v_tStockMOxFlux;
}

vector<double> Scenar_t::GetStockMOxFlux(unsigned short Z, unsigned short A, unsigned short I)
{
    IsotopicVector tIV; tIV.Add(Z,A,I,1);
    vector<double> v_tStockMOxFlux;
    v_tStockMOxFlux.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){
        double cumin = (this->v_StockMOx_CIN[i])->GetThisComposition(tIV).GetTotalMass();
        double cumout = (this->v_StockMOx_COU[i])->GetThisComposition(tIV).GetTotalMass();
        v_tStockMOxFlux.push_back(cumin - cumout);
    }
    return v_tStockMOxFlux;
}

vector<double> Scenar_t::GetStockMOxInvOfIsotope(unsigned short Z)
{
    vector<double> v_tStockInv;
    v_tStockInv.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){//
        v_tStockInv.push_back((this->v_StockMOx)[i]->GetSpeciesComposition(Z).GetTotalMass());
    }
    return v_tStockInv;
}


vector<double> Scenar_t::GetStockMOxInvOfIsotope(unsigned short Z, unsigned short A, unsigned short I)
{
    IsotopicVector tIV; tIV.Add(Z,A,I,1);
    vector<double> v_tStockInv;
    v_tStockInv.reserve(this->NTimeStep);
    //
    for (unsigned short i =0; i < this->NTimeStep; i++){//

        v_tStockInv.push_back((this->v_StockMOx)[i]->GetThisComposition(tIV).GetTotalMass());

    }
    return v_tStockInv;
}


double Scenar_t::GetStockMOxInvOfIsotopeAtTime(unsigned short Z, double time)
{
    if (this->v_StockMOx.empty()) return 0.0;
    else{
        unsigned short i = (unsigned short) (time/(this->TimeStep));
        //
        return (this->v_StockMOx)[i]->GetSpeciesComposition(Z).GetTotalMass();
    }
}


double Scenar_t::GetStockMOxInvOfIsotopeAtTime(unsigned short Z, unsigned short A, unsigned short I, double time)
{
    if (this->v_StockMOx.empty()) return 0.0;
    else{
        unsigned short i = (unsigned short) (time/(this->TimeStep));
        IsotopicVector tIV; tIV.Add(Z,A,I,1);
        //
        return (this->v_StockMOx)[i]->GetThisComposition(tIV).GetTotalMass();
    }
}

