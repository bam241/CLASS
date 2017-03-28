#ifndef SCENAR_T_HXX
#define SCENAR_T_HXX

#include <vector>

#include "ZAI.hxx"
#include "IsotopicVector.hxx"
#include "TObject.h"

//
// data structure
//

// Length of the vector for the constuctor
// this is actually equivalent to NTimeStep

class ZAI;
class Isotopicvector;
//
class Scenar_t: public TObject
{
    public:
        //
        Scenar_t();
        Scenar_t(const unsigned short vlength);
        virtual ~Scenar_t();
        //
        double BU_UOx; // [GWd/t]
        double TC_UOx;
        double BU_MOx;
        double TC_MOx;
        double TF_MOx;
        bool  IsMOxAm; //
        double Fr_MOx; // [%]
        double NB_Fuel;
        double Ks_Fuel;
        double LF_Fuel;
        unsigned short StMMOx;
        unsigned short NTimeStep;
        double TimeStep;
        vector<double> Time; // [Year]
        vector<double> Power;
        //
        bool S_IsOk;
        //
        double UOx_CT;
        double MOx_CT;
        // Cumulated Number of Load
        unsigned short UOx_NLOAD;
        unsigned short MOx_NLOAD;
        unsigned short MOx_MLOAD;

        // Facilities IV of Interest
        //
        vector<IsotopicVector*> v_PWR_UOx_BOC;
        vector<IsotopicVector*> v_PWR_MOx_BOC;
               
        vector<IsotopicVector*> v_PWR_UOx_EOC;
        vector<IsotopicVector*> v_PWR_MOx_EOC;

        vector<IsotopicVector*> v_StockUOx_CIN;
        vector<IsotopicVector*> v_StockUOx_COU;
        vector<IsotopicVector*> v_StockUOx;

        vector<IsotopicVector*> v_StockMOx_CIN;
        vector<IsotopicVector*> v_StockMOx_COU;
        vector<IsotopicVector*> v_StockMOx;
        
        // Isotopic vectors, CYCLE
        //
        vector<IsotopicVector*> v_IVTOTAL;
        vector<IsotopicVector*> v_IVINCYCLE;
        vector<IsotopicVector*> v_IVWASTE;

        // 
        // Methods
        //

        void Print();
        void Clear();
        //
        vector<double>  GetPWR_UOx_BOC(unsigned short Z); // 
        vector<double>  GetPWR_UOx_BOC(unsigned short Z, unsigned short A, unsigned short I=0);
        vector<double>  GetPWR_MOx_BOC(unsigned short Z);
        vector<double>  GetPWR_MOx_BOC(unsigned short Z, unsigned short A, unsigned short I=0);

        vector<double>  GetPWR_UOx_EOC(unsigned short Z); //
        vector<double>  GetPWR_UOx_EOC(unsigned short Z, unsigned short A, unsigned short I=0); //
        vector<double>  GetPWR_MOx_EOC(unsigned short Z);
        vector<double>  GetPWR_MOx_EOC(unsigned short Z, unsigned short A, unsigned short I=0);

        //
        vector<double>  GetTotalInvOfIsotope(unsigned short Z);
        vector<double>  GetTotalInvOfIsotope(unsigned short Z, unsigned short A, unsigned short I=0);
        double    GetTotalInvOfIsotopeAtTime(unsigned short Z, double Time);
        double    GetTotalInvOfIsotopeAtTime(unsigned short Z, unsigned short A, unsigned short I, double time);

        vector<double>  GetIncycleInvOfIsotope(unsigned short Z);
        vector<double>  GetIncycleInvOfIsotope(unsigned short Z, unsigned short A, unsigned short I=0);
        double    GetIncycleInvOfIsotopeAtTime(unsigned short Z, double Time);
        double    GetIncycleInvOfIsotopeAtTime(unsigned short Z, unsigned short A, unsigned short I, double time);

        vector<double>  GetWasteInvOfIsotope(unsigned short Z);
        vector<double>  GetWasteInvOfIsotope(unsigned short Z, unsigned short A, unsigned short I=0);
        double    GetWasteInvOfIsotopeAtTime(unsigned short Z, double Time);
        double    GetWasteInvOfIsotopeAtTime(unsigned short Z, unsigned short A, unsigned short I, double time);

        vector<double>  GetStockUOxFlux(unsigned short Z);
        vector<double>  GetStockUOxFlux(unsigned short Z, unsigned short A, unsigned short I=0);
        vector<double>  GetStockUOxInvOfIsotope(unsigned short Z);
        vector<double>  GetStockUOxInvOfIsotope(unsigned short Z, unsigned short A, unsigned short I=0);
        double  GetStockUOxInvOfIsotopeAtTime(unsigned short Z, double time);
        double  GetStockUOxInvOfIsotopeAtTime(unsigned short Z, unsigned short A, unsigned short I, double time);

        vector<double>  GetStockMOxFlux(unsigned short Z);
        vector<double>  GetStockMOxFlux(unsigned short Z, unsigned short A, unsigned short I=0);
        vector<double>  GetStockMOxInvOfIsotope(unsigned short Z);
        vector<double>  GetStockMOxInvOfIsotope(unsigned short Z, unsigned short A, unsigned short I=0);
        double  GetStockMOxInvOfIsotopeAtTime(unsigned short Z, double time);
        double  GetStockMOxInvOfIsotopeAtTime(unsigned short Z, unsigned short A, unsigned short I, double time);

        //
        ClassDef (Scenar_t,1)

};

#endif
