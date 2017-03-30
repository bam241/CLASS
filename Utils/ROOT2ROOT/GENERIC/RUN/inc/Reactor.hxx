#include <string>
#include <iostream>

using namespace std;

class Reactor : public Facility
{

    public : 

    // Constructor
    Reactor()
    // Destructor
    ~Reactor()

    // Get and Set Methods

    string GetFuelType() {return fFuelType;}
    void SetFuelType(string FuelType) {fFuelType=FuelType;}

    vector <double> GetPower() {return fPower;} // Power as a function of the time
    void SetPower(vector<double> Power) {fPower=Power;} // Power as a function of the time

    vector <double> GetEnergy() {return fEnergy;} // Energy as a function of the time
    void Setenergy(vector<double> energy) {fEnergy=Energy;} // Energy as a function of the time

    private :

    string fFuelType;

    vector <double> fPower; // Power as a function of the time
    vector <double> fEnergy; // Energy as a function of the time

    vector <double> fBurnup;  // BU as a function of the time

    double fStartingTime;
    double fShutDownTime;

    double fCycleTime;
    double fThermalPower;
    double fLoadFactor;
    double fHMMass;

    vector <IsotopicVector> fIVAtBOC; // IV @ BOC as a function of the time
};
