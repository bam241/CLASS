#include <string>
#include <iostream>

using namespace std;



class Run 
{

    public :
  
    // Contructors
    Run();
    // Destructors
    ~Run();

    // Get and Set methods
    map<string, double> GetRunParameter() {return fRunParameter;} // The map of user input parameters
    void SetRunParameter(map<string, double> RunParameter) {fRunParameter=RunParameter;} // The map of user input parameters

    int GetRunID() {return fRunID;} // Integer ID for this run
    void SetRunID(int RunID) {fRunID=RunID;} // Integer ID for this run

    bool GetRunIsOK() {return fIsOK;} // Tell is Run is OK
    void SetRunIsOK(bool IsOK) {fIsOK = IsOK;} // Tell is Run is OK

    int GetNTimeStep() {return fNTimeStep;} // Number of Time Steps
    void SetNTimeStep(int NTimeStep) {fNTimeStep=NTimeStep;} // Number of Time Steps

    vector<double> GetTimeStepVector {return fTimeStepVector;} // Vector of time steps
    void SetTimeStepVector(vector<double> TimeStepVector) {fTimeStepVector = TimeStepVector;} // Vector of time steps

    private : 

    map<string, double> fRunParameter; // The map of user input parameters

    int fRunID; // Integer ID for this run

    bool fIsOK; // Tell is Run is OK

    int fNTimeStep; // Number of Time Steps
    vector<double> fTimeStepVector; // Vector of time steps

    // Number of Facilities
    int fNumberOfFacilities;
    int fNumberOfRectaors;
    int fNumberOfStorages;
    int fNumberOfPools;
    int fNumberOfFabPlants;
    int fNumberOfSepPlants;

    // Names of Facilities
    vector<string> fNamesOfReactors;
    vector<string> fNamesOfStorages;
    vector<string> fNamesOfPools;
    vector<string> fNamesOfFabPlants;
    vector<string> fNamesOfSepPlants;
};
