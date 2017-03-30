#include <string>
#include <iostream>
#include <IsotopicVector.hxx>

using namespace std;

class Facility 
{

    public : 

    // Constructors
    Facility();
    // estructor
    ~Facility();

    string GetName() {return fName;}
    void SetName(string Name) {fName = Name;}

    private : 

    string fName;
    IsotopicVector fIV; // The IV vector according to the time



   };
