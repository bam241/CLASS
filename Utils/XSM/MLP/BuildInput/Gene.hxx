#include "ZAI.hxx"
#include <TGraph.h>
#include <TGraph.h>

#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <map>
#include <iostream>
#include <iomanip>
using namespace std;

string dtoa(double num)
{
	ostringstream os(ostringstream::out);
	os<<setprecision(3)<<num;
	return os.str();
}

string FilePath;
string DataPath;

map<ZAI, double> ZAImass;
vector<string> JobName;
vector<string> GoodJobName;

vector<double> fTime;						//Time vector of the depletion calculation (second)
vector<ZAI> fAllNuclei;						//All the nuclei present in the fuel

vector< map < ZAI, vector<double> > > fXSFis;	// map of fission cross section fXSFis[NumberOfTheEvolution][ZAI][TimeStep]
vector< map < ZAI, vector<double> > > fXSCap;
vector< map < ZAI, vector<double> > > fXSN2N;

vector<IsotopicVector> fActinideCompoInit;	//Fresh fuel composition

int fNOfTimeStep=0; //number of time step in the Evolution

string fEvolutionDataFolder = "";

bool fIsAllNucleiAlreadyFill=false;			

void	InitMass();
void	CheckJob();
void	ReadAndFill(string jobname);
void    DumpForTestingNeuron(string filename);
void    DumpInputNeuron(string filename);