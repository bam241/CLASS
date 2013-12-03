#include "CLASSWin.hxx"



#include <TApplication.h>
#include <string>
#include <stdlib.h>

using namespace std;

int main(int argc, char** argv)
{
	if(argc<2)
	{
		cerr<<endl<<"Usage: CLASGui FileName1 FileName2 ... "<<endl;
		exit(0);
	}

	vector<string> VFileName;
	int NumberOfFile=argc-1
	;
	for(int i=0;i<NumberOfFile;i++)
		VFileName.push_back(string(argv[i+1]));

	
	CLASSRead* DataRead = new CLASSRead(VFileName[0]);
	for (int i = 1; i < (int)VFileName.size(); i++)
	{
		cout << "FIle "<<  i << endl;
		DataRead->AddFile(VFileName[i]);
	}
	DataRead->ReadName();
	DataRead->ReadZAI();

	
	argc=1; //avoid to change directory by root TApplication...
	TApplication theApp("App", &argc, argv);
	MainWin *win=new MainWin(DataRead);

	
	theApp.Run();
	return 0;
}


//g++ -O3 `root-config --cflags` -o Inventory Inventory.cxx  `root-config --glibs` 