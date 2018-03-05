#ifndef _READER_HXX
#define _READER_HXX

#include <iostream>
#include <list>
#include <map>
#include <vector>

#include "TMVA/MethodCuts.h"
#include "TMVA/Reader.h"
#include "TMVA/Tools.h"
#include "TString.h"
#include "TTree.h"

using namespace std;

class CLASSReader {
 public:
  CLASSReader();
  CLASSReader(const std::vector<std::string>&);
  template <typename T>
  CLASSReader(const std::map<T, std::string>&);

  ~CLASSReader();

  void AddVariable(const std::string&);
  void SetInputData(vector<float> inputdata);

  TMVA::IMethod* BookMVA(const std::string&, const std::string&);
  const std::vector<float>& EvaluateRegression(const std::string&);

  std::size_t GetNVariables() const {
    return freader->DataInfo().GetNVariables();
  }
  std::vector<TString> GetListOfVariables() const {
    return freader->DataInfo().GetListOfVariables();
  }

 private:
  TMVA::Reader* freader;
  std::list<float> finputTMVA;  // changer vector en list !!!
};

template <typename T>
CLASSReader::CLASSReader(const std::map<T, std::string>& a)
    : freader(new TMVA::Reader("silent")), finputTMVA(a.size()) {
  std::list<float>::iterator l_it = finputTMVA.begin();
  for (typename std::map<T, std::string>::const_iterator m_it = a.begin();
       m_it != a.end(); ++m_it, ++l_it) {
    freader->AddVariable(m_it->second.c_str(), &(*l_it));
  }
}

#endif
