#include "CLASSReader.hxx"


//____________________________________________________________________________
CLASSReader::CLASSReader () :
	freader( new TMVA::Reader("silent") )
{ ; }

//____________________________________________________________________________
CLASSReader::CLASSReader ( const std::vector<std::string> & names ) :
	freader( new TMVA::Reader("silent") ) , finputTMVA( names.size() )
{
	std::list<float>::iterator l_it = finputTMVA.begin();
	for ( std::vector<std::string>::const_iterator v_it = names.begin() ; v_it != names.end() ; ++v_it, ++l_it )
	{
		freader->AddVariable( v_it->c_str() , &(*l_it) );
	}
}

//____________________________________________________________________________
CLASSReader::~CLASSReader ()
{
	delete freader;
}

//____________________________________________________________________________
void CLASSReader::AddVariable ( const std::string & name )
{
	finputTMVA.push_back( 0 );
	freader->AddVariable( name.c_str() , &finputTMVA.back() );
	
	const std::vector<TString> names = freader->DataInfo().GetListOfVariables();
}

//____________________________________________________________________________
void CLASSReader::SetInputData ( vector<float> input )
{
	const std::vector<TString> names = freader->DataInfo().GetListOfVariables();
	const std::size_t N = names.size();
	
	std::list<float>::iterator l_it = finputTMVA.begin();
	for ( std::size_t i = 0 ; i!=N ; ++i, ++l_it )
		{ *l_it = input[i]; }
}


//____________________________________________________________________________
TMVA::IMethod * CLASSReader::BookMVA ( const std::string & methodTag , const std::string & weightfile )
{
	return freader->BookMVA( methodTag.c_str() , weightfile.c_str() );
}

//____________________________________________________________________________
const std::vector<float> & CLASSReader::EvaluateRegression ( const std::string & methodTag )
{
	return freader->EvaluateRegression( methodTag.c_str() );
}

