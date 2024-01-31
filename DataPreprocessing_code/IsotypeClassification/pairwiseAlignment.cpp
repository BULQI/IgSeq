#include "pairwiseAlignment.hpp"
#include <iostream>
using namespace std;


//pairwiseAlignment();
PairwiseAlignment::PairwiseAlignment(const SequenceString* _pattern, const SequenceString* _subject, 
		    const ScoreMatrix* _m, const double& _gopen, 
				     const double& _gextension, const double& _scale, const short& _typeOfGapModel):
  c_pattern( _pattern), c_subject (_subject), c_sm( _m), c_gapOpen(_gopen), 
  c_gapExtension(_gextension), c_scale( _scale), 
  c_traceback_table(NULL), c_gm(NULL)
{
  //need to allocate the 
  //this->align();
  c_optimalIndex[0]=0;
  c_optimalIndex[1]=0;
  
  switch(_typeOfGapModel)
    {
    case 1: //454 markov chain model
      c_gm=new MarkovChainGapModel_454(_gopen, _gextension, _pattern->GetSequence(), _subject->GetSequence());
      break;
      
    case 0:
    default:
      c_gm=new AffineGapModel(_gopen, _gextension);
      break;
    }
}
  
PairwiseAlignment::~PairwiseAlignment()
{
  c_pattern=NULL;
  c_subject=NULL;
  c_sm=NULL;
  
  if(c_traceback_table!=NULL)
    delete c_traceback_table;

  if(c_gm!=NULL)
    delete c_gm;
  
}

double PairwiseAlignment::GetScore()
{
  return c_score;
}
AlignmentString PairwiseAlignment::GetAlignment()
{
  return c_alignment;
}

void PairwiseAlignment::traceBack()
{
  //emptyhere need to fill
  cout<<"doing traceback in pairwiseAlignment"<<endl;
}

