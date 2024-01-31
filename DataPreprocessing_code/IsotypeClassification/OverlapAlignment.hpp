#ifndef OVERLAPALIGNMENT_HPP
#define OVERLAPALIGNMENT_HPP
#include "pairwiseAlignment.hpp"
//#include <vector>
using namespace std;

//overlap alignment following R style, doing global alignment with without gap penalty on both ends
class OverlapAlignment: public PairwiseAlignment
{
public:
  OverlapAlignment(SequenceString* _pattern, SequenceString* _subject, 
		 ScoreMatrix* _m=&nuc44, const double& _gopen=-8, 
		 const double& _gextension=-5, const double& _scale=1,
		 const short& _typeOfGapModel=1);
  
  virtual ~OverlapAlignment();

  //const double* GetScoreArr();
  //const AlignmentString* GetAlignmentArr();
  //void alignLM();
protected:
  virtual void align();
  virtual void traceBack();
};


#endif
