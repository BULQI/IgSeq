#ifndef PAIRWISEALIGNMENT_HPP
#define PAIRWISEALIGNMENT_HPP

#include "score.hpp"
#include "Accessory/SequenceString.hpp"
#include "AlignmentString.hpp"
#include "AffineGapModel.hpp"
#include "MarkovChainGapModel_454.hpp"
#include "TracebackTable.hpp"

//this file is used to define a interface/abstract class for pairwise alignment classes
//LocalAlign
//overlapAlign
//globalAlign will be implemented in the future


//Abstract class
class PairwiseAlignment
{
  //public
public:
  //pairwiseAlignment();
  //shit!! it seems we now keep a pointer to an outside sequenstirng for pattern and subject. so we have to be
  //really carefully for the scope and don't let the pattern and subject to get out of the scope earlier.!!!
  PairwiseAlignment(const SequenceString* _pattern, const SequenceString* _subject, 
		    const ScoreMatrix* _m=&nuc44, const double& _gopen=-8, 
		    const double& _gextension=-5, const double& _scale=1,const short& _typeOfGapModel=1);
  
  virtual ~PairwiseAlignment()=0;

  double GetScore();
  AlignmentString GetAlignment();

protected:
  PairwiseAlignment(){};
  virtual void align()=0;
  virtual void traceBack();

  const SequenceString* c_pattern;//this is follwing R style pairwiseAlignment
  const SequenceString* c_subject;//this is following R style pairwiseAlignment
  const ScoreMatrix* c_sm;
  double c_gapOpen;
  double c_gapExtension;
  double c_scale;
  
  AlignmentString c_alignment;
  double c_score;
  unsigned int c_optimalIndex[2];
  
  //double* c_dp_table;
  TracebackTable* c_traceback_table;
  GapModel* c_gm;
};

#endif
