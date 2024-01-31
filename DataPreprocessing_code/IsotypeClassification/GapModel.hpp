#ifndef GAPMODEL_HPP
#define GAPMODEL_HPP

//#include "pairwiseAlignment.hpp"
#include "TracebackTable.hpp"

//this is the abstract class base class for gap models

class GapModel
{
public:
  GapModel();
  virtual ~GapModel()=0;

  //input: _patternGap, bool, true for the gaps on the pattern string; false otherwise
  virtual double GapValue(const TracebackTable* _tbTable, const unsigned int& _patternIndex, const unsigned int& _subjectIndex, const bool& _patternGap, 
			  const double& _prevEntryValue,
			  double& _MaxGapValue, unsigned int& _MaxGapIndex) const =0;
  
};


#endif
