#include "AffineGapModel.hpp"
#include <iostream>

using namespace std;

AffineGapModel::AffineGapModel(const double& _gopen, const double& _gextension):
  c_gopen(_gopen), c_gextension(_gextension)
{
  if(c_gopen>0)
    c_gopen=-1*c_gopen;
  if(c_gextension>0)
    c_gextension=-1*c_gextension;
  //empty
}

AffineGapModel::~AffineGapModel()
{
  //empty
}

//input
//_tbTable--the tracebackTable so far.
//_patternIndex, _colIndex -- this is the "coordinate" of current entry,
//_rowGap -- this is the flag to indicate whether we are doing row gap (true) or col gap (false).
//_prevEntryValue -- this is the score of the previous entry in the same row (for row gap) or in the same col (for col gap)
//_MaxGapValue -- the running Max Gap value
//_MaxGapValue-- the index of the Max Gap starting point
//
//Output
//_MaxGapValue
//_MaxGapIndex
double AffineGapModel::GapValue(const TracebackTable* _tbTable, const unsigned int& _patternIndex, 
			const unsigned int& _subjectIndex, const bool& _patternGap,
			const double& _prevEntryValue,
			  double& _MaxGapValue, unsigned int& _MaxGapIndex) const
{
  //are we doing it by row of by or by col
  //for affine gap value, we have to know this, we only need to 
  //get the maxiGapvalue and compare with the newly opened one to get 
  //the best ones
  
  //newly opened gap
  double newGapOpenValue=_prevEntryValue+c_gopen+c_gextension;
  double extendedGapValue=_MaxGapValue+c_gextension;

  if(newGapOpenValue>=extendedGapValue)
    {
      _MaxGapValue=newGapOpenValue;
      if(_patternGap)//we are doing row gaps/pattern gaps
		{
		  //_colIndex--;
		  _MaxGapIndex=_subjectIndex-1;
		  
		}
      else  //row gaps/subject gaps
		{
		  //_rowIndex--;
		  _MaxGapIndex=_patternIndex-1;
		}
      
    }
  else
    {
      _MaxGapValue=extendedGapValue;
    }
  //_tbTable[_rowIndex+_colIndex*
  return _MaxGapValue;
}

double AffineGapModel::GetGapOpenValue() const
{
  return c_gopen;
}
double AffineGapModel::GetGapExtensionValue() const
{
  return c_gextension;
}
void AffineGapModel::SetGapOpenValue(const double& _gopen)
{
  c_gopen=_gopen;
  if(c_gopen>0)
    c_gopen=-1*c_gopen;
}
void AffineGapModel::SetGapExtensionValue(const double& _gextension)
{
  c_gextension=_gextension;
  if(c_gextension>0)
    c_gextension=-1*c_gextension;
}




