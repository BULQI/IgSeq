#ifndef AFFINEGAPMODEL_HPP
#define AFFINEGAPMODEL_HPP

//#include "pairwiseAlignment.hpp"
#include "GapModel.hpp"

class AffineGapModel:public GapModel
{

public:
  AffineGapModel(const double& _gopen, const double& _gextension);

  virtual ~AffineGapModel();
  virtual double GapValue(const TracebackTable* _tbTable, const unsigned int& _patternIndex, const unsigned int& _subjectIndex, const bool& _patternGap,
			  const double& _prevEntryValue,
			  double& _MaxGapValue, unsigned int& _MaxGapIndex) const;

  double GetGapOpenValue() const;
  double GetGapExtensionValue() const;
  void SetGapOpenValue(const double& _gopen);
  void SetGapExtensionValue(const double& _gextension);

protected:
  double c_gopen;
  double c_gextension;
  
};
#endif 
