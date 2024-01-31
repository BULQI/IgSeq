#ifndef MARKOVCHAINGAPMODEL_HPP
#define MARKOVCHAINGAPMODEL_HPP

#include <cmath>
#include "AffineGapModel.hpp"

class MarkovChainGapModel_454:public AffineGapModel
{
public:
  MarkovChainGapModel_454(const double& _gopen, const double& _gextension,const string& _patternStr, const string& _subjectStr);
  virtual ~MarkovChainGapModel_454();
  virtual double GapValue(const TracebackTable* _tbTable, const unsigned int& _patternIndex, const unsigned int& _subjectIndex, const bool& _patternGap,
			  const double& _prevEntryValue,
			  double& _MaxGapValue, unsigned int& _MaxGapIndex) const;
protected:
  string c_patternString;
  string c_subjectString;

  double GetMarkovChainGapValue(const double& _regularGapOpen, const double& _regularGapExtension, const unsigned int& _runLen, const unsigned int& _gapLen, const bool& _extensionFlag) const; 
};
#endif
