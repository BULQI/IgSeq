#include "AlignmentString.hpp"
//#include <string>
#include <sstream>
using namespace std;

AlignmentString::AlignmentString():c_pattern(""), c_subject(""),c_pattern_start(-1), c_pattern_end(-1),
   c_subject_start(-1), c_subject_end(-1),
				   c_pattern_w_gap( ""), c_subject_w_gap( ""), c_score(0)
{
  //empty
}
AlignmentString::AlignmentString(const string& _pattern, const int& p_start, const int& p_end,
		  const string& _subject, const int& c_start, const int& c_end,
				 const string& _pattern_wg, const string& _subject_wg, const double& _score):c_pattern(_pattern), 
   c_subject( _subject), c_pattern_start( p_start), c_pattern_end( p_end),
   c_subject_start(c_start), c_subject_end( c_end),
													     c_pattern_w_gap( _pattern_wg), c_subject_w_gap( _subject_wg),c_score(_score)
{
  //empty
}
AlignmentString::AlignmentString(const AlignmentString& a):c_pattern(a.c_pattern), 
   c_subject( a.c_subject), c_pattern_start( a.c_pattern_start), c_pattern_end( a.c_pattern_end),
   c_subject_start(a.c_subject_start), c_subject_end( a.c_subject_end),
													     c_pattern_w_gap( a.c_pattern_w_gap), c_subject_w_gap( a.c_subject_w_gap),c_score(a.c_score)
{
	//empty
}
  
  //true to set the sub string with gap (true) or without gap(fasle)
void AlignmentString::SetPattern(const string& _pattern, const bool& gapFlag)
{
  if(gapFlag)
    c_pattern_w_gap=_pattern;
  else
    c_pattern=_pattern;

}
void AlignmentString::SetSubject(const string& _subject, const bool& gapFlag)
{
  if(gapFlag)
    c_subject_w_gap=_subject;
  else
    c_subject=_subject;
}
  
string AlignmentString::GetPattern(const bool& gapFlag) const
{
  if(gapFlag)
    return c_pattern_w_gap;
  else
    return c_pattern;
}
string AlignmentString::GetSubject(const bool& gapFlag) const
{
  if(gapFlag)
    return c_subject_w_gap;
  else
    return c_subject;
}
  
void AlignmentString::SetPatternIndex(const int& _start, const int& _end)
{
  c_pattern_start=_start;
  c_pattern_end=_end;
}
void AlignmentString::SetSubjectIndex(const int& _start, const int& _end)
{
  c_subject_start=_start;
  c_subject_end=_end;
}

unsigned int AlignmentString::GetPatternIndexStart() const
{
  return c_pattern_start;
}
unsigned int AlignmentString::GetPatternIndexEnd() const
{
  return c_pattern_end;
}
unsigned int AlignmentString::GetSubjectIndexStart() const
{
  return c_subject_start;
}
unsigned int AlignmentString::GetSubjectIndexEnd() const
{
  return c_subject_end;
}

string AlignmentString::toString() const
{
  ostringstream oss("");
  ostringstream intoss("");
  intoss<<"["<<c_pattern_start<<"]";
  string tempIntStr=intoss.str();

  ostringstream intossSub("");
  intossSub<<"["<<c_subject_start<<"]";
  string tempIntStrSub=intossSub.str();  

  //check to see the length and prefix with space if necessary
  while(tempIntStr.length()>tempIntStrSub.length())
    {
      tempIntStrSub=" "+tempIntStrSub;
    }
  while(tempIntStr.length()<tempIntStrSub.length())
    {
      tempIntStr=" "+tempIntStr;
    }
  oss<<"Alignment Strings \nPattern:"<<tempIntStr<<c_pattern_w_gap<<"\n";
  oss<<"Subject:"<<tempIntStrSub<<c_subject_w_gap<<"\n";
  oss<<"Score:"<<c_score<<"\n";

  return oss.str();
}


void AlignmentString::SetScore(const double& _s)
{
  c_score=_s;
}
double AlignmentString::GetScore() const
{
  return c_score;
}
