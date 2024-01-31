#include "TracebackTable.hpp"
#include <iostream>

using namespace std;

TracebackTableEntry::TracebackTableEntry():c_link(ZERO), c_numOfIndels(0),c_numOfIndels_s(0), c_used(false)
{
  //empty
}

TracebackTableEntry::TracebackTableEntry(const LinkBack& _link, 
			const unsigned int& _numOfIndels, const unsigned int& _numOfIndels_s ):
  c_link(_link), c_numOfIndels(_numOfIndels), c_numOfIndels_s(_numOfIndels_s),c_used(false)
{
  //empty
} 
  
LinkBack TracebackTableEntry::GetLinks() const
{
  return c_link;
}
unsigned int TracebackTableEntry::GetNumOfIndels() const
{
  return c_numOfIndels;
}
unsigned int TracebackTableEntry::GetNumOfIndels_s() const
{
  return c_numOfIndels_s;
}  
void TracebackTableEntry::SetLinks(const LinkBack& _link)
{
  c_link=_link;
}
void TracebackTableEntry::SetNumOfIndels(const unsigned int& _numOfIndels)
{
  c_numOfIndels=_numOfIndels;
}
void TracebackTableEntry::SetNumOfIndels_s(const unsigned int& _numOfIndels)
{
  c_numOfIndels_s=_numOfIndels;
}
bool TracebackTableEntry::GetPathUsageState() const
{
	return c_used;
}
void TracebackTableEntry::SetPathUsageState( const bool& s) 
{
	c_used=s;
}

//
TracebackTable::TracebackTable(const SequenceString* _pattern, const SequenceString* _subject, const short& _typeOfAlignment):c_tbt(NULL), 
			     c_lenOfPattern(_pattern->GetLength()), c_lenOfSubject(_subject->GetLength() )
{
  //now allocate memory
  c_tbt=new TracebackTableEntry[(c_lenOfPattern+1)*(c_lenOfSubject+1)];
  
  
  c_tbt[0].SetLinks(ZERO);
  switch(_typeOfAlignment)
    {
    case 0: //global
      for(unsigned int i=1;i<=c_lenOfPattern;i++)
	{
	  c_tbt[i+0*(c_lenOfPattern+1)].SetLinks(LEFT);
	  c_tbt[i+0*(c_lenOfPattern+1)].SetNumOfIndels_s(i);
	  c_tbt[i+0*(c_lenOfPattern+1)].SetNumOfIndels(0);
	}
      for(unsigned int j=1;j<=c_lenOfSubject;j++)
	{
	  c_tbt[0+j*(c_lenOfPattern+1)].SetLinks(UP);
	  c_tbt[0+j*(c_lenOfPattern+1)].SetNumOfIndels(j);
	  c_tbt[0+j*(c_lenOfPattern+1)].SetNumOfIndels_s(0);
	}
      break;
    case 1: //overlap
    case 2: //local
      for(unsigned int i=1;i<=c_lenOfPattern;i++)
	{
	  c_tbt[i+0*(c_lenOfPattern+1)].SetLinks(ZERO);
	  c_tbt[i+0*(c_lenOfPattern+1)].SetNumOfIndels(0);
	  c_tbt[i+0*(c_lenOfPattern+1)].SetNumOfIndels_s(0);
	}
      for(unsigned int j=1;j<=c_lenOfSubject;j++)
	{
	  c_tbt[0+j*(c_lenOfPattern+1)].SetLinks(ZERO);
	  c_tbt[0+j*(c_lenOfPattern+1)].SetNumOfIndels(0);
	  c_tbt[0+j*(c_lenOfPattern+1)].SetNumOfIndels_s(0);
	}
      break;
    default:
      cout<<"ERROR: not recognized alignment type. Only support global, local and overlap so far"<<endl;
      break;
    }
  
}
TracebackTable::~TracebackTable()
{
  if(c_tbt!=NULL)
    delete [] c_tbt; 
}
void TracebackTable::SetTableEntry(const unsigned int& _patternIndex, const unsigned int& _subjectIndex, const LinkBack& _link, 
		const unsigned int& _numOfIndels , const unsigned& _numOfIndels_s)
{
  c_tbt[_patternIndex+_subjectIndex*(c_lenOfPattern+1)].SetLinks(_link);
  c_tbt[_patternIndex+_subjectIndex*(c_lenOfPattern+1)].SetNumOfIndels(_numOfIndels);
  c_tbt[_patternIndex+_subjectIndex*(c_lenOfPattern+1)].SetNumOfIndels_s(_numOfIndels_s);
}

LinkBack TracebackTable::GetLink(const unsigned int& _patternIndex, const unsigned int& _subjectIndex) const
{
  return c_tbt[_patternIndex+_subjectIndex*(c_lenOfPattern+1)].GetLinks();
}
//this is used for both global and overlap gap. for localalignment this is for pattern gap 
unsigned int TracebackTable::GetNumOfIndels(const unsigned int& _patternIndex, const unsigned int& _subjectIndex) const
{
  return c_tbt[_patternIndex+_subjectIndex*(c_lenOfPattern+1)].GetNumOfIndels();
}
//for local alignment to get --subject-- gap. added lately and used to distinguish the gaps 
unsigned int TracebackTable::GetNumOfIndels_s(const unsigned int& _patternIndex, const unsigned int& _subjectIndex) const
{
  return c_tbt[_patternIndex+_subjectIndex*(c_lenOfPattern+1)].GetNumOfIndels_s();
}

bool TracebackTable::GetPathUsageState(const unsigned int& _patternIndex, const unsigned int& _subjectIndex) const
{
	return c_tbt[_patternIndex+_subjectIndex*(c_lenOfPattern+1)].GetPathUsageState();
}

void TracebackTable::SetPathUsageState(const unsigned int& _patternIndex, const unsigned int& _subjectIndex, bool s)
{
	c_tbt[_patternIndex+_subjectIndex*(c_lenOfPattern+1)].SetPathUsageState(s);
}