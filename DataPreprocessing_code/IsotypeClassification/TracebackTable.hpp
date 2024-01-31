#ifndef TRACEBACKTABLE_HPP
#define TRACEBACKTABLE_HPP

#include "Accessory/SequenceString.hpp"

//this is struct to record the link back for tracing
//zero is for local align to indicate this is a zero node for terminating a path
enum LinkBack
  {
    UP, LEFT, UPLEFT, ZERO, UP_UPLEFT, LEFT_UPLEFT, 
		UP_LEFT //special case where could happend and both cases join the identical path.
				//see "GGGATATATA" vs"GGGTATATAT"
		, UP_LEFT_UPLEFT		
		//probably there will never be UP_LEFT_UPLEFT three way link back !!!!need to confirm.
		//----confirmed, we do have UP_LEFT_UPLEFT, three way link;
  };  //there probably no such thing called left_upleft_up (3 ways).!!!need to take care there more later.


class TracebackTableEntry
{
public:
  TracebackTableEntry();
  TracebackTableEntry(const LinkBack& _link, const unsigned int& _numOfIndels=0, const unsigned int& _numOfIndels_s= 0); 
  
  LinkBack GetLinks() const;
  unsigned int GetNumOfIndels() const; //get gap on pattern 
  unsigned GetNumOfIndels_s() const; //get gap on subject 
  
  void SetLinks(const LinkBack& _link);
  void SetNumOfIndels(const unsigned int& _numOfIndels); //for pattern gap 
  void SetNumOfIndels_s(const unsigned int & _numOfIndels);//for subject gap 

  bool GetPathUsageState() const;
  void SetPathUsageState( const bool& s=true) ;
private:
  LinkBack c_link;//this is pointer pointing to the one leading to this, could be left, up, upleft, zero
  unsigned int c_numOfIndels;//this is only works for left or up(indels),using to indicate how indels leads to this,not necessarily only 1.
  unsigned int c_numOfIndels_s; //this is added 1/26/2020 for added second numOfIndels because it could possibly 
						//both pattern gap and subject gap. so we using c_numOfIndels_s for subject gap and the first for the pattern gap.
						//--- gaps on subject. 
						//the first one is the gap on pattern 
  bool c_used; //this is a new field added on 1/19/2020 in order to indicate whehter this current entry has been 
			//used by previous path. use this one by traceback to avoid mainly the case wehre
			//one entry can participate at multiple path. then when the node is used by
			//one large value of path, it will be exclude from the later path.
			//this way, there no overlapping in path. 
};

class TracebackTable
{
public:
  //typeOfAlignment, 0, global;1, overlap; 2, local.
  TracebackTable( const SequenceString* _pattern, const SequenceString* _subject, const short& _typeOfAlignment=0);
  ~TracebackTable();
  void SetTableEntry(const unsigned int& _patternIndex, const unsigned int& _subjectIndex, const LinkBack& _link, 
		const unsigned int& _numOfIndels, const unsigned& _numOfIndels_s=0);

  LinkBack GetLink(const unsigned int& _patternIndex, const unsigned int& _subjectIndex) const;
  bool GetPathUsageState(const unsigned int& _patternIndex, const unsigned int& _subjectIndex) const;
  void SetPathUsageState(const unsigned int& _patternIndex, const unsigned int& _subjectIndex, bool s=true);
  
  //this was used for both global and overlap align. 
	//for localalignment we use this for gap on pattern.
  unsigned int GetNumOfIndels(const unsigned int& _patternIndex, const unsigned int& _subjectIndex) const ;
  //now this was added recently, for localalignment since we can have both link for the same node.
  //so we have to remember things differently
  unsigned int GetNumOfIndels_s(const unsigned int& _patternIndex, const unsigned int& _subjectIndex) const ;
  
private:
  TracebackTableEntry* c_tbt;
  unsigned int c_lenOfPattern;
  unsigned int c_lenOfSubject;
};

#endif
