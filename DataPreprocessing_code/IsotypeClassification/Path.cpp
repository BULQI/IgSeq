#include "Path.hpp"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <functional>
#include <sstream>
#include <numeric>

#include "GapModel.hpp"
#include "AffineGapModel.hpp"
#include "MarkovChainGapModel_454.hpp"

#define DEBUG_D

//***********For Path class
Path::Path():c_optimalValue(0), c_traced(false)
{
  c_optimalIndex[0]=0;
  c_optimalIndex[1]=0;
  c_startIndex[0]=0;
  c_startIndex[1]=0;
  
}
Path::Path( const unsigned int _OptimalIndex[2], const unsigned int _startIndex [2], const double& _optimalValue ):
  c_optimalValue(_optimalValue), c_traced(false)
{
  //empty here;
  c_optimalIndex[0]=_OptimalIndex[0];
  c_optimalIndex[1]=_OptimalIndex[1];
  c_startIndex[0]=_startIndex[0];
  c_startIndex[1]=_startIndex[1];
}
void Path::SetOptimalIndex(unsigned int _OptimalIndex[2])
{
    c_optimalIndex[0]=_OptimalIndex[0];
  c_optimalIndex[1]=_OptimalIndex[1];
  
}
void Path::SetStartIndex(unsigned int _startIndex [2])
{
  c_startIndex[0]=_startIndex[0];
  c_startIndex[1]=_startIndex[1];
}
void Path::SetOptimalValue(const double& _optimalValue)
{
  c_optimalValue=_optimalValue;
}
const unsigned int* Path::GetOptimalIndex() const
{
  return c_optimalIndex;
}
const unsigned  int* Path::GetStartIndex() const
{
  return c_startIndex;
}
double Path::GetOptimalValue() const
{
  return c_optimalValue;
}

bool Path::isTraced() const
{
	return c_traced;
}

void Path::setTraced() 
{
	this->c_traced=true; 
}

//constructor empty
PathElementEntry::PathElementEntry()
{
	//empty, all vectors, no need to initialize
}
//destructor no need to call ??	
PathElementEntry::~PathElementEntry()
{
	//empty, no need to destruct
}
	
unsigned PathElementEntry::GetNumberOfPathes() const
{
	return c_pathID.size();
}

unsigned PathElementEntry::GetPathID(unsigned index) const
{
	return c_pathID.at(index);
}
double PathElementEntry::GetPathOptimalScore(unsigned index) const
{
	return c_optimalScore.at(index);
}

unsigned PathElementEntry::GetPathOptimalIndexPattern(unsigned index) const
{
	return c_optimalIndexPattern.at(index);
}

unsigned PathElementEntry::GetPathOptimalIndexSubject(unsigned index) const
{
	return c_optimalIndexSubject.at(index);
}

	//this here we add one single path to the current node pathelement
	//
	//in here, we have to set up thing at the same time together
	// input index: which path this is for the entry, like first, second, third....
	//			pathID: the path ID in the path vector
	//		optimalScore: the optimal score so far for the path
	//		optimalIndexPattern/subject : the index for the optimal score entry;
	//
	//This one is checking the dupliation.
void PathElementEntry::AddPathInfoNoCheckingDuplication(const unsigned& index, const unsigned& pathID, const double& optimalScore, 
			const unsigned& optimalIndexPattern, const unsigned& optimalIndexSubject)
{
//------check for duplication--------
	if(index>=c_pathID.size())
	{
		//this is a new one, so we add
		c_pathID.push_back(pathID);
		c_optimalScore.push_back(optimalScore);
		c_optimalIndexPattern.push_back(optimalIndexPattern);
		c_optimalIndexSubject.push_back(optimalIndexSubject);
	}
	else  //---- check for duplication--------
	
	{
		//reset the value 
		c_pathID.at(index)=pathID;
		c_optimalScore.at(index)=optimalScore;
		c_optimalIndexPattern.at(index)=optimalIndexPattern;
		c_optimalIndexSubject.at(index)=optimalIndexSubject;
	}
}
//copy constructor, create from new, so no need to do checking for duplicate 
//Do deep copy. the vector assignment is to do the deep copy
PathElementEntry::PathElementEntry(const PathElementEntry& pe)
{
	//deep copy
	c_pathID=pe.c_pathID;
	c_optimalScore=pe.c_optimalScore;
	c_optimalIndexPattern=pe.c_optimalIndexPattern;
	c_optimalIndexSubject=pe.c_optimalIndexSubject;
}
//add the input node path info to this current with checking duplication.
//assuming the incoming one there are no duplicated pathes
//---not valid below. it should check for duplication. don't mind about multiple through the same node.
//we should not do check for duplicates, since there are chances that
//the same path goes through the same nodes multiple times, and we 
//want to include all passes there in the record.
void PathElementEntry::AddPathInfo(const PathElementEntry& pe)
{
	//cout<<"current path element entray:"<<toString()<<endl;
	unsigned input_numPath=pe.GetNumberOfPathes();
	unsigned input_pathID;
	unsigned curr_numPath=this->c_pathID.size();
	//cout<<"strating calling adding in"<<endl;
	//cout<<"input_numPath:"<<input_numPath<<"; path id:"<<input_pathID<<";current num path:"<<curr_numPath<<endl;
	for(unsigned i=0;i<input_numPath;i++)
	{
		//cout<<"node i input:"<<i<<endl; 
		input_pathID=pe.GetPathID(i);
		//check whether 
		bool duplicated=false;
		//need to check for duplicate only the current one, but not new incoming one.
		for(unsigned j=0;j<curr_numPath;j++)
		{
			//cout<<"node j current:"<<j<<endl;
			if(input_pathID == this->GetPathID(j))
			{
				//cout<<"updating........"<<endl;
				duplicated=true;  //not valid now----we made a trick here, disabling duplicate checking
								//not valid now----change this back to true and uncomment the following two 
								//not valid now----lines to make it checking duplicates 
				//update
				this->AddPathInfoNoCheckingDuplication(j, input_pathID, 
				pe.GetPathOptimalScore(i), pe.GetPathOptimalIndexPattern(i), pe.GetPathOptimalIndexSubject(i));
				//cout<<"updating....done>>>>>>>"<<endl;
				break;
			}
		}
		if(duplicated)
		{//don't do anything
			
			continue;
		}
		//cout<<"done for one input"<<endl;
		//otherwise we need to add this current node path information.
		this->c_pathID.push_back(input_pathID);
		this->c_optimalScore.push_back(pe.GetPathOptimalScore(i));
		this->c_optimalIndexPattern.push_back(pe.GetPathOptimalIndexPattern(i));
		this->c_optimalIndexSubject.push_back(pe.GetPathOptimalIndexSubject(i));
	}
}	
//printing out the 
string PathElementEntry::toString() const
{
	ostringstream s;
	s<<"...Path Entry...(";
	s<<this->GetNumberOfPathes();
	s<<" pathes):"<<endl;
	//s<<"\n------";

	for(unsigned i=0;i<this->GetNumberOfPathes();i++)
	{
		s<<"\t------";
		s<<"(pathID:"<<this->GetPathID(i);
		s<<"; Optimal score:"<<this->GetPathOptimalScore(i);
		s<<"; pattern index:"<<this->GetPathOptimalIndexPattern(i);
		s<<"; subject index:"<<this->GetPathOptimalIndexSubject(i);
		s<<")\n";
	}
	return s.str();
}

//remove a path (indicated by "index") from the records.
//input: index is the location of the record to be removed, not
//		the pathid.
void PathElementEntry::RemovePath(const unsigned index)
{
	//double check for the correctness
	if(index>=c_pathID.size())
	{
		cerr<<"ERROR: the index is out of bound in PathElementEntry::RemovePath"<<endl;
	}
	c_pathID.erase(c_pathID.begin()+index);
	c_optimalScore.erase(c_optimalScore.begin()+index);//used to indicate up to this current entry the optimal score for the specific path
	c_optimalIndexPattern.erase(c_optimalIndexPattern.begin()+index);//used to indicate where the optimal score on pattern is.
	c_optimalIndexSubject.erase(c_optimalIndexSubject.begin()+index);//same thing, but on subject.

}
//based on the pathid, look up the index of the path in the recored (c_pathID). the record is for each node!!! there are might be multiple path through the node.
unsigned PathElementEntry::LookUpPathIndex(const unsigned& pathID) const 
{
	for(unsigned i=0;i<c_pathID.size();i++)
	{
		if(pathID==c_pathID.at(i))
			return i;
	}
	
	return (unsigned)-1;
}
 
