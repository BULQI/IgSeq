#include "LocalAlignment_CT.hpp"
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
#include "Path.hpp"

//#define DEBUG_D
//#define DEBUG

using namespace std;
/*static bool comparePathElement(Path* p1, Path* p2)
{
  return p1->GetOptimalValue() > p2->GetOptimalValue();
}*/

//In this function, we try to combine the pathid array using 
//the trick to bind the path to the pathID data.
//a normal compare function used by the sort template is taking
//two arguments, see above regular comare function for comparing the pathelement vector. 
//this current design enables us to compare the path value but only sort the 
//pathID vector without affecting the pathelement array. so we can easily/quickly 
//access to them using the pathid. pathid in this case is the index of the 
//pathelement vector. 
bool comparePathElementIndex(const unsigned & id1, const unsigned& id2, vector<Path*>& ptb)
{
	return ptb.at(id1)->GetOptimalValue() > ptb.at(id2)->GetOptimalValue();
}

LocalAlignment_CT::LocalAlignment_CT()
{
	//disable by "protect" it.
}

LocalAlignment_CT::LocalAlignment_CT(SequenceString* _pattern, SequenceString* _subject, 
			       const ScoreMatrix* _m, const double& _gopen, 
			       const double& _gextension, const double& _scale, const int& _numOfAlignments, const short& _typeOfGapModel):
  PairwiseAlignment(_pattern, _subject, _m, _gopen, _gextension, _scale,_typeOfGapModel),c_numOfAlignments(_numOfAlignments),
  c_alignmentArr(NULL), c_scoreArr(NULL), c_PathElementTable(NULL)
{
  //now we need to
  c_alignmentArr=new AlignmentString[this->c_numOfAlignments];
  c_scoreArr=new double[this->c_numOfAlignments];
  //dp_table=NULL;
  //traceback_table=NULL;
#ifdef DEBUG
  cout<<"LocalAlignment_CT: Start doing align()......"<<endl;
#endif
  this->align();
#ifdef DEBUG
  cout<<"LocalAlignment_CT: start doing Trackback()...."<<endl;
 
  //flush(cout);
#endif
  this->traceBackMultiple();
#ifdef DEBUG
  cout<<"LocalAlignment_CT: done!!!....."<<endl;
#endif
}
  
LocalAlignment_CT::~LocalAlignment_CT()
{
  //the base class destructor is called automatically
  //so trace table is taken care of there 
  //	delete [] c_traceback_table;
  
  
  //here we only need to take care other issues
  if(c_alignmentArr!=NULL)
    delete[] c_alignmentArr;
  if(c_scoreArr!=NULL)
    delete[] c_scoreArr;

  
  if (c_PathElementTable!=NULL)
  {
	  /*
	  //need to deep destructing
	  if(c_pattern!=NULL&&c_subject!=NULL)
	  {
		//size mut be 
		for(unsigned i =0;i<c_pattern->GetLength()+1;i++)
		{
			for(unsigned j =0;j<c_pattern->GetLength()+1;j++)
			{
				if(c_PathElementTable[i+(j*c_pattern->GetLength())]!=NULL)
				{
					delete c_PathElementTable[i+(j*c_pattern->GetLength())];
				}
			}
		}
	  }*/
	  delete [] c_PathElementTable;
  }
  //taking care of Path vector
  for(unsigned int i=0;i<c_path_vec.size();i++)
    {
      //cout<<"deleting path vecs"<<endl;
	  if(c_path_vec.at(i)!=NULL)
		delete c_path_vec.at(i);
    }
}

//this is 
double* LocalAlignment_CT::GetScoreArr()
{
  return c_scoreArr;
}

unsigned int LocalAlignment_CT::GetNumberOfAlignments()
{
  return c_numOfAlignments;
}

AlignmentString* LocalAlignment_CT::GetAlignmentArr()
{

  return c_alignmentArr;
}

/*
//this code is not working because of later changes on the variables and etc. but the idea is there.
void LocalAlignment_CT::align() ----this one is working one, but keeping track of a full dp table.is replaced by the new low memory one.
{
  	//Create empty dynamic programming table
  unsigned lenP=c_pattern->GetLength();
  unsigned lenS=c_subject->GetLength();
  c_dp_table=new double[(lenP+1)*(lenS+1)];//one extra on each dimension, to deal with the beging of the row and column
  c_traceback_table=new LinkBack[(lenP+1)*(lenS+1)];

  c_score=-10E9;
  

  //intialize the dp table, the first row=0 and first column=0
  //be really careful that the dp table is a linear array, we need to 
  //mannually manage it
  //the row is Pattern, column is Subject. 
  for(unsigned int i=0;i<=lenP;++i)
    {
      cout<<"doing dp table......"<<i<<endl;
      c_dp_table[i]=0;//this is the first row, like c_dp_table[0][i]
    }
  for(unsigned int i=0;i<=lenS;++i)
    {
      c_dp_table[i*lenP]=0;//this is the first column, like c_dp_table[i][0]=0;
    }

  cout<<"successfully created the empty tables and now go to get the score!\n";

  //score table with S-W
  double compval = 0;
  string strP=c_pattern->GetSequence();
  string strS=c_subject->GetSequence();


  //we have to be very careful here
  //the dimension of the dp table and the strings are not same.
  //dp table is one row/col more than the strings, since we need to 
  //have the zero/starting row or column.so that means, when we dealing with
  //dp table, it is from 1->lenP or 1->lenS, but the string, starts from zero,
  //so we are doing things i-1
  for(unsigned int i = 1; i <= lenP; ++i) 
    {	//for all values of strA
		
      for(unsigned int j = 1; j <= lenS; ++j) 
	{	//for all values of strB
	  cout<<"****doing round ("<<i<<","<<j<<")."<<endl;	
	  //MATCH
	  //if(strP.at(i-1) == strS.at(j-1)) 
	  //{				//if current sequence values are the same
	  cout<<"\tcalling score matrix:score("<<strP.at(i-1)<<","<<strS.at(j-1)<<")="<<c_sm->GetScore(strP.at(i-1),strS.at(j-1))<<endl;
	  cout<<"\t\tdp table is dp("<<i-1<<","<<j-1<<")="<<c_dp_table[i-1+(j-1)*lenP]<<endl;
	  compval = (c_dp_table[i-1+(j-1)*lenP] + c_sm->GetScore(strP.at(i-1),strS.at(j-1)));	//compval = diagonal + match score
	  //}
			
	  c_traceback_table[i+j*lenP]=UPLEFT;
	  for(int k = i-1; k > 0; --k) 
	    {		//check all sub rows
				
	      if(compval < ((c_dp_table[k+j*lenP]) + (c_gapOpen + (c_gapExtension *(i- k))))) 
		{	    //if cell above has a greater value 
					
		  compval = ((c_dp_table[k+j*lenP]) + (c_gapOpen + (c_gapExtension *(i- k))));		//set compval to that square
		  c_traceback_table[i+j*lenP]=LEFT;
		}
	    }

	  for(int k=j-1; k>0; --k) 
	    {		//check all sub columns
				
	      if(compval < ((c_dp_table[i+k*lenP]) + (c_gapOpen + (c_gapExtension *(j- k))))) 
		{	
		  //if square to the left has the highest valu
					
		  compval = ((c_dp_table[i+k*lenP]) + (c_gapOpen + (c_gapExtension *(j- k))));    //set compval to that square
		  c_traceback_table[i+j*lenP]=UP;
		}
	    }		
			
	  if(compval - 0<1E-9) 
	    {
	      compval = 0;
	      c_traceback_table[i+j*lenP]=ZERO;
	    }
	  
	  c_dp_table[i+j*lenP] = compval;	//set current cell to highest possible score and move on
	  cout<<"\t\trunning score:"<<compval<<endl;

	  //find a best one so far
	  if(c_score<compval)
	    {
	      c_score=compval;
	      c_optimalIndex[0]=i;
	      c_optimalIndex[1]=j;
	    }
	}
    }
  
}//end of the localAlign
*/

void LocalAlignment_CT::traceBackMultiple()
{
  //decide how many we will be tracing back 
  //unsigned int actualNumOfAlignments;
#ifdef DEBUG
  cout<<"+++++++++doing the multiple trce bck "<<endl;
  cout<<"c_path_vec.size():"<<c_path_vec.size()<<endl;
	flush(cout);
#endif 
  if(c_numOfAlignments>c_path_vec.size())
    {//can only do what there are.
      c_numOfAlignments=c_path_vec.size();
    }
  if(c_numOfAlignments==0)
    c_numOfAlignments=1;

  if(c_path_vec.size()==0)
	c_numOfAlignments=0;
  if(c_numOfAlignments==0)
	{
		c_alignmentArr=NULL;
		c_scoreArr=NULL;
		return ;
	}
  c_alignmentArr=new AlignmentString[c_numOfAlignments];
  c_scoreArr=new double[c_numOfAlignments];

  //showing the vector information
 #ifdef DEBUG
  cout<<"******************SHOWING PATH INFO**************"<<endl;
  cout<<"\ttotal num of path:"<<c_path_vec.size()<<endl;
  for(unsigned int i=0;i<c_path_vec.size();i++)
    {
      cout<<i<<"/"<<c_path_vec.size()<<": start at ("<<c_path_vec.at(i)->GetStartIndex()[0]<<","
	  <<c_path_vec.at(i)->GetStartIndex()[1]<<") and end at ("<<c_path_vec.at(i)->GetOptimalIndex()[0]
	  <<","<<c_path_vec.at(i)->GetOptimalIndex()[1]<<"); score of :("<<c_path_vec.at(i)->GetOptimalValue()<<";\n";
    }
#endif
  //now sort the vector of PathID, we will not change the order of the c_path_vec
  //since the 
  c_pathID_vec.resize(c_pathID_vec.size());
  iota(begin(c_pathID_vec), end(c_pathID_vec),0);
  sort(c_pathID_vec.begin(), c_pathID_vec.end(), 
		bind(comparePathElementIndex, placeholders::_1, placeholders::_2, c_path_vec));
	//cout<<"maxi one :"<<c_path_vec.at(c_pathID_vec.at(0))->GetOptimalValue()<<endl;
	//cout<<"maxi one2 :"<<c_path_vec.at(c_pathID_vec.at(1))->GetOptimalValue()<<endl;
//showing the vector information
/*  
  cout<<"******************SHOWING PATH INFO**************"<<endl;
  cout<<"\ttotal num of path:"<<c_path_vec.size()<<endl;
  for(unsigned int i=0;i<c_pathID_vec.size();i++)
    {
      cout<<i<<"/"<<c_path_vec.size()<<": start at ("<<c_path_vec.at(c_pathID_vec.at(i))->GetStartIndex()[0]<<","
	  <<c_path_vec.at(c_pathID_vec.at(i))->GetStartIndex()[1]<<") and end at ("<<c_path_vec.at(c_pathID_vec.at(i))->GetOptimalIndex()[0]
	  <<","<<c_path_vec.at(c_pathID_vec.at(i))->GetOptimalIndex()[1]<<"); score of :"<<c_path_vec.at(c_pathID_vec.at(i))->GetOptimalValue()<<";\n";
    }
  */
//cout<<"\tefore looping......."<<endl;
//flush(cout);

  for(unsigned int i=0;i<c_numOfAlignments;i++)
    {
#ifdef DEBUG
		cout<<"i:"<<i<<"..."<<endl;
#endif
      //prepare the input and output for trackack() function
      c_optimalIndex[0]=c_path_vec.at(c_pathID_vec.at(i))->GetOptimalIndex()[0];
      c_optimalIndex[1]=c_path_vec.at(c_pathID_vec.at(i))->GetOptimalIndex()[1];
	  c_current_path_tk=c_pathID_vec.at(i); //path ID for the one currently being checked/traced back 

#ifdef DEBUG	  
	  cout<<"//*************showing path info************"<<endl;
	  cout<<"Path start at ("<<c_path_vec.at(c_pathID_vec.at(i))->GetStartIndex()[0]<<","
	  <<c_path_vec.at(c_pathID_vec.at(i))->GetStartIndex()[1]<<") and end at ("<<c_path_vec.at(c_pathID_vec.at(i))->GetOptimalIndex()[0]
	  <<","<<c_path_vec.at(c_pathID_vec.at(i))->GetOptimalIndex()[1]<<"); score of :("<<c_path_vec.at(c_pathID_vec.at(i))->GetOptimalValue()<<";\n";
      
	  cout<<"before doing trace back;"<<endl;
	  //flush(cout);
#endif
      //now call the traceBack to do the job
      traceBack();
      
	  //for the new algorithm, we need to sort again after each trace
	  //the good news is that we only need to sort the rest instead of the one we have done with tracing back;
	  if(i<c_numOfAlignments-1)
	  {
		sort(c_pathID_vec.begin()+i+1,c_pathID_vec.end(),  
				bind(comparePathElementIndex, placeholders::_1, placeholders::_2, c_path_vec));
	  }
      //now record the output of trace back
      c_alignmentArr[i].SetPattern(c_alignment.GetPattern(true), true);
      c_alignmentArr[i].SetPattern(c_alignment.GetPattern(false), false);
  
      c_alignmentArr[i].SetSubject(c_alignment.GetSubject(true), true);
      c_alignmentArr[i].SetSubject(c_alignment.GetSubject(false), false);
  
      c_alignmentArr[i].SetPatternIndex(c_alignment.GetPatternIndexStart(), c_alignment.GetPatternIndexEnd());
      c_alignmentArr[i].SetSubjectIndex(c_alignment.GetSubjectIndexStart(), c_alignment.GetSubjectIndexEnd());
      c_alignmentArr[i].SetScore(c_path_vec.at(c_pathID_vec.at(i))->GetOptimalValue());
      c_scoreArr[i]=c_path_vec.at(c_pathID_vec.at(i))->GetOptimalValue();
    }

  //after set up everything, we want to set the c_alignment to the best alignment one,
  //we don't need this for localAlign, but for compatibility purpose only.
	c_alignment.SetPattern(c_alignmentArr[0].GetPattern(true), true);
      c_alignment.SetPattern(c_alignmentArr[0].GetPattern(false), false);
  
      c_alignment.SetSubject(c_alignmentArr[0].GetSubject(true), true);
      c_alignment.SetSubject(c_alignmentArr[0].GetSubject(false), false);
  
      c_alignment.SetPatternIndex(c_alignmentArr[0].GetPatternIndexStart(), c_alignmentArr[0].GetPatternIndexEnd());
      c_alignment.SetSubjectIndex(c_alignmentArr[0].GetSubjectIndexStart(), c_alignmentArr[0].GetSubjectIndexEnd());
      c_alignment.SetScore(c_path_vec.at(0)->GetOptimalValue());
      c_score=c_path_vec.at(c_pathID_vec.at(0))->GetOptimalValue();
	  
	//cout<<"+++++++++done the multiple trace back"<<endl;
	flush(cout);
  
  //done!!!
}

//----!!!NO OBsoleted NO, we should use the one for no checking the duplicates
//see below.  (see below function for reasoning. )
//this is used to update the vector of pathId so that we can remember
//to update the pathes in it later. we only record the path Id in this 
//vector. note, this pathid record is not the one in the object/class
//we only need this temporarily for check each path (traceback())
void insert_checkingDuplicate(vector<unsigned>& v, const unsigned& id
			)
{
	bool already_in=false;
	for(unsigned i=0 ; i < v.size(); ++i)
	{
		if(id==v[i])
		{
			//we already has it in the vector,
			//so break out
			already_in=true;
			break;
		}
	}
	if(!already_in)
	{
		v.push_back(id);
	}
}
//-----
//this is used to update the vector of pathId so that we can remember
//to update the pathes in it later. we only record the path Id in this 
//vector. note, this pathid record is not the one in the object/class
//we only need this temporarily for check each path (traceback())
//in this case, we don't do checking for duplicates. this is a good way, 
//since the path might coming from different ways, so we could end up with 
//multiple count of the same path through the nodes. so keep it as a 
//count.
//----Updated again, we need to check for duplicate. this is different from 
//add pathes to the nodes, this is used to indicated which path we need to 
//check/update. although it is likely to have same path goes through the same 
//nodes multiple time, we only need to check the path onece. since we 
//will check all passes of the path .
//also when we add the pathes to the nodes we check duplicate, (see 
// AddPathInfo() function)
//
void insert_NoCheckingDuplicate(vector<unsigned>& v, const unsigned& id
			)
{
		v.push_back(id);
}

//in this function, we check the current node that is being used by the current 
//path and then look for the other pathes that is going through this nodes.
//here we simply write them down to a temporary vector for later update. By later
//we mean after trace back this current path.
//we set it up as a object function, since we need to get access to the information
//about the class objec, the c_pattern length and the pathElementTable. 
void LocalAlignment_CT::recordPathInfoToBeRemove(const unsigned& _patternIndex, 
				const unsigned & _subjectIndex
			, const unsigned& c_current_path_tb,
					//outputs 
				vector<unsigned>& vec_pathToBeChecked 
				)
{//in here, we simply go through the current node and get the path which is not the current one 
//and save it in the vector  
	 unsigned lenP=c_pattern->GetLength();
	 //unsigned lenS=c_subject->GetLength();
	 //get the path ids from the current node
	 unsigned numOfPathes=c_PathElementTable[_patternIndex+_subjectIndex*(lenP+1)].GetNumberOfPathes();
	 
	 for(unsigned i=0;i<numOfPathes;i++)
	 {
		 unsigned temp_PathID=c_PathElementTable[_patternIndex+_subjectIndex*(lenP+1)].GetPathID(i);
		 if(temp_PathID!=c_current_path_tb)
		 {
			 
			 //vec_pathToBeChecked.push_back(running_PathID);
			 //vec_pathToBeChecked_StartIndex_Pattern.push_back(_patternIndex);
			 //vec_pathToBeChecked_StartIndex_Subject.push_back(_subjectIndex);
			 if(!c_path_vec.at(temp_PathID)->isTraced())  //we will not update the path that has been traced back previousely, only the ones has not been done with tracking back. 
				insert_checkingDuplicate(vec_pathToBeChecked, temp_PathID);
		 }
	 }
}

//this is called each time we trace back one path. The rationale is that after the tracing back of 
//of one path (current_pathID), we have used some nodes in the table. these nodes can have 
//multiple links, so in this way, some other path become dead and not valid anymore. For 
//those pathes that were involved we have "remember" them in a vector (vec_pathToBeChecked).
//now we need to go through the path vector to update in case 
//some node has been taken by the current path. 
//this is the driving function will call
//the recursive function to update the whole path for those remebered during the tracing back. 
//input: 
//	current_pathID, the ID for the currently traced back.
//	vec_pathToBeChecked, the pathes that has been "remembered" in the previous tracing back.
//					//
void LocalAlignment_CT::updatePathForMultipleLinks(const unsigned& current_pathID, 
			vector<unsigned>& vec_pathToBeChecked)
{
	//go through the paths 
	unsigned running_pathID; //to path to be check for updating
	unsigned running_starting_index_s;
	unsigned running_starting_index_p;
	double running_optimalScore;
	unsigned lenP=c_pattern->GetLength();
#ifdef DEBUG
	cout<<"UpdatedPathForMultipleLinks:current path ID:"<<current_pathID<<endl;
	cout<<"total number of pathes to be checked:"<<vec_pathToBeChecked.size()<<endl;
#endif
	for(unsigned i=0;i<vec_pathToBeChecked.size();i++)
	{
#ifdef DEBUG
		cout<<" checking ith path:"<<i<<endl;
#endif
		//following the path and starting at the path starting indexes.
		running_pathID=vec_pathToBeChecked.at(i);
#ifdef DEBUG		
		cout<<" checking path id :"<<running_pathID<<endl;
		cout<<"//*************showing path info************"<<endl;
	  cout<<"Path start at ("<<c_path_vec.at(running_pathID)->GetStartIndex()[0]<<","
	  <<c_path_vec.at(running_pathID)->GetStartIndex()[1]<<") and end at ("<<c_path_vec.at(running_pathID)->GetOptimalIndex()[0]
	  <<","<<c_path_vec.at(running_pathID)->GetOptimalIndex()[1]<<"); score of :("<<c_path_vec.at(running_pathID)->GetOptimalValue()<<";\n";
#endif      
		running_starting_index_p=c_path_vec.at(running_pathID)->GetStartIndex()[0];
		running_starting_index_s=c_path_vec.at(running_pathID)->GetStartIndex()[1];
		
		//get the the current one best score so far, need to check the pathelementEntry table for 
		//the specific node
		PathElementEntry* p_pee=&(c_PathElementTable[running_starting_index_p + running_starting_index_s*(lenP+1)]);
		unsigned index_pee= p_pee->LookUpPathIndex(running_pathID);
		running_optimalScore= p_pee->GetPathOptimalScore(index_pee);//<-IMPORTANT, we need to get the current one optimal score so far.
#ifdef DEBUG		
		cout<<"running_optimalScore:"<<running_optimalScore<<endl;
		cout<<"calling the recurvise working horse...."<<endl;
#endif
		//now calling to chedk the nodes, recursively, assuming the initial node of the checking path is good. 
		//if we have done good job, we will never have the initial node of a path to be part of the previous path.!!!
		//this is why we set removeThis node to be false. it has to be ZERO link back, which is a starting node.
		//also we need to reset the running_path to be zero, and optimal index now is the starting pointer
		//c_path_vec.at(running_PathID)-
		running_optimalScore=0;
		updatePathForMultipleLinks_Node(running_pathID, running_starting_index_p, running_starting_index_s,
				running_optimalScore, false, ZERO);
#ifdef DEBUG
		cout<<"done.........."<<endl;
#endif
	}
	
}

//this function is used to check whether a path (denoted by a PathID) goes throught
//a node (denoted by index_p and index_s). A path not going through a node because it 
//is go through by other path uniquely with higher score during align() and it could also be caused due to 
//it is a node with multiple pathes, but were taken by previous passing of updating. 
// PathID is the same as the location/index of the path in the original c_path_vec formed after align();
// Updated 8/27/2021, add codition so a zero node is always returning true, meaning the "current" path through zero nodes
bool LocalAlignment_CT::isPathThruThisNode(const unsigned& pathID, const unsigned& index_p, const unsigned& index_s) const 
{
    //new add new condition, if this node is a zero node (LINKBACK), then this node is taken!!! return true.
    if(c_traceback_table->GetLink(index_p, index_s)==ZERO)
    {
        return true;
    }
	unsigned lenP=c_pattern->GetLength();
	PathElementEntry* p_pee=&(c_PathElementTable[index_p + index_s*(lenP+1)]);
	//go throught c_pathID_vec and look for the pathID
	for(unsigned i=0;i<p_pee->GetNumberOfPathes();i++)
	{
		if(p_pee->GetPathID(i)==pathID)
		{
			return true;
		}
	}
	return false;
}

//function to remove the pathID form the pathelementEntry table pathID_vec for the specific node
//   note, so that the node not in the path anymore, mainly because the nodes has been taken and no more
//available now.

void LocalAlignment_CT::removePathFromNode(const unsigned& pathID, const unsigned& index_p, const unsigned& index_s)
{
	unsigned lenP=c_pattern->GetLength();
	
	PathElementEntry* p_pee=&(c_PathElementTable[index_p + index_s*(lenP+1)]);
	//go throught c_pathID_vec and look for the pathID
	for(unsigned i=0;i<p_pee->GetNumberOfPathes();i++)
	{
		if(p_pee->GetPathID(i)==pathID)
		{
			p_pee->RemovePath(i);
			return ;
		}
	}
	return ;
}


//
//----Recursive--- function
//input: removThisPath, boolean, used to indicate that we need to remove the this path from this node
//		and also all the nodes following this current node, since the road is cut off 
//		anyway.
//	PathID, this current path being checked.
//	optimal score is the best score so far for the path being checked.  
//  link, this is the link of the caller. together with this current node link it is used to indicate the current node has to have the correct link back in order to ex
//			extend the path. the case might be even the node is on the path valid, but 
//			it is not linked from the immdiate node calling the recursive funciton. 
//			in this case, we still don't count it as a valid one. The caller will only
//			input 3 cases, upleft, up and left for positions. inside this function we 
//			will check for more complicated ones. 
void LocalAlignment_CT::updatePathForMultipleLinks_Node(const unsigned& pathID, const unsigned& index_p, const unsigned& index_s
				, const double& optimalScore, const bool& removeThisPath, const LinkBack& link)
			
{
	//check whether we are done for the current  path
	//		1) the p index can not go further (reaching the end)
	//		2) the s index can not go further (reaching the end)
	//				either one is good enough to return, since we are doing the local alignment 
	//				if we are reaching the side (either side) we will not go further to be
	//				better scored. so we can stop.
	//      3) we are reaching a ZERO node.
	//		4) the currentPath don't go through the node being checked. then we will stop. 
	//		5) the node has been used.???, we will remove the path from the pathElementTable entry from this node, but still will go down the links/path to update the downstream node for this path. 
#ifdef DEBUG	
	cout<<"%%%%%checking node("<<index_p <<","<<index_s<<");"<<endl;
	//Case 1 and 2
	cout<<"c_pattern->GetLength():"<<c_pattern->GetLength();
#endif
	if(index_p>c_pattern->GetLength()||index_s>c_subject->GetLength())
	{
#ifdef DEBUG
			cout<<"returning due to out of string (reaching side)"<<endl;
#endif
		return ;
	}
	//Case 3
	//cout<<"1:ZERO :"<<ZERO<<";current link:"<<c_traceback_table->GetLink(index_p, index_s)<<endl;
	if(c_traceback_table->GetLink(index_p, index_s)==ZERO)
	{
#ifdef DEBUG
		cout<<"returning duo to ZERO out of score"<<endl;
#endif
		return;
	}
	unsigned lenP=c_pattern->GetLength();
	//Case 4
	if(!isPathThruThisNode(pathID, index_p, index_s))
	{
#ifdef DEBUG
		cout<<c_PathElementTable[index_p+index_s*(1+lenP)].toString()<<endl;
		cout<<"returning due to no pass throug"<<endl;
#endif
		return ;
	}
	
	//case 6
	LinkBack current_node_link=c_traceback_table->GetLink(index_p, index_s);
	//cout<<"ZERO :"<<ZERO<<";current link:"<<c_traceback_table->GetLink(index_p, index_s)<<endl;
	//checking the links and
	
		switch(link)
		{
			case ZERO:
				if(current_node_link!=UPLEFT)
				{
#ifdef DEBUG
					cout<<"start one should UPLEFT, so returning...SOMETHING WRONG> it is impossilbe"<<endl;
#endif 
					return ;
				}
				break;
			case UP:
				if(current_node_link!=UP&&current_node_link!=UP_UPLEFT&&current_node_link!=UP_LEFT&&current_node_link!=UP_LEFT_UPLEFT)
				{
#ifdef DEBUG
					//this is no way to be a good node for this branch, break out
					cout<<"good node, but not pass through this branch, go return UP case"<<endl;
#endif
					return ;
				}
				break;
			case UPLEFT:
				if(current_node_link!=UPLEFT&&current_node_link!=UP_UPLEFT&&current_node_link!=LEFT_UPLEFT&&current_node_link!=UP_LEFT_UPLEFT)
				{
#ifdef DEBUG
					//this is no way to be a good node for this branch, break out
					cout<<"good node, but not pass through this branch, go return UPLEFT case"<<endl;
#endif
					return ;
				}
				break;
			case LEFT:
				if(current_node_link!=LEFT&&current_node_link!=LEFT_UPLEFT&&current_node_link!=UP_LEFT&&current_node_link!=UP_LEFT_UPLEFT)
				{
#ifdef DEBUG
					//this is no way to be a good node for this branch, break out
					cout<<"good node, but not pass through this branch, go return LEFT case"<<endl;
#endif
					return ;
				}
				break;
			default:
				cerr<<"this is impossible,quit"<<endl;
				exit(-1);
				break;
		}

	//for case 5 and other cases, we still will 	
	//if we are here, do the job
	double running_optimalScore;
	//unsigned running_optimalIndex[2];
	bool running_removeThisPath=false;
	//unsigned lenP=c_pattern->GetLength();
	//Path* current_path=c_path_vec.at(pathID);
	//update information 
	if(c_traceback_table->GetPathUsageState(index_p, index_s)||removeThisPath) //case #5 
	{  //used, so we need to get ride of the path from this node
#ifdef DEBUG
		cout<<"this node has been used by other higher priority path, but search is going on"<<endl;
#endif
		running_removeThisPath=true; 
		//remove the pathID from this nodes, for one try , we might come back again if path passes here many times. 
		removePathFromNode(pathID, index_p, index_s); 
#ifdef DEBUG
		cout<<c_PathElementTable[index_p+index_s*(1+lenP)].toString()<<endl;
#endif
	}
	else //this node has Not been used, so we are ok.
	{
#ifdef DEBUG
		cout<<"****good nodes for this path"<<endl;
#endif
		//update
		PathElementEntry* p_pee=&(c_PathElementTable[index_p+index_s*(lenP+1)]);
		unsigned index_pee=p_pee->LookUpPathIndex(pathID);
		running_optimalScore=p_pee->GetPathOptimalScore(index_pee);
		if(running_optimalScore>optimalScore)  //best one so far, update
		{
			//optimalScore=running_optimalScore;
			//update information,
			c_path_vec.at(pathID)->SetOptimalValue(running_optimalScore);
			unsigned temp_OptimalIndex[2];
			temp_OptimalIndex[0]=index_p; temp_OptimalIndex[1]=index_s;
			c_path_vec.at(pathID)->SetOptimalIndex(temp_OptimalIndex);
		}
		else  //no better than the optimal value of the path so far, remember the best one, and write it to the pathelemententry
		{//and get ready to call the next node along the path.
			running_optimalScore=optimalScore;
			//running_optimalIndex[0]=c_PathElementTable[index_p+index_s*(lenP+1)].GetPathOptimalIndexPattern(index_pee);
			//running_optimalIndex[1]=c_PathElementTable[index_p+index_s*(lenP+1)].GetPathOptimalIndexSubject(index_pee);
		}
	} //OK we are done with good job, 

#ifdef DEBUG	
	cout<<"//*************showing path info************"<<pathID <<endl;
	  cout<<"Path start at ("<<c_path_vec.at(pathID)->GetStartIndex()[0]<<","
	  <<c_path_vec.at(pathID)->GetStartIndex()[1]<<") and end at ("<<c_path_vec.at(pathID)->GetOptimalIndex()[0]
	  <<","<<c_path_vec.at(pathID)->GetOptimalIndex()[1]<<"); score of :("<<c_path_vec.at(pathID)->GetOptimalValue()<<";\n";
#endif      
	//update index and go next
	//go right 
	
	updatePathForMultipleLinks_Node(pathID, index_p+1, index_s
				, running_optimalScore, running_removeThisPath, LEFT);
	//go right down
	updatePathForMultipleLinks_Node(pathID, index_p+1, index_s+1
				, running_optimalScore, running_removeThisPath, UPLEFT);
	
	//go down
	updatePathForMultipleLinks_Node(pathID, index_p, index_s+1
				, running_optimalScore, running_removeThisPath, UP);
	
	//good 
	return ;
}

void LocalAlignment_CT::traceBack()
{
	//note the c_current_path_tb was set in tracebackMultiple to pass the information here.
	//
	
  //vector<unsigned int[2]> OptimalValueIndex; //the index that current best score resides
  //vector<unsigned int[2]> startIndex;//where the 
  //vector<double[2]> curr_and_max;

  //cout<<"doing trace back in the local alignment"<<endl;
  //flush(cout);
  //starting from optimalIndex going backing
  unsigned int i=c_optimalIndex[0];
  unsigned int j=c_optimalIndex[1];
  
  vector<unsigned int> vec_pathToBeChecked;//use this to remeber what to check after this path has been traced back;
  //vector<unsigned> vec_pathToBeChecked_StartIndex_Pattern;
  //vector<unsigned vec_pathToBeChecked_StartIndex_Subject;
#ifdef DEBUG
  cout<<"---------Start tracing back---------------"<<endl;
  cout<<"starting from ("<<i<<","<<j<<")"<<endl;
#endif
  string c_pattern_wg;//aligned string with gap
  string c_subject_wg;//aligend string with gap
  //cout<<"length of patten :"<<c_pattern->GetLength()<<endl;
  /*if(c_traceback_table[0+4*(c_pattern->GetLength()+1)].GetLinks()==ZERO)
    cout<<"======C_tracetable entry (0,4) is ZERO"<<endl;
  else
    cout<<"======C_tracetable entry (0,4) is NON-ZERO"<<endl;  
  */
  // while(c_traceback_table[i+j*(c_pattern->GetLength()+1)].GetLinks()!=ZERO)  //keep going till we reach a zero
  while(c_traceback_table->GetLink(i,j)!=ZERO)  //keep going till we reach a zero
    {
#ifdef DEBUG
      cout<<"\t("<<i<<","<<j<<"):"<<endl;
#endif
      unsigned int currentIndex;
      //check the link table, decide where to go
      //switch(c_traceback_table[i+j*(c_pattern->GetLength()+1)].GetLinks())
	  if(c_traceback_table->GetPathUsageState(i,j))
	  {
		  cout<<"ERROR:This one has been used by other path. something wrong!!"<<endl;
		  exit(-1);
	  }
      switch(c_traceback_table->GetLink(i,j))
		{
		case UP:
#ifdef DEBUG
			cout<<"UP:"<<endl;
#endif
		  //now we have to check how many indels
		  currentIndex=j;
		  for(unsigned int k =0;k<c_traceback_table->GetNumOfIndels(i, currentIndex);k++)//c_traceback_table[i+currentIndex*(c_pattern->GetLength()+1)].GetNumOfIndels();k++)
			{
			  c_pattern_wg="-"+c_pattern_wg;
		  
			  c_subject_wg=c_subject->GetSequence().at(j-1)+c_subject_wg;
			  c_traceback_table->SetPathUsageState(i,j,true);//this one has been used.
			  j--;
			}
		  break;
		case LEFT:
#ifdef DEBUG
		  cout<<"LEFT"<<endl;
#endif
		  //need to check how many indels
		  currentIndex=i;
		  for(unsigned int k=0;k<c_traceback_table->GetNumOfIndels_s(currentIndex,j);k++)//c_traceback_table[currentIndex+j*(c_pattern->GetLength()+1)].GetNumOfIndels();k++)
			{
			  c_pattern_wg=c_pattern->GetSequence().at(i-1)+c_pattern_wg;
			  c_subject_wg="-"+c_subject_wg;
			  c_traceback_table->SetPathUsageState(i,j,true);//this one has been used.
			  i--;
			}
		  break;
		case UPLEFT:
#ifdef DEBUG
		  cout<<"UPLEFT"<<endl;
#endif
		  c_pattern_wg=c_pattern->GetSequence().at(i-1)+c_pattern_wg;
		  c_subject_wg=c_subject->GetSequence().at(j-1)+c_subject_wg;
		  c_traceback_table->SetPathUsageState(i,j,true);//this one has been used.
		  i--;
		  j--;
		  break;
		case UP_UPLEFT: 
			//we take upleft first, and also we have to check for the path 
#ifdef DEBUG
		  cout<<"UP_UPLEFT"<<endl;
#endif
		  recordPathInfoToBeRemove(i,j, c_current_path_tk, vec_pathToBeChecked 
					);
		  if(isPathThruThisNode(c_current_path_tk, i-1, j-1))
		  {
#ifdef DEBUG
			  cout<<"\taking  UPLEFT route..."<<endl;
#endif
			  c_pattern_wg=c_pattern->GetSequence().at(i-1)+c_pattern_wg;
			  c_subject_wg=c_subject->GetSequence().at(j-1)+c_subject_wg;
			  c_traceback_table->SetPathUsageState(i,j,true);//this one has been used.
			  
			  i--;
			  j--;
		  }
		  else //take the up path, this should be fine, it must have path through this 
		  {
#ifdef DEBUG
			  cout<<"\taking  UP route"<<endl;
#endif              
			currentIndex=j;
			for(unsigned int k =0;k<c_traceback_table->GetNumOfIndels(i, currentIndex);k++)//c_traceback_table[i+currentIndex*(c_pattern->GetLength()+1)].GetNumOfIndels();k++)
			{
			  c_pattern_wg="-"+c_pattern_wg;
		  
			  c_subject_wg=c_subject->GetSequence().at(j-1)+c_subject_wg;
			  c_traceback_table->SetPathUsageState(i,j,true);//this one has been used.
			  j--;
			}
		  }
		  break;
		case LEFT_UPLEFT: 
#ifdef DEBUG
		  cout<<"LEFT_UPLEFT"<<endl;
#endif
			recordPathInfoToBeRemove(i,j, c_current_path_tk, vec_pathToBeChecked 
					);
			if(isPathThruThisNode(c_current_path_tk, i-1, j-1))
			{
			//take upleft first, if it allows
#ifdef DEBUG
			  cout<<"\t taking UPLEFT route"<<endl;
#endif
			  c_pattern_wg=c_pattern->GetSequence().at(i-1)+c_pattern_wg;
			  c_subject_wg=c_subject->GetSequence().at(j-1)+c_subject_wg;
			  c_traceback_table->SetPathUsageState(i,j,true);//this one has been used.
			  
			  i--;
			  j--;
			}
			else //take the left path, this should be fine, it must have a path through this
			{
#ifdef DEBUG
			  cout<<"\t taking LEFT route"<<endl;
#endif
				currentIndex=i;
			  for(unsigned int k=0;k<c_traceback_table->GetNumOfIndels_s(currentIndex,j);k++)//c_traceback_table[currentIndex+j*(c_pattern->GetLength()+1)].GetNumOfIndels();k++)
				{
				  c_pattern_wg=c_pattern->GetSequence().at(i-1)+c_pattern_wg;
				  c_subject_wg="-"+c_subject_wg;
				  c_traceback_table->SetPathUsageState(i,j,true);//this one has been used.
				  i--;
				}
			}
		  break;
		case UP_LEFT: //in this case, we take left link nodes first (subject gap), 
#ifdef DEBUG
		  cout<<"UP_LEFT"<<endl;
#endif
		  //need to check how many indels, we need to remeber each one 
		  //if they are mutiple 
		  
		  recordPathInfoToBeRemove(i,j, c_current_path_tk, vec_pathToBeChecked
				);
		 //take the left first, if it is possible
		 if(isPathThruThisNode(c_current_path_tk, i, j-1))
		 {
#ifdef DEBUG
			  cout<<"\t taking LEFT route"<<endl;
#endif             
			 currentIndex=i;
		  for(unsigned int k=0;k<c_traceback_table->GetNumOfIndels_s(currentIndex,j);k++)//c_traceback_table[currentIndex+j*(c_pattern->GetLength()+1)].GetNumOfIndels();k++)
			{
			  c_pattern_wg=c_pattern->GetSequence().at(i-1)+c_pattern_wg;
			  c_subject_wg="-"+c_subject_wg;
			  c_traceback_table->SetPathUsageState(i,j,true);//this one has been used.
			  i--;
			}
		 }
		 else //so take the up then 
		 {
#ifdef DEBUG
			  cout<<"\t taking UP route"<<endl;
#endif
			 currentIndex=j;
			for(unsigned int k =0;k<c_traceback_table->GetNumOfIndels(i, currentIndex);k++)//c_traceback_table[i+currentIndex*(c_pattern->GetLength()+1)].GetNumOfIndels();k++)
			{
			  c_pattern_wg="-"+c_pattern_wg;
		  
			  c_subject_wg=c_subject->GetSequence().at(j-1)+c_subject_wg;
			  c_traceback_table->SetPathUsageState(i,j,true);//this one has been used.
			  j--;
			}
		 }
		  break;
		case UP_LEFT_UPLEFT:
#ifdef DEBUG
			  cout<<"UP_LEFT_UPLEFT"<<endl;
#endif
			recordPathInfoToBeRemove(i,j, c_current_path_tk, vec_pathToBeChecked 
					);
			if(isPathThruThisNode(c_current_path_tk, i-1, j-1))
			{
			//take upleft first, if it allows
#ifdef DEBUG
			  cout<<"taking UPLEFT route..."<<endl;
#endif
			  c_pattern_wg=c_pattern->GetSequence().at(i-1)+c_pattern_wg;
			  c_subject_wg=c_subject->GetSequence().at(j-1)+c_subject_wg;
			  c_traceback_table->SetPathUsageState(i,j,true);//this one has been used.
			  
			  i--;
			  j--;
			}
			else //take the left path, this should be fine, it must have a path through this
			{
				if(isPathThruThisNode(c_current_path_tk, i, j-1))
				{
#ifdef DEBUG
			  cout<<"taking LEFT route..."<<endl;
#endif                    
					currentIndex=i;
				  for(unsigned int k=0;k<c_traceback_table->GetNumOfIndels_s(currentIndex,j);k++)//c_traceback_table[currentIndex+j*(c_pattern->GetLength()+1)].GetNumOfIndels();k++)
					{
					  c_pattern_wg=c_pattern->GetSequence().at(i-1)+c_pattern_wg;
					  c_subject_wg="-"+c_subject_wg;
					  c_traceback_table->SetPathUsageState(i,j,true);//this one has been used.
					  i--;
					}
				}
				else
				{
#ifdef DEBUG
			  cout<<"taking UP route.........."<<endl;
#endif                    
					currentIndex=j;
					for(unsigned int k =0;k<c_traceback_table->GetNumOfIndels(i, currentIndex);k++)//c_traceback_table[i+currentIndex*(c_pattern->GetLength()+1)].GetNumOfIndels();k++)
					{
					  c_pattern_wg="-"+c_pattern_wg;
				  
					  c_subject_wg=c_subject->GetSequence().at(j-1)+c_subject_wg;
					  c_traceback_table->SetPathUsageState(i,j,true);//this one has been used.
					  j--;
					}
				}
			}
		  break;
		default:
		  cerr<<"ERROR:not defined entry in trace back table!! Exit!"<<endl;
		  exit(-1);
		  break;
		}//switch
    }  //end of while loop.
#ifdef DEBUG
	cout<<"Done with tracing for the this current path ready to update !!!"<<endl;
#endif
  //we are done
  
  c_alignment.SetPattern(c_pattern_wg, true);
  c_alignment.SetPattern(c_pattern->GetSequence().substr(i,c_optimalIndex[0]-i), false);
  
  c_alignment.SetSubject(c_subject_wg, true);
  c_alignment.SetSubject(c_subject->GetSequence().substr(j,c_optimalIndex[1]-j), false);
  
  c_alignment.SetPatternIndex(i, c_optimalIndex[0]-1);
  c_alignment.SetSubjectIndex(j,c_optimalIndex[1]-1);
  c_alignment.SetScore(c_score);
  
  //done with trace back, set it to be traced
  	c_path_vec.at(c_current_path_tk)->setTraced();
#ifdef DEBUG
	cout<<"before doing the update for mutltiple links ( recurisive one)"<<endl;
#endif
  //now we are done for this path track back,
  //no we need to go through the path list to update in case some node has been take by the current path
  updatePathForMultipleLinks(c_current_path_tk, vec_pathToBeChecked);
#ifdef DEBUG			
  cout<<"+++++++++done the trace back "<<endl;
	//flush(cout);
#endif 
}
//internal function (protected)
//used by --the align()-- function for the case where we have linkback to UPLEFT ONLY case 
//, to update/write the path information to the pathelemententry 
//nodes, and also create the new path and add to the path vector of the LocalAlignment_CT object
//i is the pattern and j is the subject , current score.

void LocalAlignment_CT::updateNodePath(const unsigned& i, const unsigned& j, const unsigned& lenP, 
		const unsigned& lenS, const double & compval)
{
	if(c_traceback_table->GetLink(i-1,j-1)==ZERO)//[(i-1)+(j-1)*(lenP+1)].GetLinks()==ZERO)
	{
	  //cout<<"\t\t&&&&&&&&&&&&ADDING A NEW PATH"<<endl;
	  //this is new path starting. since it following a zero and then a match.
	  //we need to add a entry to a path vec to keep track this new path
	  unsigned int tempIndexArr[2];tempIndexArr[0]=i;tempIndexArr[1]=j;
	  c_path_vec.push_back(new Path(tempIndexArr, tempIndexArr,compval ));
	  c_pathID_vec.push_back(c_pathID_vec.size());
	  
	  //set up path element info for the node, a new path just started,
	  //add without checking since this is a new path and there is no way we have added to 
	  //the pathelementable.
		c_PathElementTable[i+j*(lenP+1)].AddPathInfoNoCheckingDuplication(
				c_PathElementTable[i+j*(lenP+1)].GetNumberOfPathes(),c_path_vec.size()-1,
				compval,i,j);
	}
   else //not a zero, it was a larger than zero one, it has to be on a previous path
	{
	  //first get the path from the vector
	  //depends on whether there are many different path in the upleft node, 
	  //we need to check every path for the optimal value and depends
	  //on whether the optimal value is larger and update the optimal value and indexs
	  //cout<<"in calling function"<<endl;
	  Path* tempPath;
	  for(unsigned k=0;k<c_PathElementTable[i-1+(j-1)*(lenP+1)].GetNumberOfPathes();k++)
	  {
		  //cout<<"\t----k:"<<k<<endl;
		unsigned pathID  =c_PathElementTable[i-1+(j-1)*(lenP+1)].GetPathID(k); 
		  tempPath=c_path_vec.at(pathID);
		  PathElementEntry npe;
		  
		  unsigned int tempIndexArr[2];
		  const unsigned * temp=tempPath->GetOptimalIndex();
		  tempIndexArr[0]=temp[0];tempIndexArr[1]=temp[1];
		  double pathOV=tempPath->GetOptimalValue();
		  //cout<<"\t\tcalling checking"<<endl;
		  //now compare to update or not
		  if(compval>=tempPath->GetOptimalValue())
			{
			  //cout<<"\t\t\t***********ffound a maximum one for the existing path, update info::path index"<<pathElementTable[i-1+(j-1)*(lenP+1)]<<endl;
			  tempIndexArr[0]=i;tempIndexArr[1]=j;
			  tempPath->SetOptimalIndex(tempIndexArr);
			  tempPath->SetOptimalValue(compval);
			  pathOV=compval;
			  //tempPath->
			}
			//cout<<"===done"<<endl;
			//here we create a new PathElementEntry and add it to the current node (checking duplicate, in case
			//the path has been added by previous step.?? is that possible. might not be possible.
			//here the reason to create new node, because this is where we add new path information to
			//the PathElementEntry table. 
		  npe.AddPathInfoNoCheckingDuplication(0, pathID, pathOV, 
				tempIndexArr[0], tempIndexArr[1]);
		   //cout<<"\t\t\tin between"<<endl;cout<<"i:"<<i<<";j:"<<j<<"compval:"<<compval<<endl;
		  c_PathElementTable[i+j*(lenP+1)].AddPathInfo(npe);
		  //cout<<"endl"<<endl;
		  //here,we overwrote the information to the k-th nodes with "new" information.
		  //although the information might not be new.
		  //c_PathElementTable[i+j*(lenP+1)].AddPathInfoNoCheckingDuplication(k, pathID, pathOV, 
		  //		tempIndexArr[0], tempIndexArr[1]);
	  }
	}
	
}

//the variant of alignment using less memeory
//this one we don't keep the original mxn dp table,
//but instead we keep only one column (in fact for coding
//purpose, we keep two columns to make the algorithm working
//efficient;

//*********NOTE:
//     patter is on the top is the columns
//	subject is on the left is the rows
//			patter p p p............p
//			S
//			S 
//			S 
//			.
//			.
//			.
//			s

void LocalAlignment_CT::align()
{
#ifdef DEBUG
	cout<<"+++++++++doing the alignment "<<endl;
	//flush(cout);
#endif
  	//Create empty dynamic programming table
  unsigned lenP=c_pattern->GetLength();
  unsigned lenS=c_subject->GetLength();
  double* dp_table_prev_col=new double[(lenS+1)];//one extra on this, to deal with the beging of the column
  double* dp_table_curr_col=new double[(lenS+1)];//one extra on this, to deal with the beging of the column

  //GapModel* gm=new AffineGapModel(c_gapOpen, c_gapExtension);
  //GapModel* gm=new MarkovChainGapModel_454(c_gapOpen, c_gapExtension, c_pattern->GetSequence(), c_subject->GetSequence());
  //here we keep two columns, to effieciently and easily manipulate the column. it cold be only one
  c_traceback_table=new TracebackTable(c_pattern, c_subject, 2);//(lenP+1)*(lenS+1)];
  //the following is not necessary, so we are defaulting all to zero
  /*for(unsigned int i=0;i<=lenP;i++)
    {
      c_traceback_table[i+0*(lenP+1)].SetLinks(ZERO);
    }
  for(unsigned int j=0;j<=lenS;j++)
    {
      c_traceback_table[0+j*(lenP+1)].SetLinks(ZERO);
    }
  */
  //cout<<"lenP is "<<lenP<<endl;
  /*if(c_traceback_table[0+4*(lenP+1)].GetLinks()==ZERO)
    cout<<"======C_tracetable entry (0,4) is ZERO"<<endl;
  else
    cout<<"======C_tracetable entry (0,4) is NON-ZERO"<<endl;
  */
  c_score=-1E9;
  
  //intialize the dp table, the first row=0 and first column=0
  //the row is Pattern, column is Subject. 
  dp_table_curr_col[0]=0;//this is the first row

  for(unsigned int i=0;i<=lenS;++i)
    {
      dp_table_prev_col[i]=0;//this is the first column
    }

  //cout<<"successfully created the empty tables and now go to get the score!\n";

  //score table with S-W
  double compval = 0;
  string strP=c_pattern->GetSequence();
  string strS=c_subject->GetSequence();
  
  double* maximumGapValue_subject=new double[lenS+1]; 
	//we still keep this same len as the curr/prev col, but we will never use the first one.
  
  unsigned int* maximumGapIndex_subject=new unsigned int[lenS+1];
	//just as above, this one is used to keep record of the maximum Gap Value so far, and the first one [0] is not used.
  
  //maximumGapValue[0]=-1E9;
  //maximumGapValue[1]=c_gapOpen+c_gapExtension*1;
  //we will start doing the job at column 1, so set intial value of gapIndex=0 to infinity
  for(unsigned int i=0;i<=lenS;i++)
    {
      maximumGapValue_subject[i]=-1E50;
      maximumGapIndex_subject[i]=-0;
      //  maximumGapIndex[i]=0;
      //maximumGapIndex[1]=1;//to the
    } 

  //for pattern gap, this maximumGapValue, a scalar instead of matrix is OK. 
  //since we keep this as a running one, like the current column
  double maximumGapValue_pattern=-1E60;
  unsigned int maximumGapIndex_pattern=0;

//update 1/19/2020 by add the pathElementTable as a class object, so we don't use it now
  //a new table used to keep track how each entry in a dp table or tracebacktable 
  //related to the path element in the vector path
//  unsigned int* pathElementTable=new unsigned int[(lenP+1)*(lenS+1)];   //<-- use class variable now.
  //again it has one extra one in each dimesnion, first ones are not used, only for a purpose of simplicity.
  //also not every elements in the talbe is assigned, ZERO link ones are not set, only the larger than zero one 
  //is set related to a new path
	//pathElementTable now is a class object, but we still initialize it here, same size, but so far
	//we don't have each individual one initialized, need to initialize at the moment to 
	//use it. 
	c_PathElementTable=new PathElementEntry [(lenP+1)*(lenS+1)];
  
  
  //we have to be very careful here
  //the dimension of the dp table and the strings are not same.
  //dp table is one row/col more than the strings, since we need to 
  //have the zero/starting row or column.so that means, when we dealing with
  //dp table, it is from 1->lenP or 1->lenS, but the string, starts from zero,
  //so we are doing things i-1
  for(unsigned int i = 1; i <= lenP; ++i) //for all values of strA
    {	//now we are starting a new column, we need to
      //1)first exchange the curr and prev col(actually, it is better to do this one at end of each loop, otherwise it will be trouble in the beginning of the first loop
      //2)reset the curr_col first one to zero, it should be zero anyway.
      
      dp_table_curr_col[0]=0; //local alignment so set to be zero 
      
      maximumGapValue_pattern=-1E60;
      maximumGapIndex_pattern=0;

      //now, we go through each element and do the job
      for(unsigned int j = 1; j <= lenS; ++j) 
	{	//for all values of strB
#ifdef DEBUG
	  cout<<"****doing round ("<<i<<","<<j<<")."<<endl;
#endif 	  
	  //MATCH
	  //if(strP.at(i-1) == strS.at(j-1)) 
	  //{				//if current sequence values are the same
	  //cout<<"\tcalling score matrix:score("<<strP.at(i-1)<<","<<strS.at(j-1)<<")="<<c_sm->GetScore(strP.at(i-1),strS.at(j-1))<<endl;
	  //cout<<"\t\tdp table is dp("<<i-1<<","<<j-1<<")="<<dp_table_prev_col[j-1]<<endl;

	  //***********do gap first*********
	  //here to make the affine linear model works, we need to keep a running max gap value for row across,
	  //since don't keep all the rows in the memory,
	  //but for the column across, we are fine, since we keep all the columns till this elements

	  //this is the one for all the sub rows, we need to keep a maximum one so far and update it with the information from this round.
	  //now we don't have to go through every entry, we only need to compare the Max one with this current one and keep track it.
	  //for(int k = i-1; k > 0; --k) 
	  // {		//check all sub rows
	  //cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&calling gm gapValue"<<endl;
	  //update the subject gap, row gap 
	  c_gm->GapValue(c_traceback_table, i, j, false, 
				dp_table_prev_col[j], maximumGapValue_subject[j], maximumGapIndex_subject[j]);	

	  /*double openNewGapValue=dp_table_prev_col[j]+c_gapOpen + c_gapExtension;
	  double maxGapExtendedValue=maximumGapValue[j]+c_gapExtension;
	  //check to update
	  if(openNewGapValue>=maxGapExtendedValue)
	    {
	      maximumGapValue[j]=openNewGapValue;
	      maximumGapIndex[j]=i-1;
	    }
	  else
	    {
	      maximumGapValue[j]=maxGapExtendedValue;
	      //the maximumGapInde[j] unchaned, keep the same
	    }
	  */
	  //if(compval < maximumGapValue[j]) 
	  //  {	    //if cell above has a greater value 
	      
	  compval = maximumGapValue_subject[j];		//set compval to that square
	  //c_traceback_table[i+j*(lenP+1)].SetLinks(LEFT);
	  //c_traceback_table[i+j*(lenP+1)].SetNumOfIndels(i-maximumGapIndex[j]);
	  //assume now this is the best score and set it as if this is a true good one, we will update this in the next two steps
	  //very specifically we need to reset the pattern gap to be zero, if this is a true best score. but we can not simply do this now. it will
	  //interfer the later steps.
	  c_traceback_table->SetTableEntry(i, j, LEFT, c_traceback_table->GetNumOfIndels(i,j),i-maximumGapIndex_subject[j]);
	  //  }
	  //cout<<"maximumGapValue[j]:"<<maximumGapValue_subject[j]<<",";
	  //cout<<"campval after rowGap:"<<compval<<";";

	  //do pattern gap now.
	  //this is the pattern/column gap, we keep it same as we are doing with the whole dp table
	  c_gm->GapValue(c_traceback_table, i,j,true, dp_table_curr_col[j-1], maximumGapValue_pattern, maximumGapIndex_pattern);

	  /*	  //this is the column across, we keep it same as we are doing with the whole dp table
	  for(int k=j-1; k>0; --k) 
	    {		//check all sub columns
				
	      if(compval < ((dp_table_curr_col[k]) + (c_gapOpen + (c_gapExtension *(j- k))))) 
		{	
		  //if square to the left has the highest valu
					
		  compval = ((dp_table_curr_col[k]) + (c_gapOpen + (c_gapExtension *(j- k))));    //set compval to that square
		  //c_traceback_table[i+j*(lenP+1)].SetLinks(UP);
		  //c_traceback_table[i+j*(lenP+1)].SetNumOfIndels(j-k);
		  c_traceback_table->SetTableEntry(i,j, UP, j-k);
		}
		}*/	
	  if(
			(compval-maximumGapValue_pattern) < -1E-12 //this gap pattern is bigger, allow a small error 
		)
	    {
	      //set compval to that square, UP only
	      compval=maximumGapValue_pattern;
	      c_traceback_table->SetTableEntry(i,j,UP, j-maximumGapIndex_pattern,0);//now this is a pattern gap, so set the subject pattern to be zero 
	    }
	  else
	  {
		  if(abs(compval - maximumGapValue_pattern)<1E-12) //special case two gaps EQUAL, set both route UP and LEFT,
						//note, in here we assume the subject gap of NumOfIndels is set correctly, so we simply read it and reset it accordingly
		  {
			  //
			  c_traceback_table->SetTableEntry(i,j,UP_LEFT, 
						 j-maximumGapIndex_pattern,c_traceback_table->GetNumOfIndels_s(i,j));

		  }
		  else  //the left gap value is the largest so far. set it back to LEFT, LEFT only!!
				//we do this because we want to reset pattern gap to be zero.<-----, do we have to do this????
		  {
			  c_traceback_table->SetTableEntry(i, j, LEFT, 0,i-maximumGapIndex_subject[j]);
		  }
	  }

	  //cout<<"maximumGapValue[j]:"<<maximumGapValue_subject[j]<<",";
	  //cout<<"campval after rowGap:"<<compval<<";";
	  //cout<<"compval afer col gap:"<<compval<<endl;	

	  //***************do Match/Mismatch now***********
	  if( //match and mismatch from upleft is bigger. so we need take this path, score of (i-1),(j-1), because the string is starting at zero and dp table is starting at 1
			(compval - (dp_table_prev_col[(j-1)] + c_sm->GetScore(strP.at(i-1),strS.at(j-1)))) < -1E-12
	    )
	    {
	      compval = (dp_table_prev_col[(j-1)] + c_sm->GetScore(strP.at(i-1),strS.at(j-1)));	//compval = diagonal + match score
	      //}
	      //cout<<"\t\tcompval after match/mismatch:"<<compval<<";";
	      c_traceback_table->SetTableEntry(i,j, UPLEFT,0,0);//c_traceback_table[i+j*(lenP+1)].SetLinks(UPLEFT);//for this one, we don't have to set the #numOfIndels, since there is none
	    }
	  else //in this case we have a tie, meaning we can go up and upleft or go left and upleft. this only(?)happens when
						//we have a gap value and a mismatch penalty equals.
	  {
		 if(abs(compval-(dp_table_prev_col[(j-1)] + c_sm->GetScore(strP.at(i-1),strS.at(j-1))))
					< 1E-12
			) //equal case, this is subtle, we currently only set it up
		 {
			 if(c_traceback_table->GetLink(i,j)==UP)
			 {
				c_traceback_table->SetTableEntry(i,j, UP_UPLEFT, 
							c_traceback_table->GetNumOfIndels(i,j),0);
			 }
			 if(c_traceback_table->GetLink(i,j)==LEFT)
			 {
				 
				c_traceback_table->SetTableEntry(i,j, LEFT_UPLEFT, 
							0,c_traceback_table->GetNumOfIndels_s(i,j));
			 }
			 if(c_traceback_table->GetLink(i,j)==UP_LEFT)
			 {
				 //
				 c_traceback_table->SetTableEntry(i,j, UP_LEFT_UPLEFT, 
							c_traceback_table->GetNumOfIndels(i,j),c_traceback_table->GetNumOfIndels_s(i,j));
			 }
		 } 
		 //else //do nothing, in this case, the compval > alignment score. we will check whether compval<0 next 
	  }
	  
	  //here, we try to update information after compare values 	
	  if((compval-0) < 1E-15 ) //found a zero or smaller value
	    {
	      compval = 0;
	      //The following is not necessary, since everything so far is default to ZERO
	      //c_traceback_table[i+j*(lenP+1)]=ZERO;
	      c_traceback_table->SetTableEntry(i,j, ZERO,0);
	      //we don't set numOfIndels;keep default
	    }
	  else //found a larger than zero value (no value could , this is where we need to check whether this is good starting for a new path
	    {
			//cout<<"Before+++:";
			//cout<<c_PathElementTable[i+j*(lenP+1)].toString()<<endl;

			LinkBack lb=c_traceback_table->GetLink(i,j);
			switch(lb)
			{
				case UPLEFT:
				//cout<<"\t|| --link: UPLEFT || ";fflush(stdout);
				//it is possible that this is a new path start. 
				//there is no way that a path starting from an insertion, it has to be a match
				this->updateNodePath(i,j,lenP, lenS, compval);
				/*	if(c_traceback_table->GetLink(i-1,j-1)==ZERO)//[(i-1)+(j-1)*(lenP+1)].GetLinks()==ZERO)
						{
						  //cout<<"\t\t&&&&&&&&&&&&ADDING A NEW PATH"<<endl;
						  //this is new path starting. since it following a zero and then a match.
						  //we need to add a entry to a path vec to keep track this new path
						  unsigned int tempIndexArr[2];tempIndexArr[0]=i;tempIndexArr[1]=j;
						  c_path_vec.push_back(new Path(tempIndexArr, tempIndexArr,compval ));
						  //set up path element info for the node, a new path just started,
						c_PathElementTable[i+j*(lenP+1)].AddPathInfoNoCheckingDuplication(
									c_PathElementTable[i+j*(lenP+1)].GetNumberOfPathes(),c_path_vec.size()-1,
									compval,i,j);
						}
					  else //not a zero, it was a larger than zero one, it has to be on a path
						{
						  //first get the path from the vector
						  //depends on whether there are many different path in the upleft node, 
						  //we need to check every path inform.
						  
						  Path* tempPath;
						  for(unsigned i=0;i<c_PathElementTable[i-1+(j-1)*(lenP+1)].GetNumberOfPathes();i++)
						  {
							unsigned pathID  =c_PathElementTable[i-1+(j-1)*(lenP+1)].GetPathID(i); 
							  tempPath=c_path_vec.at(pathID);
							  PathElementEntry npe;
							  
							  unsigned int tempIndexArr[2];
							  const unsigned * temp=tempPath->GetOptimalIndex();
							  tempIndexArr[0]=temp[0];tempIndexArr[1]=temp[1];
							  double pathOV=tempPath->GetOptimalValue();
							  //now compare to update or not
							  if(compval>tempPath->GetOptimalValue())
								{
								  //cout<<"\t\t\t***********ffound a maximum one for the existing path, update info::path index"<<pathElementTable[i-1+(j-1)*(lenP+1)]<<endl;
								  tempIndexArr[0]=i;tempIndexArr[1]=j;
								  tempPath->SetOptimalIndex(tempIndexArr);
								  tempPath->SetOptimalValue(compval);
								  pathOV=compval;
								  //tempPath->
								}
								
							  npe.AddPathInfoNoCheckingDuplication(0, pathID, pathOV, 
									tempIndexArr[0], tempIndexArr[1]);
							   
							  c_PathElementTable[i+j*(lenP+1)].AddPathInfo(npe);
						  }
						}*/
					break;
				case UP:
				//cout<<"\t|| --link: UP || ";fflush(stdout);
				//for not a UPLEFT entry, means this one is a indel, the score is descreasing, well so 1)it has to follow 
				//a path 2)it is no way to the maximum one so far, we only need to update the pathElement table information, don't have
				//update the path entry
					c_PathElementTable[i+j*(lenP+1)].AddPathInfo(c_PathElementTable[i+(j-c_traceback_table->GetNumOfIndels(i,j))*(lenP+1)]);
					break;
				case LEFT:
				//cout<<"\t|| --link: LEFT || ";fflush(stdout);
				//for not a UPLEFT entry, means this one is a indel, the score is descreasing, well so 1)it has to follow 
				//a path 2)it is no way to the maximum one so far, we only need to update the pathElement table information, don't have
				//update the path entry
					c_PathElementTable[i+j*(lenP+1)].AddPathInfo(c_PathElementTable[(i-c_traceback_table->GetNumOfIndels_s(i,j))+j*(lenP+1)]);	//pathElementTable[i+j*(lenP+1)]=pathElementTable[(i-c_traceback_table->GetNumOfIndels(i,j))+j*(lenP+1)];
					break;
				case UP_UPLEFT:
				//cout<<"\t|| --link: UP_UPLEFT || ";fflush(stdout);
				//in this case we have the entry pointing back with two path up and upLEFT, meaning from up by pattern gap and
				//from upleft by mismatch lead to identical score, so we keep both and update both path 
				//this could happen when the mismatch and gap penalty score are identical.
				   //doing up link first
				   //cout<<"i:"<<i<<";j:"<<j<<"compval:"<<compval<<endl;
					//cout<<"num of inde:"<<c_traceback_table->GetNumOfIndels_s(i,j)<<endl;
				   c_PathElementTable[i+j*(lenP+1)].AddPathInfo(c_PathElementTable[i+(j-c_traceback_table->GetNumOfIndels(i,j))*(lenP+1)]);
				//take care of the UPLEFT case
				//cout<<"in between"<<endl;
					this->updateNodePath(i,j,lenP, lenS, compval);
					//cout<<"done for this one"<<endl;
					break;
				case LEFT_UPLEFT:
				//cout<<"\t|| --link: LEFT_UPLEFT || ";fflush(stdout);
				//in this case we have the entry pointing back with two path left and upLEFT, meaning from up by subject gap and
				//from upleft by mismatch lead to identical score, so we keep both and update both path 
				//this could happen when the mismatch and gap penalty score are identical.
					//add path information to node for UP link
					
					c_PathElementTable[i+j*(lenP+1)].AddPathInfo(c_PathElementTable[(i-c_traceback_table->GetNumOfIndels_s(i,j))+j*(lenP+1)]);
					//cout<<"i:"<<i<<";j:"<<j<<endl;
					//adding path information to node for upleft
					this->updateNodePath(i,j,lenP, lenS, compval);
					
					break;
				case UP_LEFT:
				//cout<<"\t|| --link: UP_&_LEFT || ";fflush(stdout);
				//in this special case we have the entry pointing back with two linkbacks left and up, meaning from up by subject gap and
				//from left by pattern gap lead to identical score and also identical paths??? (possible, but
				//not 100% sure), so we keep both linkback and update path (also check to make sure, we did not 
				//duplicate things. 
				//this could happen when the mismatch and gap penalty score are identical.
				//we basically add both up and left node to this current one
				//this simply 
					//add path information to node for UP link
					c_PathElementTable[i+j*(lenP+1)].AddPathInfo(c_PathElementTable[i+(j-c_traceback_table->GetNumOfIndels(i,j))*(lenP+1)]);
					//add path information to node for UP link
					c_PathElementTable[i+j*(lenP+1)].AddPathInfo(c_PathElementTable[(i-c_traceback_table->GetNumOfIndels_s(i,j))+j*(lenP+1)]);
					break;
				case UP_LEFT_UPLEFT:
				//cout<<"\t|| --link: UP_&_LEFT_&_UPLEFT || ";fflush(stdout);
				//in this special case we have the entry pointing back with THREE linkbacks left and up, meaning from up by subject gap and
				//from left by pattern gap lead to identical score and also upleft by mismatch
				//identical paths ??(possible, but
				//not 100% sure), so we keep all linkback and update path (also check to make sure, we did not 
				//duplicate things. 
				//this could happen when the mismatch and gap penalty score are identical.
				//we basically add both up and left node to this current one
				//this simply 
					//add path information to node for UP link
					c_PathElementTable[i+j*(lenP+1)].AddPathInfo(c_PathElementTable[i+(j-c_traceback_table->GetNumOfIndels(i,j))*(lenP+1)]);
					//add path information to node for UP link
					c_PathElementTable[i+j*(lenP+1)].AddPathInfo(c_PathElementTable[(i-c_traceback_table->GetNumOfIndels_s(i,j))+j*(lenP+1)]);
					//adding path information to node for upleft
					this->updateNodePath(i,j,lenP, lenS, compval);
					break;
				default :
					cerr<<"ERROR for setting up the link back entry"<<endl; 
					exit(-2);
					break;
			}

/*				  //check its
				if(c_traceback_table->GetLink(i,j)==UPLEFT)//[i+j*(lenP+1)].GetLinks()==UPLEFT)
				{
				  
				}
				else  //for not a UPLEFT entry, means this one is a indel, the score is descreasing, well so 1)it has to follow 
				//a path 2)it is no way to the maximum one so far, we only need to update the pathElement table information, don't have
				//update the path entry
				{
				  LinkBack tempLink=c_traceback_table->GetLink(i,j);//[i+j*(lenP+1)].GetLinks();//it can not be ZERO, can not be UPPERLEFT
				  
				  if(tempLink==UP)
					{
					  //pathElementTable[i+j*(lenP+1)]=pathElementTable[i+(j-c_traceback_table[i+j*(lenP+1)].GetNumOfIndels())*(lenP+1)];
					  
					}
				  else
					{
					  if(tempLink==LEFT)
						{
						  //pathElementTable[i+j*(lenP+1)]=pathElementTable[(i-c_traceback_table[i+j*(lenP+1)].GetNumOfIndels())+j*(lenP+1)];
						  
						}
					}
				}
			}*/
		}
	  //cout<<"This current element is from path index:"<<pathElementTable[i+j*(lenP+1)];
	  /*
	  switch (c_traceback_table->GetLink(i,j))
	    {
	    case UP:
	      cout<<"\t\tlink is UP;#indels is"<<c_traceback_table->GetNumOfIndels(i,j)<<endl;
	      break;
	    case LEFT:
	      cout<<"\t\tlink is LEFT;#indels is"<<c_traceback_table->GetNumOfIndels(i,j)<<endl;
	      break;
	    case UPLEFT:
	      cout<<"\t\tlink is UPLEFT"<<endl;
	      break;
	    case ZERO:
	      cout<<"\t\tlink is ZERO"<<endl;
	      break;
	    default:
	      cout<<"\t\tlink is not found"<<endl;
	      break;
	      
	    }
	  */
	  dp_table_curr_col[j] = compval;	//set current cell to highest possible score and move on
#ifdef DEBUG
	  cout<<"\t\trunning score:"<<compval<<endl;
#endif 
	  //find a best one so far, acutally these values here are not very useful in LocalAlignment_CT, 
	  //they are usful for doing pairwise and overlap, since there are only one best value and, 
	  //but in local alignment we recording many best vlues for each path.
	  //we still keep track of them anyway.
	  if(c_score<compval)
	    {
	      c_score=compval;
	      c_optimalIndex[0]=i;
	      c_optimalIndex[1]=j;
	    }
#ifdef DEBUG
	cout<<"\t\tafter+++:";
					cout<<c_PathElementTable[i+j*(lenP+1)].toString()<<endl;
					cout<<"\tlink:"<<c_traceback_table->GetLink(i,j)<<endl;
#endif
	}//end of col
      
      //reset the prev one to curr col and make it ready for next loop
      double* tempP=dp_table_curr_col;
      dp_table_curr_col=dp_table_prev_col;
      dp_table_prev_col=tempP;
    }//end of row
	//cout<<"end of the align()"<<endl;
	//printPathEntryTable();
	//cout<<"---total number of pathes:"<<c_pathID_vec.size()<<endl;
  //clean up
  delete[] dp_table_prev_col;
  delete[] dp_table_curr_col;
  delete[] maximumGapValue_subject;
  delete[] maximumGapIndex_subject;
  
  //cout<<"+++++++++done with the alignment "<<endl;
  //	flush(cout);
  //delete gm;
  //traceback_table will be deleted upon destruction
}//end of the localAlign

void LocalAlignment_CT::printPathEntryTable() const
{
	//printPathEntryTable

	unsigned lenP=this->c_pattern->GetLength();
	unsigned lenS=this->c_subject->GetLength();
	cout<<"----PathEntryTable("<<lenP<<"x"<<lenS<<")----"<<endl;
	PathElementEntry* p_pee;
	for(unsigned i=1;i<lenP+1;i++)
	{
		for(unsigned j=1;j<lenS+1;j++)
		{
			
			cout<<"**Entry("<<i<<","<<j<<")--";
			p_pee=&(c_PathElementTable[i+j*(lenP+1)]);
			cout<<p_pee->toString();
		}
	}
	cout<<"-------------------------------"<<endl;
}


