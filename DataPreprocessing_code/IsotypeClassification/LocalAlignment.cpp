#include "LocalAlignment.hpp"
#include <stdlib.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>

#include "Path.hpp"
#include "GapModel.hpp"
#include "AffineGapModel.hpp"
#include "MarkovChainGapModel_454.hpp"

#define DEBUG_D

using namespace std;
static bool comparePathElement(Path* p1, Path* p2)
{
  return p1->GetOptimalValue() > p2->GetOptimalValue();
}

LocalAlignment::LocalAlignment(const SequenceString* _pattern, const SequenceString* _subject, 
			       const ScoreMatrix* _m, const double& _gopen, 
			       const double& _gextension, const double& _scale, const int& _numOfAlignments, const short& _typeOfGapModel):
  PairwiseAlignment(_pattern, _subject, _m, _gopen, _gextension, _scale,_typeOfGapModel),c_numOfAlignments(_numOfAlignments),
  c_alignmentArr(NULL), c_scoreArr(NULL)
{
  //now we need to
  c_alignmentArr=new AlignmentString[this->c_numOfAlignments];
  c_scoreArr=new double[this->c_numOfAlignments];
  //dp_table=NULL;
  //traceback_table=NULL;
  this->align();
  //cout<<"in between"<<endl;
  //flush(cout);
  this->traceBackMultiple();
}
  
LocalAlignment::~LocalAlignment()
{
  //the base class destructor is called automatically

  //here we only need to take care other issues
  if(c_alignmentArr!=NULL)
    delete[] c_alignmentArr;
  if(c_scoreArr!=NULL)
    delete[] c_scoreArr;
  //taking care of Path vector
  for(unsigned int i=0;i<c_path_vec.size();i++)
    {
      //cout<<"deleting path vecs"<<endl;
	  if(c_path_vec.at(i)!=NULL)
			delete c_path_vec.at(i);
    }
}

//this one now return a deep copy of score array. 
//the outer caller need to deinitialize memory  
double* LocalAlignment::GetScoreArr() const
{
	double * temp=new double [this->c_numOfAlignments];
	memcpy(temp, c_scoreArr, c_numOfAlignments*sizeof(double));
	
  return temp;
}

//this one return a pointer to the object internal alignment array . 
//use with caution.
double* LocalAlignment::GetScoreArr_p()
{
	double * temp=new double [this->c_numOfAlignments];
	memcpy(temp, c_scoreArr, c_numOfAlignments*sizeof(double));
	
  return temp;
}

//this one return a pointer to the object internal alignment array . 
//use with caution.
AlignmentString* LocalAlignment::GetAlignmentArr_p()
{

  return c_alignmentArr;
}

//this one now return a deep copy of score array. 
//the outer caller need to deinitialize memory  
unsigned int LocalAlignment::GetNumberOfAlignments() const
{
  return c_numOfAlignments;
}

//this one now return a deep copy of score array. 
//the outer caller need to deinitialize memory  
AlignmentString* LocalAlignment::GetAlignmentArr() const
{
	AlignmentString* temp=new AlignmentString[c_numOfAlignments];
	//memcpy(temp, this->c_alignmentArr, c_numOfAlignments*sizeof(AlignmentString));
	for(unsigned i=0;i<c_numOfAlignments;i++)
		temp[i]=this->c_alignmentArr[i];
  return temp;
}

/*
//this code is not working because of later changes on the variables and etc. but the idea is there.
void LocalAlignment::align() ----this one is working one, but keeping track of a full dp table.is replaced by the new low memory one.
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

void LocalAlignment::traceBackMultiple()
{
  //decide how many we will be tracing back 
  //unsigned int actualNumOfAlignments;
  //cout<<"+++++++++doing the multiple trce bck "<<endl;
  //cout<<"c_path_vec.size():"<<c_path_vec.size()<<endl;
	flush(cout);
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
  /*
  cout<<"******************SHOWING PATH INFO**************"<<endl;
  cout<<"\ttotal num of path:"<<c_path_vec.size()<<endl;
  for(unsigned int i=0;i<c_path_vec.size();i++)
    {
      cout<<i<<"/"<<c_path_vec.size()<<": start at ("<<c_path_vec.at(i)->GetStartIndex()[0]<<","
	  <<c_path_vec.at(i)->GetStartIndex()[1]<<") and end at ("<<c_path_vec.at(i)->GetOptimalIndex()[0]
	  <<","<<c_path_vec.at(i)->GetOptimalIndex()[1]<<");\n";
    }
  */
  //now sort the vector
  sort(c_path_vec.begin(), c_path_vec.end(), comparePathElement);
  
//cout<<"\tefore looping......."<<endl;
//flush(cout);

  for(unsigned int i=0;i<c_numOfAlignments;i++)
    {
	//cout<<"i:"<<i<<"..."<<endl;
      //prepare the input and output
      c_optimalIndex[0]=c_path_vec.at(i)->GetOptimalIndex()[0];
      c_optimalIndex[1]=c_path_vec.at(i)->GetOptimalIndex()[1];
      //cout<<"before doing trace back;"<<endl;
	  //flush(cout);
      //now call the traceBack to do the job
      traceBack();
      
      //now record the output of trace back
      c_alignmentArr[i].SetPattern(c_alignment.GetPattern(true), true);
      c_alignmentArr[i].SetPattern(c_alignment.GetPattern(false), false);
  
      c_alignmentArr[i].SetSubject(c_alignment.GetSubject(true), true);
      c_alignmentArr[i].SetSubject(c_alignment.GetSubject(false), false);
  
      c_alignmentArr[i].SetPatternIndex(c_alignment.GetPatternIndexStart(), c_alignment.GetPatternIndexEnd());
      c_alignmentArr[i].SetSubjectIndex(c_alignment.GetSubjectIndexStart(), c_alignment.GetSubjectIndexEnd());
      c_alignmentArr[i].SetScore(c_path_vec.at(i)->GetOptimalValue());
      c_scoreArr[i]=c_path_vec.at(i)->GetOptimalValue();
    }

  //after set up everything, we want to set the c_alignment to the best alignment one
c_alignment.SetPattern(c_alignmentArr[0].GetPattern(true), true);
      c_alignment.SetPattern(c_alignmentArr[0].GetPattern(false), false);
  
      c_alignment.SetSubject(c_alignmentArr[0].GetSubject(true), true);
      c_alignment.SetSubject(c_alignmentArr[0].GetSubject(false), false);
  
      c_alignment.SetPatternIndex(c_alignmentArr[0].GetPatternIndexStart(), c_alignmentArr[0].GetPatternIndexEnd());
      c_alignment.SetSubjectIndex(c_alignmentArr[0].GetSubjectIndexStart(), c_alignmentArr[0].GetSubjectIndexEnd());
      c_alignment.SetScore(c_path_vec.at(0)->GetOptimalValue());
      c_score=c_path_vec.at(0)->GetOptimalValue();
	  
	//  cout<<"+++++++++done the multiple trace back"<<endl;
	//flush(cout);
  
  //done!!!
}

void LocalAlignment::traceBack()
{
  //vector<unsigned int[2]> OptimalValueIndex; //the index that current best score resides
  //vector<unsigned int[2]> startIndex;//where the 
  //vector<double[2]> curr_and_max;

  //cout<<"doing trace back in the local alignment"<<endl;
  //flush(cout);
  //starting from optimalIndex going backing
  unsigned int i=c_optimalIndex[0];
  unsigned int j=c_optimalIndex[1];

  //cout<<"starting from ("<<i<<","<<j<<")"<<endl;
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
      //cout<<"\t("<<i<<","<<j<<"):";
      unsigned int currentIndex;
      //check the link table, decide where to go
      //switch(c_traceback_table[i+j*(c_pattern->GetLength()+1)].GetLinks())
      switch(c_traceback_table->GetLink(i,j))
	{
	case UP:
	  //cout<<"UP:"<<endl;
	  //now we have to check how many indels
	  currentIndex=j;
	  for(unsigned int k =0;k<c_traceback_table->GetNumOfIndels(i, currentIndex);k++)//c_traceback_table[i+currentIndex*(c_pattern->GetLength()+1)].GetNumOfIndels();k++)
	    {
	      c_pattern_wg="-"+c_pattern_wg;
	  
	      c_subject_wg=c_subject->GetSequence().at(j-1)+c_subject_wg;
	      j--;
	    }
	  break;
	case LEFT:
	  //cout<<"LEFT"<<endl;
	  //need to check how many indels
	  currentIndex=i;
	  for(unsigned int k=0;k<c_traceback_table->GetNumOfIndels_s(currentIndex,j);k++)//c_traceback_table[currentIndex+j*(c_pattern->GetLength()+1)].GetNumOfIndels();k++)
	    {

	      c_pattern_wg=c_pattern->GetSequence().at(i-1)+c_pattern_wg;
	      c_subject_wg="-"+c_subject_wg;
	      i--;
	    }
	  break;
	case UPLEFT:
	  //cout<<"UPLEFT"<<endl;
	  c_pattern_wg=c_pattern->GetSequence().at(i-1)+c_pattern_wg;
	  c_subject_wg=c_subject->GetSequence().at(j-1)+c_subject_wg;
	  i--;
	  j--;
	  break;
	default:
	  cerr<<"ERROR:not defined entry in trace back table!! Exit!"<<endl;
	  exit(-1);
	  break;
	}
    }
  //cout<<"Done with alingment!!!"<<endl;
  //we are done
  c_alignment.SetPattern(c_pattern_wg, true);
  c_alignment.SetPattern(c_pattern->GetSequence().substr(i,c_optimalIndex[0]-i), false);
  
  c_alignment.SetSubject(c_subject_wg, true);
  c_alignment.SetSubject(c_subject->GetSequence().substr(j,c_optimalIndex[1]-j), false);
  
  c_alignment.SetPatternIndex(i, c_optimalIndex[0]-1);
  c_alignment.SetSubjectIndex(j,c_optimalIndex[1]-1);
  c_alignment.SetScore(c_score);
  
  //cout<<"+++++++++done the trace back "<<endl;
	//flush(cout);
  
}


//the variant of alignment using less memeory
//this one we don't keep the original mxn dp table,
//but instead we keep only one column (in fact for coding
//purpose, we keep two columns to make the algorithm working
//efficient;
void LocalAlignment::align()
{
	//cout<<"+++++++++doing the alignment "<<endl;
	//flush(cout);
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

  //a new table used to keep track how each entry in a dp table or tracebacktable 
  //related to the path element in the vector path
  unsigned int* pathElementTable=new unsigned int[(lenP+1)*(lenS+1)];
  //again it has one extra one in each dimesnion, first ones are not used, only for a purpose of simplicity.
  //also not every elements in the talbe is assigned, ZERO link ones are not set, only the larger than zero one 
  //is set related to a new path

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
      
      
      dp_table_curr_col[0]=0;
      
      maximumGapValue_pattern=-1E60;
      maximumGapIndex_pattern=0;

      //now, we go through each element and do the job
      for(unsigned int j = 1; j <= lenS; ++j) 
	{	//for all values of strB
	  //cout<<"****doing round ("<<i<<","<<j<<")."<<endl;	
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
	  c_gm->GapValue(c_traceback_table, i,j,false, dp_table_prev_col[j], maximumGapValue_subject[j], maximumGapIndex_subject[j]);	

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
	  c_traceback_table->SetTableEntry(i, j, LEFT, c_traceback_table->GetNumOfIndels(i,j), i-maximumGapIndex_subject[j]);
	  //  }
	  //cout<<"maximumGapValue[j]:"<<maximumGapValue_subject[j]<<",";
	  //cout<<"campval after rowGap:"<<compval<<";";

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
	      //set compval to that square
	      compval=maximumGapValue_pattern;
	      c_traceback_table->SetTableEntry(i,j,UP, j-maximumGapIndex_pattern,0);
	    }	
		else
	  {
		  /*if(abs(compval - maximumGapValue_pattern)<1E-12) //special case two gaps EQUAL, set both route UP and LEFT,
						//note, in here we assume the subject gap of NumOfIndels is set correctly, so we simply read it and reset it accordingly
		  {
			  //
			  c_traceback_table->SetTableEntry(i,j,LEFT, 
						 j-maximumGapIndex_pattern,c_traceback_table->GetNumOfIndels_s(i,j));

		  }
		  else  //the left gap value is the largest so far. set it back to LEFT, LEFT only!!
				//we do this because we want to reset pattern gap to be zero.<-----, do we have to do this????
		  {*/
			  c_traceback_table->SetTableEntry(i, j, LEFT, 0,i-maximumGapIndex_subject[j]);
		  //}
	  }
	  //cout<<"maximumGapValue[j]:"<<maximumGapValue_subject[j]<<",";
	  //cout<<"campval after rowGap:"<<compval<<";";
	  //cout<<"compval afer col gap:"<<compval<<endl;	

	  //***************do Match/Mismatch now***********
	  if(compval - (dp_table_prev_col[(j-1)] + c_sm->GetScore(strP.at(i-1),strS.at(j-1))) <  -1E-12)
	    {
	      compval = (dp_table_prev_col[(j-1)] + c_sm->GetScore(strP.at(i-1),strS.at(j-1)));	//compval = diagonal + match score
	      //}
	      //cout<<"\t\tcompval after match/mismatch:"<<compval<<";";
	      c_traceback_table->SetTableEntry(i,j, UPLEFT,0,0);//c_traceback_table[i+j*(lenP+1)].SetLinks(UPLEFT);//for this one, we don't have to set the #numOfIndels, since there is none
	    }
	  
	  //here, we try to update information after compare values 	
	  if((compval-0) < 1E-15) //found a zero or smaller value
	    {
	      compval = 0;
	      //The following is not necessary, since everything so far is default to ZERO
	      //c_traceback_table[i+j*(lenP+1)]=ZERO;
	      c_traceback_table->SetTableEntry(i,j, ZERO,0,0);
	      //we don't set numOfIndels;keep default
	    }
	  else //found a larger than zero value, this is where we need to check whether this is good starting for a new path
	    {
	      //check its
	      if(c_traceback_table->GetLink(i,j)==UPLEFT)//[i+j*(lenP+1)].GetLinks()==UPLEFT)//it is possible that this is a new path start. 
		//there is no way that a path starting from an insertion, it has to be a match
		{
		  if(c_traceback_table->GetLink(i-1,j-1)==ZERO)//[(i-1)+(j-1)*(lenP+1)].GetLinks()==ZERO)
		    {
		      //cout<<"\t\t&&&&&&&&&&&&ADDING A NEW PATH"<<endl;
		      //this is new path starting. since it following a zero and then a match.
		      //we need to add a entry to a path vec to keep track this new path
		      unsigned int tempIndexArr[2];tempIndexArr[0]=i;tempIndexArr[1]=j;
		      c_path_vec.push_back(new Path(tempIndexArr, tempIndexArr,compval ));
		      pathElementTable[i+j*(lenP+1)]=c_path_vec.size()-1;
		    }
		  else //not a zero, it was a larger than zero one, it has to be on a path
		    {
		      //first get the path from the vector
		     
		      Path* tempPath=c_path_vec.at(pathElementTable[i-1+(j-1)*(lenP+1)]);
		      //now compare to update or not
		      if(compval>tempPath->GetOptimalValue())
			{
			  //cout<<"\t\t\t***********ffound a maximum one for the existing path, update info::path index"<<pathElementTable[i-1+(j-1)*(lenP+1)]<<endl;
			  unsigned int tempIndexArr[2];tempIndexArr[0]=i;tempIndexArr[1]=j;
			  tempPath->SetOptimalIndex(tempIndexArr);
			  tempPath->SetOptimalValue(compval);
			  //tempPath->
			}
		      pathElementTable[i+j*(lenP+1)]=pathElementTable[i-1+(j-1)*(lenP+1)];
		    }
		}
	      else  //for not a UPLEFT entry, means this one is a indel, the score is descreasing, well so 1)it has to follow 
		//a path 2)it is no way to the maximum one so far, we only need to update the pathElement table information, don't have
		//update the path entry
		{
		  LinkBack tempLink=c_traceback_table->GetLink(i,j);//[i+j*(lenP+1)].GetLinks();//it can not be ZERO, can not be UPPERLEFT
		  
		  if(tempLink==UP)
		    {
		      //pathElementTable[i+j*(lenP+1)]=pathElementTable[i+(j-c_traceback_table[i+j*(lenP+1)].GetNumOfIndels())*(lenP+1)];
		      pathElementTable[i+j*(lenP+1)]=pathElementTable[i+(j-c_traceback_table->GetNumOfIndels(i,j))*(lenP+1)];
		    }
		  else
		    {
		      if(tempLink==LEFT)
			{
			  //pathElementTable[i+j*(lenP+1)]=pathElementTable[(i-c_traceback_table[i+j*(lenP+1)].GetNumOfIndels())+j*(lenP+1)];
			  pathElementTable[i+j*(lenP+1)]=pathElementTable[(i-c_traceback_table->GetNumOfIndels_s(i,j))+j*(lenP+1)];
			}
		    }
		}  
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
	  //cout<<"\t\trunning score:"<<compval<<endl;

	  //find a best one so far
	  if(c_score<compval)
	    {
	      c_score=compval;
	      c_optimalIndex[0]=i;
	      c_optimalIndex[1]=j;
	    }
	}//end of col
      
      //reset the prev one to curr col and make it ready for next loop
      double* tempP=dp_table_curr_col;
      dp_table_curr_col=dp_table_prev_col;
      dp_table_prev_col=tempP;
    }//end of row

  //clean up
  delete[] dp_table_prev_col;
  delete[] dp_table_curr_col;
  delete[] maximumGapValue_subject;
  delete[] maximumGapIndex_subject;
  delete[] pathElementTable;
  
  //cout<<"+++++++++done with the alignment "<<endl;
  //	flush(cout);
  //delete gm;
  //traceback_table will be deleted upon destruction
}//end of the localAlign

//***********For Path class
/*
Path::Path():c_optimalValue(0)
{
  c_optimalIndex[0]=0;
  c_optimalIndex[1]=0;
  c_startIndex[0]=0;
  c_startIndex[1]=0;
  
}
Path::Path( const unsigned int _OptimalIndex[2], const unsigned int _startIndex [2], const double& _optimalValue ):
  c_optimalValue(_optimalValue)
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
const unsigned int* Path::GetOptimalIndex()
{
  return c_optimalIndex;
}
const unsigned  int* Path::GetStartIndex()
{
  return c_startIndex;
}
double Path::GetOptimalValue()
{
  return c_optimalValue;
}
*/