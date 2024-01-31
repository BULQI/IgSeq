#include "OverlapAlignment.hpp"
#include <stdlib.h>
#include <iostream>
#include <vector>

#include "GapModel.hpp"
#include "AffineGapModel.hpp"
#include "MarkovChainGapModel_454.hpp"

using namespace std;

OverlapAlignment::OverlapAlignment(SequenceString* _pattern, SequenceString* _subject, 
			       ScoreMatrix* _m, const double& _gopen, 
			       const double& _gextension, const double& _scale, const short& _typeOfGapModel):
  PairwiseAlignment(_pattern, _subject, _m, _gopen, _gextension, _scale, _typeOfGapModel)
{
  //now we need to
  //c_alignmentArr=new AlignmentString[this->c_numOfAlignments];
  //c_scoreArr=new double[this->c_numOfAlignments];
  //dp_table=NULL;
  //traceback_table=NULL;
  this->align();
  this->traceBack();
}
  
OverlapAlignment::~OverlapAlignment()
{
  //the base class destructor is called automatically

  //here we only need to take care other issues
  
}

void OverlapAlignment::traceBack()
{
  //vector<unsigned int[2]> OptimalValueIndex; //the index that current best score resides
  //vector<unsigned int[2]> startIndex;//where the 
  //vector<double[2]> curr_and_max;

  //cout<<"doing trace back in the local alignment"<<endl;
  
  //starting from optimalIndex going backing
  unsigned int i=c_optimalIndex[0];
  unsigned int j=c_optimalIndex[1];

  //cout<<"starting from ("<<i<<","<<j<<")"<<endl;
  string c_pattern_wg;//aligned string with gap
  string c_subject_wg;//aligend string with gap

  
  //the while criteria, could be this or i==0&&j==0
  //while(c_traceback_table[i+j*(c_pattern->GetLength()+1)].GetLinks()!=ZERO)  //keep going till we reach a zero
  while (c_traceback_table->GetLink(i, j)!=ZERO)
    {
      //cout<<"\t("<<i<<","<<j<<"):";
      unsigned int currentIndex;
      //check the link table, decide where to go
      switch(c_traceback_table->GetLink(i, j))//(c_traceback_table[i+j*(c_pattern->GetLength()+1)].GetLinks())
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
	  for(unsigned int k=0;k<c_traceback_table->GetNumOfIndels(currentIndex, j);k++)//c_traceback_table[currentIndex+j*(c_pattern->GetLength()+1)].GetNumOfIndels();k++)
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
  
}


//the variant of alignment using less memeory
//this one we don't keep the original mxn dp table,
//but instead we keep only one column (in fact for coding
//purpose, we keep two columns to make the algorithm working
//efficient;
void OverlapAlignment::align()
{
  	//Create empty dynamic programming table
  unsigned lenP=c_pattern->GetLength();
  unsigned lenS=c_subject->GetLength();
  double* dp_table_prev_col=new double[(lenS+1)];//one extra on this, to deal with the beging of the column
  double* dp_table_curr_col=new double[(lenS+1)];//one extra on this, to deal with the beging of the column

  //GapModel* gm=new AffineGapModel(c_gapOpen, c_gapExtension);
  //GapModel* gm=new MarkovChainGapModel_454(c_gapOpen, c_gapExtension, c_pattern->GetSequence(), c_subject->GetSequence());

  //here we keep two columns, to effieciently and easily manipulate the column. it cold be only one
  c_traceback_table=new TracebackTable(c_pattern, c_subject, 1);//[(lenP+1)*(lenS+1)];
  //the following is to set all the first rows and cols to zero. 
  /*c_traceback_table[0].SetLinks(ZERO);
  for(unsigned int i=1;i<=lenP;i++)
    {
      c_traceback_table[i+0*(lenP+1)].SetLinks(ZERO);
      c_traceback_table[i+0*(lenP+1)].SetNumOfIndels(0);
    }
  for(unsigned int j=1;j<=lenS;j++)
    {
      c_traceback_table[0+j*(lenP+1)].SetLinks(ZERO);
      c_traceback_table[0+j*(lenP+1)].SetNumOfIndels(0);
    }
  */

  c_score=-1E9;
  
  //intialize the dp table, the first row=gapsExtended and first column=gapsExtended
  //the row is Pattern, column is Subject. 
  dp_table_curr_col[0]=0;//this is the first row

  dp_table_prev_col[0]=0;
  for(unsigned int i=1;i<=lenS;++i)
    {
      dp_table_prev_col[i]=0;//this is the first column
    }

  //cout<<"successfully created the empty tables and now go to get the score!\n";

  //score table with S-W
  double compval = 0;
  string strP=c_pattern->GetSequence();
  string strS=c_subject->GetSequence();
  
  double* maximumGapValue_subject=new double[lenS+1]; //we still keep this same len as the curr/prev col, but we will never use the first one.
  unsigned int* maximumGapIndex_subject=new unsigned int[lenS+1];//just as above, this one is used to keep record of the maximum Gap Value so far, and the first one [0] is not used.
  
  //maximumGapValue[0]=-1E9;
  //maximumGapValue[1]=c_gapOpen+c_gapExtension*1;
  //we will start doing the job at column 1, so set intial value of gapIndex=0 to infinity
  for(unsigned int i=0;i<=lenS;i++)
    {
      maximumGapValue_subject[i]=-1E60;
      maximumGapIndex_subject[i]= 0;
      //  maximumGapIndex[i]=0;
      //maximumGapIndex[1]=1;//to the
    } 

  //for pattern gap, this maximumGapValue, a scalar is OK. since we keep this as a running one, like the current column
  double maximumGapValue_pattern=-1E60;
  unsigned int maximumGapIndex_pattern=0;


  double maximumValueLastRowCol=-1E70;
  unsigned int maximumIndexLastRowCol[]={0,0};

  //we have to be very careful here
  //the dimension of the dp table and the strings are not same.
  //dp table is one row/col more than the strings, since we need to 
  //have the zero/starting row or column.so that means, when we dealing with
  //dp table, it is from 1->lenP or 1->lenS, but the string, starts from zero,
  //so we are doing strings i-1
  for(unsigned int i = 1; i <= lenP; ++i) //for all values of strA
    {	//now we are starting a new column, we need to
      //1)first exchange the curr and prev col(actually, it is better to do this one at end of each loop, otherwise it will be trouble in the beginning of the first loop
      //2)reset the curr_col first one to gap, it should be zero anyway.
      
      
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
	      }*/
	  compval=maximumGapValue_subject[j];
	  //if(compval < maximumGapValue[j]) 
	  //  {	    //if cell above has a greater value 
	      
	  //  compval = maximumGapValue[j];		//set compval to that square
	  c_traceback_table->SetTableEntry(i, j, LEFT, i-maximumGapIndex_subject[j]);//[i+j*(lenP+1)].SetLinks(LEFT);
	  //c_traceback_table[i+j*(lenP+1)].SetNumOfIndels(i-maximumGapIndex[j]);
	      // }
	      //cout<<"maximumGapValue[j]:"<<maximumGapValue[j]<<",";
	      //cout<<"campval after rowGap:"<<compval<<";";

	  //this is the pattern/column gap, we keep it same as we are doing with the whole dp table
	  c_gm->GapValue(c_traceback_table, i,j,true, dp_table_curr_col[j-1], maximumGapValue_pattern, maximumGapIndex_pattern);
	  /*for(int k=j-1; k>0; --k) 
	    {		//check all sub columns
				
	      if(compval < ((dp_table_curr_col[k]) + (c_gapOpen + (c_gapExtension *(j- k))))) 
		{	
		  //if square to the left has the highest valu
					
		  compval = ((dp_table_curr_col[k]) + (c_gapOpen + (c_gapExtension *(j- k))));    //set compval to that square
		  //c_traceback_table[i+j*(lenP+1)].SetLinks(UP);
		  //c_traceback_table[i+j*(lenP+1)].SetNumOfIndels(j-k);
		  c_traceback_table->SetTableEntry(i,j, UP, j-k);
		}
		}	*/
	  if(compval<maximumGapValue_pattern)
	    {
	      //set compval to that square
	      compval=maximumGapValue_pattern;
	      c_traceback_table->SetTableEntry(i,j,UP, j-maximumGapIndex_pattern);
	    }	
	  //cout<<"compval afer col gap:"<<compval<<endl;

	  //*********do Macth/Mismatch now
	  if(compval < (dp_table_prev_col[(j-1)] + c_sm->GetScore(strP.at(i-1),strS.at(j-1)))) 
	  {	    //if cell upleft leads to a greater value 
	      
	  //  compval = maximumGapValue[j];		//set compval to that square
	    compval = (dp_table_prev_col[(j-1)] + c_sm->GetScore(strP.at(i-1),strS.at(j-1)));
	    //c_traceback_table[i+j*(lenP+1)].SetLinks(UPLEFT);
	    //c_traceback_table[i+j*(lenP+1)].SetNumOfIndels(0);
	    c_traceback_table->SetTableEntry(i,j, UPLEFT, 0);
	  }
	  //compval = diagonal + match score
	  //}
	  //cout<<"\t\tcompval after match/mismatch:"<<compval<<";";
	  //c_traceback_table[i+j*(lenP+1)].SetLinks(UPLEFT);//for this one, we don't have to set the #numOfIndels, since there is none

	  
	  /*if((compval-0) < 1E-10) 
	    {
	      compval = 0;
	      //The following is not necessary, since everything so far is default to ZERO
	      c_traceback_table[i+j*lenP]=ZERO;
	      //we don't set numOfIndels;keep default
	    }
	  */
	  /*witch (c_traceback_table[i+j*(lenP+1)].GetLinks())
	    {
	    case UP:
	      cout<<"\t\tlink is UP;#indels is"<<c_traceback_table[i+j*(lenP+1)].GetNumOfIndels()<<endl;
	      break;
	    case LEFT:
	      cout<<"\t\tlink is LEFT;#indels is"<<c_traceback_table[i+j*(lenP+1)].GetNumOfIndels()<<endl;
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

	  //c_score=compval;
	  //find a best one so far
	  /*if(c_score<compval)
	    {
	      
	      
	      }*/

	}//end of col
      //at the end of row, we need to keep track of biggest one on the last row.
      if(maximumValueLastRowCol<dp_table_curr_col[lenS])
	{
	  //cout<<"&&&&&&found a maximum score ("<<i<<","<<lenS<<"0="<<dp_table_curr_col[lenS]<<endl;
	  maximumValueLastRowCol=dp_table_curr_col[lenS];
	  maximumIndexLastRowCol[0]=i;
	  maximumIndexLastRowCol[1]=lenS;
	}
      //reset the prev one to curr col and make it ready for next loop
      double* tempP=dp_table_curr_col;
      dp_table_curr_col=dp_table_prev_col;
      dp_table_prev_col=tempP;
    }//end of row

  //now need to go through the col to get the biggest one
  for(unsigned int k=1;k<=lenS;k++)
    {
      if(maximumValueLastRowCol<dp_table_prev_col[k]) //using prev col becase we have switched the prev and curr col
	{
	  //cout<<"&&&&&&found a maximum score ("<<lenP<<","<<k<<"0="<<dp_table_prev_col[k]<<endl;
	  maximumValueLastRowCol=dp_table_prev_col[k];
	  maximumIndexLastRowCol[0]=lenP;
	  maximumIndexLastRowCol[1]=k;
	}
    }
  //for the global alignment,we always traceback from the lowerright corner
  c_optimalIndex[0]=maximumIndexLastRowCol[0];
  c_optimalIndex[1]=maximumIndexLastRowCol[1];
  
  c_score=maximumValueLastRowCol;

  //clean up
  delete[] dp_table_prev_col;
  delete[] dp_table_curr_col;
   delete [] maximumGapValue_subject;
  delete [] maximumGapIndex_subject;

  //delete gm;
  //traceback_table will be deleted upon destruction
}//end of the localAlign
