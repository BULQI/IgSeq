#include "MarkovChainGapModel_454.hpp"
#include <iostream>
using namespace std;

//**************************************
//issue: currently when a regular gap following a MarkovChain Gap, there is no gap opening penalty
//   ???????????do we need to set a new issue to account for this penalty!!!!!!!
//**********************************


//this model is for 454 sequencing reading, to model the umbiguity of stretch of identical nucleotides run. it is using affine model too

MarkovChainGapModel_454::MarkovChainGapModel_454(const double& _gopen, const double& _gextension
						 ,const string& _patternStr, const string& _subjectStr
						 )
:AffineGapModel(_gopen, _gextension), c_patternString(_patternStr), c_subjectString(_subjectStr)
{
  //empty
}

MarkovChainGapModel_454::~MarkovChainGapModel_454()
{
  //empty
}
double MarkovChainGapModel_454::GapValue(const TracebackTable* _tbTable, const unsigned int& _patternIndex, const unsigned int& _subjectIndex, const bool& _patternGap,
					 const double& _prevEntryValue,
					 double& _MaxGapValue, unsigned int& _MaxGapIndex) const
{
  //to start, there are only two possible cases: 1) the extended from the MaxGapIndex;2)a new open gap
  //for both cases, we need to check whether this is a stretch of identical nucleotides and then decide how we 
  //are going to apply cost.
  //first take care of extended from the MaxGapIndex
  //check to see whether we need to use the Markov Chain Model
  //cout<<"\t######calling markovChainModel:"<<endl;
  char curr_char;
  bool regular_cost_flag=false;
  double extendedGapValue;
  
  //do gaps for extendedGapValue, in here, we assume it is extended from the maxigap value one, all gaps from 
  //current one till the maxiGapIndex one.
  if(_patternGap)//doing row gap now
    {
      //cout<<"\t######gap on pattern"<<endl;
      curr_char=c_subjectString[_subjectIndex-1];
      //?????????will have to use the linkback table to check (???????????????????)
      //no, since we are considering gaps on pattern string, we only check the subjectstring char till the GapMaxIndex, and keep the pattern char same
      bool zeroOnes=true;
      for(unsigned int i=1;_subjectIndex-i>=_MaxGapIndex&&_subjectIndex-i!=0;i++)
	{
	  zeroOnes=false;
	  //cout<<"c_patternString["<<_patternIndex-1<<"]:"<<c_patternString[_patternIndex-1]<<"c_subjectString["
	  //	  <<_subjectIndex-i-1<<"]:"<<c_subjectString[_subjectIndex-i-1]<<";curr_char:"<<curr_char<<endl;

	    if(c_subjectString[_subjectIndex-i-1]!=curr_char)
	    {
	      
	      regular_cost_flag=true;
	      break;
	    }
	    
	    //for the MaxGapIndex ones, we need to make sure the chars are identical, they are match
	    if(_subjectIndex-i==_MaxGapIndex)
		{
		  if(c_subjectString[_subjectIndex-i-1]!=c_patternString[_patternIndex-1])
		    {
		      regular_cost_flag=true;
		      break;
		    }
		}
	  
	}
      if(regular_cost_flag||zeroOnes) //regular
	{
	  //cout<<"\t######regular affine cost"<<endl;
	  extendedGapValue=_MaxGapValue+c_gextension;
	}
      else //Markov chain
	{
	  //cout<<"\t######Markov chain affine cost"<<endl;
	  //need to check how long the run is
	  unsigned int runLen=0;
	  //if we are here, _patternIndex-_MaxGapIndex won't be zero, and _subjectIndex-(_patternIndex-_MaxGapIndex)won't either
	  /*while(c_subjectString[_subjectIndex-(_patternIndex-_MaxGapIndex)-runLen-1]!=curr_char||c_patternString[_MaxGapIndex-runLen-1]!=curr_char)
	    {
	      
	      runLen++;
	      if(_subjectIndex-(_patternIndex-_MaxGapIndex)-runLen-1<0&&_MaxGapIndex-runLen-1<0)
		{
		  break;
		}
		}*/
	  //try to trace back using trace back table in order to figure out the length of the run of the identical nucleotide
	  while(_tbTable->GetLink(_patternIndex-runLen,_MaxGapIndex-runLen )==UPLEFT)
	    {//here we don't have to check whether the index will get negative, since the zero ones won't pass the check of while clause.
	      if(c_subjectString[_MaxGapIndex-runLen-1]==curr_char&&c_patternString[_patternIndex-runLen-1]==curr_char)
		//(c_subjectString[_subjectIndex-(_patternIndex-_MaxGapIndex)-runLen-1]==curr_char||c_patternString[_MaxGapIndex-runLen-1]==curr_char)
		runLen++;
	      else
		break;
	    }
	  //get the gap len
	  unsigned int gapLen=_subjectIndex-_MaxGapIndex;

	  //get the cost
	  extendedGapValue=_MaxGapValue+GetMarkovChainGapValue(c_gopen, c_gextension, runLen, gapLen, true);	  
	}
    }
  else //doing subject gap, gap on subject strings
    {
      //cout<<"\t#####gap on subject:"<<endl;
      curr_char=c_patternString[_patternIndex-1];
      bool zeroOnes=true;
      for(unsigned int i=1;_patternIndex-i>=_MaxGapIndex&&_patternIndex-i!=0;i++)
	{
	  //cout<<"c_patternString["<<_patternIndex-i-1<<"]:"<<c_patternString[_patternIndex-i-1]<<"c_subjectString["
	  //	  <<_subjectIndex-1<<"]:"<<c_subjectString[_subjectIndex-1]<<";curr_char:"<<curr_char<<endl;
	  zeroOnes=false;
	  if(c_patternString[_patternIndex-i-1]!=curr_char)
	    {
	      regular_cost_flag=true;
	      break;
	    }
	  
	}
      if(regular_cost_flag||zeroOnes) //regular
	{
	  //cout<<"\t######regular affine cost"<<endl;
	  extendedGapValue=_MaxGapValue+c_gextension;
	}
      else //Markov chain
	{
	  //cout<<"\t######Markov Chain affine cost"<<endl;
	  //need to check how long the run is
	  unsigned int runLen=0;
	  
	  while(_tbTable->GetLink(_MaxGapIndex-runLen, _subjectIndex-runLen )==UPLEFT)
	    {//here we don't have to check whether the index will get negative, since the zero ones won't pass the check of while clause.
	      if(c_patternString[_MaxGapIndex-runLen-1]==curr_char&&c_subjectString[_subjectIndex-runLen-1]==curr_char)
		runLen++;
	      else
		break;
	    }
	  //get the gap len
	  unsigned int gapLen=_patternIndex-_MaxGapIndex;

	  //get the cost
	  extendedGapValue=_MaxGapValue+GetMarkovChainGapValue(c_gopen, c_gextension, runLen, gapLen, true);	  
	}
    }//patterngap or subjectgap

  //now we want to solve the problem of figure out the open new gap value
  unsigned int runLen=0;
  if(_patternGap)//pattern gap
    {
      while(_tbTable->GetLink(_patternIndex-runLen,_subjectIndex-runLen-1 )==UPLEFT)
	    {//here we don't have to check whether the index will get negative, since the zero ones won't pass the check of while clause.
	      //cout<<"c_patternString["<<_patternIndex-runLen-1<<"]:"<<c_patternString[_patternIndex-runLen-1]<<"c_subjectString["
	      //	  <<_subjectIndex-runLen-1-1<<"]:"<<c_subjectString[_subjectIndex-runLen-1-1]<<";curr_char:"<<curr_char<<endl;

	      //it is important to mention that string character index is one less the index in the dynamic programming table
	      if(c_patternString[_patternIndex-runLen-1]==curr_char&&c_subjectString[_subjectIndex-runLen-1-1]==curr_char)
		runLen++;
	      else
		break;
	    }
    }
  else
    {
      while(_tbTable->GetLink(_patternIndex-runLen-1,_subjectIndex-runLen )==UPLEFT)
	    {//here we don't have to check whether the index will get negative, since the zero ones won't pass the check of while clause.
	      if(c_patternString[_patternIndex-runLen-1-1]==curr_char&&c_subjectString[_subjectIndex-runLen-1]==curr_char)
		runLen++;
	      else
		break;
	    }
    }
  //cout<<"\t^^^^^^^^^^runlen:"<<runLen<<endl;
  double newOpenGapValue;
  if(runLen==0) //regular one
    {
      newOpenGapValue=_prevEntryValue+c_gopen+c_gextension;
    }
  else
    {
      newOpenGapValue =_prevEntryValue+GetMarkovChainGapValue(c_gopen, c_gextension, runLen, 1, false);
    }
  
  if(extendedGapValue > newOpenGapValue)
    {
      _MaxGapValue=extendedGapValue;
      //max gap index keeps the same
      //return extendedGapValue;
    }
  else  //new gap, so we only need to point to one previous one
    {
      if(_patternGap) //gap on patter string
	{
	  _MaxGapIndex=_subjectIndex-1;
	}
      else
	{
	  _MaxGapIndex=_patternIndex-1;
	}
      _MaxGapValue=newOpenGapValue;
      //return newOpenGapValue;
    }
  //cout<<"\t^^^^^^^^^^MaxGapValue:"<<_MaxGapValue<<";MaxGapIndex:"<<_MaxGapIndex<<endl;
  return _MaxGapValue;
}

double MarkovChainGapModel_454::GetMarkovChainGapValue(const double& _regularGapOpen, const double& _regularGapExtension, const unsigned int& _runLen, const unsigned int& _gapLen, const bool& _extensionFlag) const
{
  int etf=1;
  if(_extensionFlag)
    etf=0;
  return (etf*_regularGapOpen+_regularGapExtension)*(1/(1+exp(-20.0*(((double)_gapLen)/(_gapLen+_runLen)-0.5))));
} 
