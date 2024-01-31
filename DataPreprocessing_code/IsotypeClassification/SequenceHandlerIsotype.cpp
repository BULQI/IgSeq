#include "SequenceHandler.hpp"
#include "Accessory/FastaHandler.hpp"
#include <fstream>
#include <iostream>
#include <stdlib.h>
#include <sstream>

#include "SequenceHandlerIsotype.hpp"
#include "OverlapAlignment.hpp"
#include "AlignmentString.hpp"
#include "LocalAlignment.hpp"
#include "Accessory/FastqHandler.hpp"

//#define debug 


//**********************************************

//the _adaptorname is the name for the aligned sequence we will be 
// _adaptorName is of format "adaptor A:Barcode:IgM/G/D"
static unsigned int LookUpVectorIndex(const string& _adaptorName, vector<SequenceString>& _vecPrimer)
{
  //cout<<"&&&&&&looking up"<< _adaptorName <<endl;
  //use the name as the key
  for(unsigned int i=0;i<_vecPrimer.size();i++)
    {
      //cout<<"\t%%%%%%"<<_vecPrimer.at(i).GetName()<<endl;
      unsigned found;
      found=_adaptorName.compare(_vecPrimer.at(i).GetName()); //will return -1( string::npos) if can not find it.
      //cout<<"\t%%%%%%%%%%%"<<(signed)found<<":npos: "<<string::npos<<endl;
      //if(found!=std::string::npos)
      if(found==0)
	{
	  return i;
	}
    }
  cout<<"***ERROR: can not find which catogory to store data (IgM/G/D)"<<endl;
  return -1;
}

//this is a helper function to do align with isotype constant sequence and return whether this is a valid alignement
// valid means, this current alignment is better for score, long enough, within offset, pass matchRate threshold 
//	return true if this is a good valid alignment and false otherwise
bool alignConstant( const SequenceString* seq, const SequenceString* iso, 
		const ScoreMatrix* _sm, const double& _gopen, 
		 const double& _gextension, const double& _scale,  const short& _gapModel, const mapType& _mtype,
		 const double& _matchRateThreshold, const unsigned& _minimumOverlapLength, const unsigned& _offset, 
		 //output now
		 AlignmentString* bestAlign,   double* bestScore /*either 5' or 3' best score*/
	)
	{
		//cout<<"Isotype set:"<<j<<endl;
	      //cout<<_vecIsotype.at(j).toString()<<endl;
		  //flush(cout);
	      LocalAlignment lal (seq, iso, _sm, _gopen, _gextension
						, 1 /*scale*/, 1  /*numOfLocalAlignments*/, _gapModel);
		  //cout<<"\tdone with local alignment...."<<endl;
		  //flush(cout);
		  if(lal.GetNumberOfAlignments()==0)
		  {
			  return false;//don't do anyting.
		  }
	      AlignmentString as=lal.GetAlignment();
#ifdef debug  
			cout<<"*****Seq name:"<<seq->GetName()<<endl;
			cout<<"\talignment score:"<<as.GetScore()<<endl;
			flush(cout);
#endif	      
		  //need to get the match rate
	      string strPattern=as.GetPattern(true);
	      string strSubject=as.GetSubject(true);
#ifdef debug
	      cout<<"\tstrPattern:"<<strPattern<<endl;
	      cout<<"\tstrSubject:"<<strSubject<<endl;
#endif
	      double match_rate=1-CompareStrings(strPattern, strSubject)/((double)strPattern.length());
	      //cout<<"\tcompare:"<<CompareStrings(strPattern, strSubject)<<";length():"<<strPattern.length()<<endl;
#ifdef debug	      
			cout<<"\tmatch_rate:"<<match_rate<<" match rate thresld:"<<_matchRateThreshold<<endl;
			cout<<"\tstrPattern.length:"<<strPattern.length()<<"; minimumOverlpLength:"<<_minimumOverlapLength<<endl;
			cout<<"\toffset threshold:"<<_offset<<endl;
			flush(cout);
#endif		
	      //check the best score and match rate
	      if(as.GetScore()>*bestScore&&match_rate>_matchRateThreshold&&strPattern.length()>_minimumOverlapLength)
		{
#ifdef debug
			cout<<"length:"<<strPattern.length()<<" and length thresld:"<<_minimumOverlapLength<<endl;
#endif
		  //we need to see the offset too, only important to the pattern, here the 
		  //the pattern is the long one, we align the primer sequence against
		  //the primer should be in the beginning, not too far --- 5' mapping,
		  //or the primer is in the end, not too far from the end --- 3' mapping
		  unsigned temp_offset=as.GetPatternIndexStart();
		  
		  if(_mtype==ThreePrime)
			  temp_offset=seq->GetLength()-as.GetPatternIndexEnd();
		  //cout<<"the alignment offset :"<<temp_offset<<endl;
		  if(temp_offset<_offset )
		    {//we good
#ifdef debug
		      cout<<"\t***get one bigger"<<endl;
			  cout<<"starting offset"<<as.GetPatternIndexStart()<<"and offset:"<<_offset<<endl;
#endif
		      *bestScore=as.GetScore();
		      *bestAlign=as;//with name
		      
		      return true;
		    }
		}
			return false;
	}//end of function alignconstant


void MappingIsotypes(const vector<SequenceString>& _vecSeq, /*this is the sequence data that we want to find isotypes*/
		     vector<SequenceString>& _vecIsotype, /*this is the isotype sequences*/
		     const mapType& type, /*indicating whether it is 5'prime or 3' prime*/
		     ScoreMatrix* _sm, const double& _gapOpen, const double& _gapExtension,
		     const double& _matchRateThreshold, const unsigned _minimumOverlapLength,
		     const unsigned int& _offset, //const unsigned int& _offsetReverse, 
		     const string& _R1_fname,
		     const string& _R2_fname,
		     const bool& _demux,
			 const vector<string>& _vecSeq_Q, /*read 1 sequence quality*/
			 const vector<SequenceString>& _vecSeq_R2, /*R2 sequence data that we want to find isotypes*/
			 const vector<string>& _vecSeq_Q_R2/*r2 quality.*/
		     )
{
		
//	cout<<"%%%%%%starting "<<endl;
//			cout<<"\tvecSeq_Q.size():"<<_vecSeq_Q.size()<<endl;
//			cout<<"\tvecSeq_R2.size():"<<_vecSeq_R2.size()<<endl;
//			cout<<"\tvecSeq_Q_R2.size():"<<_vecSeq_Q_R2.size()<<endl;
			
  short gapModel=0; //0, affine gaps;1, 454 gap model.
  unsigned int numOfSeqsUnit=20000;
  unsigned int timeOfWriting=1;
  int* numOfWritesDone_mapped =new int [_vecIsotype.size()+1];
  for(unsigned int i=0;i<_vecIsotype.size()+1;i++)
    {
      numOfWritesDone_mapped[i]=0;
    }

  int numOfWritesDone_unmapped=0;
  ios_base::openmode mode=ofstream::trunc; //by default, we need to clear out the file.

  //now we need to prepare the output vector, now we will prepare the outfile for each entry in _vecForward
  //they are the alleles for different isotypes (IgG/M), different subtypes (IgG1/2/3/4) and alleles IgG1*01/*02/*03, etc
  //
  vector<SequenceString>* pt_vec_mapped=new vector<SequenceString>[_vecIsotype.size()+1]; //one more for holding unknown isotype
  vector<SequenceString>* pt_vec_demux=new vector<SequenceString>[_vecIsotype.size()+1];//holding only the sequences by isotype, //but no mapping as output
  vector<string>* pt_vec_demux_Q=NULL;
  if(_vecSeq_Q.size()>0)
  {
		pt_vec_demux_Q=new vector<string>[_vecIsotype.size()+1];
  }
  vector<SequenceString>* pt_vec_demux_R2=NULL;
  if(_vecSeq_R2.size()>0)
  {
	  pt_vec_demux_R2=new vector<SequenceString>[_vecIsotype.size()+1];
  }
  vector<string>* pt_vec_demux_Q_R2=NULL;
  if(_vecSeq_Q_R2.size()>0)
  {
		pt_vec_demux_Q_R2=new vector<string>[_vecIsotype.size()+1];
  }
  
  //vector<SequenceString>* pt_vec_mapForward=new vector<SequenceString>[_vecForward.size()+1];//one more for holding unkown isotype
  //vector<SequenceString> vec_mapReverse;
  vector<SequenceString> vec_mapNone;
  vector<string> vec_mapNone_Q;
  vector<SequenceString> vec_mapNone_R2;
  vector<string> vec_mapNone_Q_R2;
  
  //the stats are arranged by isotypes(matched)
  //holding the position where within the subject (isotype) the break is.
  vector<unsigned int>* pt_vec_map_pos_end=new vector<unsigned int>[_vecIsotype.size()];
  vector<unsigned int>* pt_vec_map_pos_start=new vector<unsigned int>[_vecIsotype.size()];

  //holding the position where within the pattern (seq) the break is.
  vector<unsigned int>* pt_vec_map_pos_end_seq=new vector<unsigned int>[_vecIsotype.size()];
  vector<unsigned int>* pt_vec_map_pos_start_seq=new vector<unsigned int>[_vecIsotype.size()];

  //holding the stats of alignment overlap and match rate.
  vector<double>* pt_vec_map_match_rate=new vector<double>[_vecIsotype.size()];
  vector<unsigned int>* pt_vec_map_overlap=new vector<unsigned int>[_vecIsotype.size()];

  vector<unsigned int>* pt_vec_map_seq_len=new vector<unsigned int>[_vecIsotype.size()];

  //vector<unsigned int[2]> vec_mapCrossOver_pos_end;
  //vector<unsigned int>vec_mapBreakOutsideConstant_pos_end;

  //now, we need to prepare the output files
  AlignmentString tempAS;
  //AlignmentString* tempAS_arr=NULL;
  //unsigned int numOfLocalAlignments=1; //we most likely will do overlap alignment. so this might not be necessary
  string strPattern;
  string strSubject;
  
  double match_rate=0.0;
  cout<<"\n";
  for(unsigned int i=0;i<_vecSeq.size();i++)
    {
      //printing progress
#ifndef debug     
				if(i/numOfSeqsUnit *numOfSeqsUnit==i)
					{
#endif
	  cout<<"..."<<i<<"/"<<_vecSeq.size();

#ifndef	debug  			
						flush(cout);
	  
	  //need to write to file.
					}
#endif
      
      //here we separate the two different type of mapping, in case we need to do outputing differently
      double best5PrimeScore= -10000000;
      AlignmentString	best5PrimeAlign;
      double best3PrimeScore= -1000000;
      AlignmentString best3PrimeAlign;

      //The following two are used to remember for the best alignment for the Isotype seq index
      unsigned int 	best5PrimeIndex=0; //what is this for? to remeber index of local alignment? do we need it for overlap alignment
      unsigned int	best3PrimeIndex=0;
	
	// the below two was used to indiciate whether a valid aignment was found??
	//		if not, mapnone was recorded.
      bool 	found5PrimeFlag=false;
      bool	found3PrimeFlag=false;
	  //cout<<"sequence"<< i<<":"<<endl;
	  //    cout<<_vecSeq.at(i).toString()<<endl;
	//	  flush(cout);
      //cout<<"Map isotype set"<<endl;
	  //flush(cout);
      for(unsigned int j=0;j<_vecIsotype.size();j++)
	{
		//note, for each alignment, now when we ask for the tempAS_arr, we ask for a deep copy, meanning we get some 
		//new memory allocated, so we need to take care of this by deleting/cleaning up.
		//so that we don't allow memory leak. 
		//if(tempAS_arr!=NULL) //<-----------
		//	delete [] tempAS_arr;
#ifdef debug
		cout<<"Isotype set:"<<j<<endl;
	    cout<<"\t"<<_vecIsotype.at(j).toString()<<endl;		  
		//flush(cout);
#endif 
	  //now we need to separate out the two different mapping case, 5' and 3'
	  if(type==FivePrime)
	    {
		
		  //check to see 
	      bool currentFlag=alignConstant(&(_vecSeq.at(i)), &(_vecIsotype.at(j)), _sm, _gapOpen, _gapExtension
					,1, gapModel, FivePrime, 
					_matchRateThreshold, _minimumOverlapLength, _offset, 
		 //output now
		  &best5PrimeAlign,   &best5PrimeScore );
		  if(currentFlag)
		  {
			  found5PrimeFlag=currentFlag;
			  best5PrimeIndex=j;
		  }
	    }
	  else //in this case this is 3' mapping
	    {
#ifdef debug	      
			cout<<"map reverse set"<<endl;
#endif
	      //3' should be mapped on the end of the reads, need to reverse complement the sequence too
	      //for(unsigned int k=0;k< _vecReverse.size();k++)
	      //{
		  //cout<<"start doing k:"<<k<<endl;
		  SequenceString reverseComplementReverse=ReverseComplement(_vecIsotype.at(j));
#ifdef debug		  
					  cout<<"\t**after revcomp"<<endl;
					  
					  cout<<"\t\t "<<reverseComplementReverse.toString()<<endl;
#endif
			bool currentFlag=false;
			if(_vecSeq_R2.size()>0) //there is R2, in this case, we assume the R2 is 3' first. we can simply do the 5' mapping, no revcomp
			{
					currentFlag=alignConstant(&(_vecSeq_R2.at(i)), &(_vecIsotype.at(j)), _sm, _gapOpen, _gapExtension
							,1, gapModel, FivePrime,
							_matchRateThreshold, _minimumOverlapLength, _offset, 
				 //output now
				  &best3PrimeAlign,   &best3PrimeScore );
			}
			else{
					currentFlag=alignConstant(&(_vecSeq.at(i)), &reverseComplementReverse, _sm, _gapOpen, _gapExtension
							,1, gapModel, ThreePrime,
							_matchRateThreshold, _minimumOverlapLength, _offset, 
				 //output now
				  &best3PrimeAlign,   &best3PrimeScore );
			}
		  if(currentFlag)
		  {
			  found3PrimeFlag=currentFlag;
			  best3PrimeIndex=j;
		  }
	    }//for 3 prime mapping
	}//loop through the isotype sequences 
      
      //cout<<"Done with mapping:i="<<i<<endl;
      //###done with mapping, now we need to check out constant region
      

      //for constant region, we only need to get the constant alignment and will write where the constant breaks
      //if there are breaks within the constant region. we are anticipating more to be full length???
      
      //first, we need to check the alignment, ?????local alignment or overlap alignment????? 
      //  **********************************       ---sequence read
      //    ==^^^^^^^^^^^^===                      ---foward primer
      //             =====!!!!!!!!!!!!!==          ----reverse primer
      //
      //= : mismatch
      //* : sequence read
      //^ : forward primer match
      //! : reverse primer match
      
      //so the algorithm is we check the alignment, write down where the constant alignment ends
      //we also want to distinguish between ends within or outside the constant region
      
      //now we want to write down data by each different allele, we can always manually combine them into different categories.

      //********WE NEED THIS LATER**************************************
      /*if(!foundReverseFlag||!foundForwardFlag)
	{
	  continue;
	}
      */
      
      //write down the stats, for where the alignment breaks (the end position).
      //write to the variable of pt_vec_map_pos_end 
      //first look up what this sequence is from IgM or IgG 
      //isotype
      //Here we only list two cases where we got alignment, it could
	 //be the case that we don't have alignemtn, in which case, we simply skip the stats
      //unsigned int foundIndex=LookUpVectorIndex(_vecForward.at(bestForwardIndex).GetName(), g_vec_primer_isotype);
      unsigned int foundIndex=best5PrimeIndex;//here we don't need to look up, because we did not precess adaptor/primer, we simply get sort them by each allel
      if(found3PrimeFlag)
	{
	  foundIndex=best3PrimeIndex;
	}

      //here, to keep track of positions, because we are running local alignment, so we can simply keep 
	  //track the end position of alignment.
      //we also have to check whether the break is in the constant
      if(found5PrimeFlag)
	{	   
	  pt_vec_map_pos_end[foundIndex].push_back(best5PrimeAlign.GetSubjectIndexEnd());
	  pt_vec_map_pos_start[foundIndex].push_back(best5PrimeAlign.GetSubjectIndexStart());

	  pt_vec_map_pos_end_seq[foundIndex].push_back(best5PrimeAlign.GetPatternIndexEnd());
	  pt_vec_map_pos_start_seq[foundIndex].push_back(best5PrimeAlign.GetPatternIndexStart());
	  
	  strPattern=best5PrimeAlign.GetPattern(true);
	  strSubject=best5PrimeAlign.GetSubject(true);
	  match_rate=1-CompareStrings(strPattern, strSubject)/((double)strPattern.length());
	  
	  pt_vec_map_match_rate[foundIndex].push_back(match_rate);
	  pt_vec_map_overlap[foundIndex].push_back(strPattern.length());
	   
	  pt_vec_map_seq_len[foundIndex].push_back(_vecSeq.at(i).GetLength());
	  if(_demux)
	  {
	    pt_vec_demux[foundIndex].push_back(_vecSeq.at(i));
		if(_vecSeq_Q.size()>0)
			pt_vec_demux_Q[foundIndex].push_back(_vecSeq_Q.at(i));
		if(_vecSeq_R2.size()>0)
		{
			pt_vec_demux_R2[foundIndex].push_back(_vecSeq_R2.at(i));
			if(_vecSeq_Q_R2.size()>0)
			{
				pt_vec_demux_Q_R2[foundIndex].push_back(_vecSeq_Q_R2.at(i));
			}
		}
	  }
	}
	 
      if(found3PrimeFlag)
	{ 

		if(_vecSeq_R2.size()==0) //here we do real map3 end, otherwise we do 5' mapping on R2
		{
			  //this part, we need to test, to see whether this is correct
			  //the rationale is that we revcomp the alignment, then the break part for this revComp is that we get the beginning of the alignment and then used the total length (-1) to get the correct break position (end pos)
			  pt_vec_map_pos_end[foundIndex].push_back(_vecIsotype.at(foundIndex).GetLength()-1 - best3PrimeAlign.GetSubjectIndexStart());
			  pt_vec_map_pos_start[foundIndex].push_back(_vecIsotype.at(foundIndex).GetLength()-1 -best3PrimeAlign.GetSubjectIndexEnd());

			  pt_vec_map_pos_end_seq[foundIndex].push_back(best3PrimeAlign.GetPatternIndexEnd());
			  pt_vec_map_pos_start_seq[foundIndex].push_back(best3PrimeAlign.GetPatternIndexStart());
			  
			  strPattern=best3PrimeAlign.GetPattern(true);
			  strSubject=best3PrimeAlign.GetSubject(true);
			  match_rate=1-CompareStrings(strPattern, strSubject)/((double)strPattern.length());
			  
			  pt_vec_map_match_rate[foundIndex].push_back(match_rate);
			  pt_vec_map_overlap[foundIndex].push_back(strPattern.length());

			  pt_vec_map_seq_len[foundIndex].push_back(_vecSeq.at(i).GetLength());
		}
		else { // here we are not doing map 3, but instead do map 5 on R2
			pt_vec_map_pos_end[foundIndex].push_back( best3PrimeAlign.GetSubjectIndexEnd());
			  pt_vec_map_pos_start[foundIndex].push_back(best3PrimeAlign.GetSubjectIndexStart());

			  pt_vec_map_pos_end_seq[foundIndex].push_back(best3PrimeAlign.GetPatternIndexEnd());
			  pt_vec_map_pos_start_seq[foundIndex].push_back(best3PrimeAlign.GetPatternIndexStart());
			  
			  strPattern=best3PrimeAlign.GetPattern(true);
			  strSubject=best3PrimeAlign.GetSubject(true);
			  match_rate=1-CompareStrings(strPattern, strSubject)/((double)strPattern.length());
			  
			  pt_vec_map_match_rate[foundIndex].push_back(match_rate);
			  pt_vec_map_overlap[foundIndex].push_back(strPattern.length());

			  pt_vec_map_seq_len[foundIndex].push_back(_vecSeq_R2.at(i).GetLength());
		}
	  if(_demux)
	  {
	    pt_vec_demux[foundIndex].push_back(_vecSeq.at(i));
		if(pt_vec_demux_Q!=NULL)
			pt_vec_demux_Q[foundIndex].push_back(_vecSeq_Q.at(i));
		if(pt_vec_demux_R2!=NULL)
			pt_vec_demux_R2[foundIndex].push_back(_vecSeq_R2.at(i));
		if(pt_vec_demux_Q_R2!=NULL)
			pt_vec_demux_Q_R2[foundIndex].push_back(_vecSeq_Q_R2.at(i));
	  }
	}
      

      //now we are done with stats recording.
      
      //***************************output*****************
      //prepare the aligned for debugging purpose
      //#the forward 
      string leadingSpaceOriginal(""); //for sequence
      string leadingSpace5Prime(""); //for 5' mapping isotype
      string leadingSpace3Prime(""); //for 3' mapping isotype
      string replaceOne;
      SequenceString tempLst5Prime;  //for 5' isotype
      SequenceString tempLst3Prime;  //for 3' isotype
      SequenceString tempLstSeq;  //for sequence 
            
      //cout<<"forwardSet Output read"<<endl;
      if(found5PrimeFlag)
	{
	  //###we need to figure out the leading spaces in front of the original sequences
	  unsigned int startOriginal=best5PrimeAlign.GetPatternIndexStart();
	  unsigned int startSubject=best5PrimeAlign.GetSubjectIndexStart();
	  if(startSubject>startOriginal)
	    {
	      //leadingSpaceOriginalsg=as.character();
	      for(unsigned int p=0;p<startSubject-startOriginal;p++)
			{
			  leadingSpaceOriginal.push_back('-');
			}
	    }
	  else
	    {
	      for(unsigned int p=0;p<startOriginal-startSubject;p++)
			{
			  leadingSpace5Prime.push_back('-');
			}
	    }
	  
	  //now we need to add the aligned sequence to replace the original one
	  //#we assume the original one is longer than the aligned one, it has to be
	  // #this is the intial part
	  replaceOne=_vecIsotype.at(best5PrimeIndex).GetSequence().substr(0, best5PrimeAlign.GetSubjectIndexStart());
	  //#aligned part
	  replaceOne.append( best5PrimeAlign.GetSubject(true));
	  //#last part unaligned
	  replaceOne.append(_vecIsotype.at(best5PrimeIndex).GetSequence().substr(best5PrimeAlign.GetSubjectIndexEnd()+1)
			    );
	  	  
	  tempLst5Prime.SetSequence(leadingSpace5Prime+replaceOne);
	  tempLst5Prime.SetName(_vecIsotype.at(best5PrimeIndex).GetName());
	}
    else  //this is not aligned, might be because we do 3' or no alignment possible.
	{
	  //#no need to add leading space
	  tempLst5Prime.SetName("NoMatch");
	  tempLst5Prime.SetSequence("");
	}
      //############here to do!!!!!!!!!
                  
      //#the reverse
      //cout<<"reverse set output read"<<endl;
	if(found3PrimeFlag)  //working on the isotype part, not the sequence part
	  {
	    //#now we need to add the aligned sequence to replace the original one
	    //#we assume the original one is longer than the aligned one, it has to be
	    //	 #this is the intial part
		
	    SequenceString rc3PrimeSeq=ReverseComplement(_vecIsotype.at(best3PrimeIndex));
		if(_vecSeq_R2.size()>0) //with separated R2. we should map 5' of R2, no need to do revcomp of isotype
		{
			rc3PrimeSeq=_vecIsotype.at(best3PrimeIndex);
		}
	    replaceOne=rc3PrimeSeq.GetSequence().substr(0, best3PrimeAlign.GetSubjectIndexStart());
	    //cout<<"\t\t*****subject start:"<<bestReverseAlign.GetSubjectIndexStart()<<";subject end:"<<bestReverseAlign.GetSubjectIndexEnd()<<endl;
	    //cout<<"\t\t***replaceOne:"<<replaceOne<<endl;
	    //#aligned part
	    replaceOne.append(best3PrimeAlign.GetSubject(true));
	    //cout<<"\t\t***replaceOne22:"<<replaceOne<<endl;
	    //#last part unaligned
	    replaceOne.append( rc3PrimeSeq.GetSequence().substr(best3PrimeAlign.GetSubjectIndexEnd()+1, rc3PrimeSeq.GetLength()));
	    //cout<<"\t\t***replaceOne333:"<<replaceOne<<endl;
	    tempLst3Prime.SetName(_vecIsotype.at(best3PrimeIndex).GetName());
		
	    //tempLstR.SetSequence(replaceOne);
	    //#now we need to figure out how the leading space to put in front of reverse one
	    if(best3PrimeAlign.GetPatternIndexStart() >= best3PrimeAlign.GetSubjectIndexStart())
	      {
				unsigned int templen=best3PrimeAlign.GetPatternIndexStart()- best3PrimeAlign.GetSubjectIndexStart();
				for(unsigned int p=0;p<templen;p++)
				  {
					leadingSpace3Prime.push_back('-');
				  }
				tempLst3Prime.SetSequence(leadingSpace3Prime+replaceOne);
			}
	    else
	      {
				unsigned int templen=best3PrimeAlign.GetSubjectIndexStart()-best3PrimeAlign.GetPatternIndexStart() ;
				for(unsigned int p=0;p<templen;p++)
				{
				  leadingSpaceOriginal.push_back('-');
				}
				tempLst3Prime.SetSequence(replaceOne);
		//#here, in this case, the adaptor+primer is longer than the seqs just by alignment, then we need to simply remove some leading part of the adaptor primer
		//tempLst3Prime.SetSequence(replaceOne.substr( best3PrimeAlign.GetSubjectIndexStart()-best3PrimeAlign.GetPatternIndexStart(),
		//					     replaceOne.length()));
	      }
	  }
	else
	  {
	    tempLst3Prime.SetSequence("");
	    tempLst3Prime.SetName("NoMatch");
	  }
	//cout<<"end of reverse on, tempLstR.seq:"<<tempLstR.toString()<<endl;

	//#now we need to take care of the read sequence alignment string
	//#on the reverse part first
	//cout<<"seq string output read...."<<endl;
	//unsigned int spaceCarryOverFTR=0;//####this is the leading space for reverse adaptor primer, because the insertion in the forward alignment
	if(type==ThreePrime && _vecSeq_R2.size()>0)
	{
		//cout<<"in here 3 prime and R2 avail"<<endl;
		tempLstSeq.SetSequence(_vecSeq_R2.at(i).GetSequence());
			tempLstSeq.SetName(_vecSeq_R2.at(i).GetName());
			replaceOne=_vecSeq_R2.at(i).GetSequence();
	}
	else {
			//cout<<"set one for no 3' and R2 available"<<endl;
			tempLstSeq.SetSequence(_vecSeq.at(i).GetSequence());
			tempLstSeq.SetName(_vecSeq.at(i).GetName());
			replaceOne=_vecSeq.at(i).GetSequence();
	}
	if(found3PrimeFlag)
	{
	  //cout<<"1.."<<endl;
	  replaceOne=replaceOne.substr(0, best3PrimeAlign.GetPatternIndexStart());
	  //#aligned part
	  //cout<<"2.."<<endl;
	  replaceOne.append( best3PrimeAlign.GetPattern(true));
	  //#last part unaligned

	  //cout<<"3..."<<endl;
	  replaceOne.append( tempLstSeq.GetSequence().substr(best3PrimeAlign.GetPatternIndexEnd()+1));//nchar(ff[[i]])), sep="");
	  tempLstSeq.SetSequence(replaceOne);
	  //tempLstSeq.SetName(_vecSeq.at(i).GetName());
	}
	if(found5PrimeFlag)
	  {
	    //cout<<"1.."<<endl;
	    replaceOne=replaceOne.substr(0, best5PrimeAlign.GetPatternIndexStart());
	    //#aligned part
	    //cout<<"2.."<<endl;
	    replaceOne.append( best5PrimeAlign.GetPattern(true));
	    //#last part unaligned
	    //cout<<"3.."<<endl;
	    //cout<<"\ttempLstSeq.GetSequence().length():"<<tempLstSeq.GetSequence().length()<<endl;
	    //cout<<"\tbestForwardAlign.GetPatternIndexEnd():"<<bestForwardAlign.GetPatternIndexEnd()<<endl;
	    replaceOne.append( tempLstSeq.GetSequence().substr(best5PrimeAlign.GetPatternIndexEnd()+1));// nchar(ff[[i]])));
	    //cout<<"4.."<<endl;
	    tempLstSeq.SetSequence(replaceOne);
	    //spaceCarryOverFTR=bestForwardAlign.GetPattern(true).length()- bestForwardAlign.GetPattern(false).length();

	    //cout<<"***with gap:"<<bestForwardAlign.GetPattern(true)<<";without gap:"<<bestForwardAlign.GetPattern(false)<<endl;
	    //cout<<"length is "<<bestForwardAlign.GetPattern(true).length()<<":"<<bestForwardAlign.GetPattern(false).length()<<endl;
	    //cout<<"&&&&&&&&&&&&&carry over FTPR is :"<<spaceCarryOverFTR<<endl;
	  }
	//in case there is no matching on either side, we do nothing, since tempLstSeq has been assigned in the beginning
	
	//string tempStrSpaces;
	
	//for(unsigned int m = 0 ; m < spaceCarryOverFTR; m++)
	//  {
	//    tempStrSpaces.append("-");
	//  }
	//tempLstR.SetSequence(leadingSpaceOriginal+tempLstR.GetSequence());//<-paste(tempStr, as.character(tempLstR$seq), sep="");
	//tempLstR$seq<-paste(leadingSpaceOriginal, tempStr, sep="");
	tempLstSeq.SetSequence(leadingSpaceOriginal+tempLstSeq.GetSequence());
	//cout<<"****************tempLstSeq ----:"<<tempLstSeq.GetName()<<"||"<<tempLstSeq.GetSequence()<<endl;
	//******************write to vectors***************
	//#now put the sequences to the correct vectors
	//cout<<"\tready to output strings........"<<endl;
	vector<SequenceString> * p_vec_map;
       /*if(foundCrossOverFlag)
	  {
	    p_vec_map=&vec_mapCrossOver;
	    stringstream ss;
	    ss<<tempLstF.GetName()<<":"<<bestForwardAlign.GetSubjectIndexStart()<<":"
	      <<bestForwardAlign.GetSubjectIndexEnd()
	      <<":"<<bestReverseAlign.GetSubjectIndexStart()
	      <<":"<<bestReverseAlign.GetSubjectIndexEnd();
	    tempLstF.SetName(ss.str());
	  }
	else  //not crossover flag
	  {
	    if(foundBreakOutsideConstantFlag)
	      {
		p_vec_map=&vec_mapBreakOutsideConstant;
		stringstream ss;
		ss<<tempLstF.GetName()<<":"<<bestForwardAlign.GetSubjectIndexStart()<<":"
		  <<bestForwardAlign.GetSubjectIndexEnd();

		tempLstF.SetName(ss.str());
	      }
	    else  //'if' foundBreakOutsideConstantFlag
	    {
		if(foundForwardFlag)
		  {
		    unsigned int found=LookUpVectorIndex(tempLstF.GetName(), _vecForward);
		    if(foundReverseFlag)
		      {

			//check to see whether we need to store the data by isotype
			//if(g_by_isotype_flag)
			//	{
			
			//cout<<"\t**********writing to vec now"<<found<<endl;

			//here we always store by isotype
			if(((signed)found) != -1)//string::npos)
			  {
			    p_vec_map=&pt_vec_mapBoth[found];  //g_vec_mapBoth is initilized in the processingAdaptor function above
			  }
			else
			  {
			    cout<<"***WARNING:Can not find the forward primer name, something is wrong"<<endl;
			    p_vec_map=&pt_vec_mapBoth[_vecForward.size()];
			  }
			stringstream ss;
			ss<<tempLstF.GetName()<<":"<<bestForwardAlign.GetSubjectIndexStart()<<":"
			  <<bestForwardAlign.GetSubjectIndexEnd()
			  <<":"<<bestReverseAlign.GetSubjectIndexStart()
			  <<":"<<bestReverseAlign.GetSubjectIndexEnd();
			tempLstF.SetName(ss.str());						
		      }
		    else
		      {
			if(((signed)found) != -1)//string::npos)
			  {
			    p_vec_map=&pt_vec_mapForward[found];  //g_vec_mapBoth is initilized in the processingAdaptor function above
			  }
			else
			  {
			    cout<<"***ERROR:Can not find the forward primer name, something is wrong"<<endl;
			    p_vec_map=&pt_vec_mapForward[_vecForward.size()];
			  }
			stringstream ss;
			ss<<tempLstF.GetName()<<":"<<bestForwardAlign.GetSubjectIndexStart()<<":"
			  <<bestForwardAlign.GetSubjectIndexEnd();
			
			tempLstF.SetName(ss.str());
		      }
		  }
		else
		  {
		    if((foundReverseFlag))
		      {
			p_vec_map=&vec_mapReverse;
		      }
		      else
			p_vec_map=&vec_mapReverse;
		  }
	      }
	  }
	  
	*/
        p_vec_map=&vec_mapNone;
	if(found5PrimeFlag)
	  {
	    unsigned int found=LookUpVectorIndex(tempLst5Prime.GetName(), _vecIsotype);
	    //check to see whether we need to store the data by isotype
	    //if(g_by_isotype_flag)
	    //	{
#ifdef debug	    
	    cout<<"\t**********writing to vec now"<<found<<endl;
#endif	    
	    //here we always store by isotype
	    if(((signed)found) != -1)//string::npos)
	      {
		p_vec_map=&pt_vec_mapped[found];  //g_vec_mapBoth is initilized in the processingAdaptor function above
	      }
	    else
	      {
		cout<<"***WARNING:Can not find the forward primer name, something is wrong"<<endl;
		p_vec_map=&pt_vec_mapped[_vecIsotype.size()];
	      }
	    stringstream ss;
	    ss<<tempLst5Prime.GetName()<<":"<<best5PrimeAlign.GetSubjectIndexStart()<<":"
	      <<best5PrimeAlign.GetSubjectIndexEnd();
	    
	    tempLst5Prime.SetName(ss.str());						
	  }
	
	if(found3PrimeFlag)
	  {
	    unsigned int found=LookUpVectorIndex(tempLst3Prime.GetName(), _vecIsotype);
	    if(((signed)found) != -1)//string::npos)
	      {
		p_vec_map=&pt_vec_mapped[found];  //g_vec_mapBoth is initilized in the processingAdaptor function above
	      }
	    else
	      {
		cout<<"***ERROR:Can not find the forward primer name, something is wrong"<<endl;
		p_vec_map=&pt_vec_mapped[_vecIsotype.size()];
	      }
	    stringstream ss;
	    ss<<tempLst3Prime.GetName()<<":"<<best3PrimeAlign.GetSubjectIndexStart()<<":"
	      <<best3PrimeAlign.GetSubjectIndexEnd();
	    
	    tempLst3Prime.SetName(ss.str());
	  }
	
	//now we need to store the data into the correct vectors/files
	p_vec_map->push_back(tempLstSeq);
	//cout<<"add tempLstSeq"<<endl;
	//this following is to write the fastq or R2 if possible. 
	if(!found5PrimeFlag && !found3PrimeFlag)
	{
		if(_vecSeq_R2.size()>0 )
		{
			if(type==FivePrime)
			{//R1 is R1 and R2 is R2
				cout<<"map none: add Five Primer"<<endl;
				vec_mapNone_R2.push_back(_vecSeq_R2.at(i));
				if(_vecSeq_Q_R2.size()>0)
					vec_mapNone_Q_R2.push_back(_vecSeq_Q_R2.at(i));
				
				if(_vecSeq_Q.size()>0)
				{
					vec_mapNone_Q.push_back(_vecSeq_Q.at(i));
				}
			}
			else  //for 3' prime map and R1 R2 
			{//R1 is R2 and R2 is R1, meaning vec_mapNone is taking in R2, now this vec_mapNone_R2 holding R1 now
					if(_vecSeq_Q.size()>0)
					{
						vec_mapNone_Q.push_back(_vecSeq_Q_R2.at(i));
					}
					
					vec_mapNone_R2.push_back(_vecSeq.at(i));
					if(_vecSeq_Q_R2.size()>0)
						vec_mapNone_Q_R2.push_back(_vecSeq_Q.at(i));
			}
			
		}
		else{ //no R2 
				if(_vecSeq_Q.size()>0)
				{
					vec_mapNone_Q.push_back(_vecSeq_Q.at(i));
				}
		}
	}
	
	if(found5PrimeFlag)
	  p_vec_map->push_back(tempLst5Prime);
	if(found3PrimeFlag)
	  p_vec_map->push_back(tempLst3Prime);
		
	//cout<<"start writing output........"<<endl;
	//need to figure out how to (trunc or app) write the output at each round.
	
	if(i>=timeOfWriting*numOfSeqsUnit) //#write once every 20000 sequences
	  {
	    
	    //if(timeOfWriting>1)
	    //  {
	    //	mode=ofstream::app;
	    //  }
	    timeOfWriting++;
	    string t_fileName;
	    //cout<<"i round:"<<i<<endl;
	    //mapped
	    vector<string> header;
	    header.push_back("lengthAlign");
	    header.push_back("matchRate");
	    header.push_back("seq_start");
	    header.push_back("seq_end");
	    header.push_back("isotype_start");
	    header.push_back("isotype_end");
	    header.push_back("Seq_length");

	    bool writeHeader=true;

	    vector<vector<double> > vec_stats(7);
	    for(unsigned int s=0;s<=_vecIsotype.size();s++)
	      {
		if(pt_vec_mapped[s].size()>0)
		  {
		    if(numOfWritesDone_mapped[s]==0)
		      {
			mode=ofstream::trunc;
			writeHeader=true;
		      }
		    else
		      {
			mode=ofstream::app;
			writeHeader=false;
		      }
		    //cout<<"\t------writing files at i-----:"<<s<<endl;
		    //cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
		    if(s<_vecIsotype.size())
		      {
			t_fileName=_R1_fname+ _vecIsotype.at(s).GetName()+".fasta";

			vector<double> temp(pt_vec_map_overlap[s].begin(),pt_vec_map_overlap[s].end());
			vec_stats[0]= temp;
			temp.clear();

			vec_stats[1]=pt_vec_map_match_rate[s];

			vector<double> temp1(pt_vec_map_pos_start_seq[s].begin(),pt_vec_map_pos_start_seq[s].end());
			vec_stats[2]=temp1;
			temp1.clear();

			vector<double> temp2(pt_vec_map_pos_end_seq[s].begin(),pt_vec_map_pos_end_seq[s].end());
			vec_stats[3]=temp2;
			temp2.clear();

			vector<double> temp3(pt_vec_map_pos_start[s].begin(),pt_vec_map_pos_start[s].end());
			vec_stats[4]=temp3;
			temp3.clear();

			vector<double> temp4(pt_vec_map_pos_end[s].begin(), pt_vec_map_pos_end[s].end());
			vec_stats[5]=temp4;
			temp4.clear();

			vector<double> temp5(pt_vec_map_seq_len[s].begin(), pt_vec_map_seq_len[s].end());
			vec_stats[6]=temp5;
			temp5.clear();
			
			WriteTextTableFile(t_fileName+"_ConstantStat.txt", vec_stats, ' ', writeHeader,mode, header);
			
			pt_vec_map_overlap[s].clear();
			pt_vec_map_match_rate[s].clear();
			pt_vec_map_pos_start_seq[s].clear();
			pt_vec_map_pos_end_seq[s].clear();
			pt_vec_map_pos_start[s].clear();
			pt_vec_map_pos_end[s].clear();
			pt_vec_map_seq_len[s].clear();
		      }
		    else
		      {
			t_fileName=_R1_fname+ "notFoundIsotype.fasta";
		      }
		    WriteFasta(t_fileName, pt_vec_mapped[s],100, mode);
		    if(_demux)
		      {
				  if(_vecSeq_Q.size()<=0)
				  {
						t_fileName=_R1_fname+ _vecIsotype.at(s).GetName()+"_demux.fasta";
						WriteFasta(t_fileName, pt_vec_demux[s],100, mode);
						pt_vec_demux[s].clear();
				  }else
				  {
						t_fileName=_R1_fname+ _vecIsotype.at(s).GetName()+"_demux.fastq";
						WriteFastq(t_fileName, pt_vec_demux[s],pt_vec_demux_Q[s], mode);
						pt_vec_demux[s].clear();
						pt_vec_demux_Q[s].clear();
				  }
				  //Read2 
				  if(_vecSeq_R2.size()>0)
				  {
						if(_vecSeq_Q_R2.size()<=0)
						  {
								t_fileName=_R2_fname+ _vecIsotype.at(s).GetName()+"_demux.fasta";
								WriteFasta(t_fileName, pt_vec_demux_R2[s],100, mode);
								pt_vec_demux_R2[s].clear();
						  }else   //fastq
						  {
								t_fileName=_R2_fname+ _vecIsotype.at(s).GetName()+"_demux.fastq";
								WriteFastq(t_fileName, pt_vec_demux_R2[s],pt_vec_demux_Q_R2[s], mode);
								pt_vec_demux_R2[s].clear();
								pt_vec_demux_Q_R2[s].clear();
						  }
				  }
				  
		      }
		    //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
		    pt_vec_mapped[s].clear();
		    //cout<<"vecBoth cleared:"<<vec_mapBoth.size()<<endl;
		    //WriteTextFile(t_fileName+"_ConstantStat.txt", pt_vec_map_pos_end[s], ' ', 1,mode);
		    
		    numOfWritesDone_mapped[s]++;
		    //pt_vec_map_pos_end[s].clear();//g_vec_len_mapBoth.at(s).clear();
		    
		    /*if(g_trim_flag)
		      {
		      WriteFasta(t_fileName+"_trim.fas", g_vec_mapBoth_trim.at(s), 100,ofstream::app);
		      g_vec_mapBoth_trim.at(s).clear();
		      }*/
		  }
	      }
	    
	    //none
	    if(vec_mapNone.size()>0)
	      {
				//cout<<"\t------writing files at i-----:"<<s<<endl;
				//cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
				if(numOfWritesDone_unmapped==0)
				  {
					mode=ofstream::trunc;
				  }
				else
				  {
					mode=ofstream::app;
				  }
				  if(type==ThreePrime&&_vecSeq_R2.size()>0) //R2 is R1 and R1 is R2
				  {
					  if(_vecSeq_Q_R2.size()>0)
					  {
						t_fileName=_R2_fname+"_MapNone.fastq";
						WriteFastq(t_fileName, vec_mapNone, vec_mapNone_Q,  mode);
						
						t_fileName=_R1_fname+"_MapNone.fastq";
						WriteFastq(t_fileName, vec_mapNone_R2, vec_mapNone_Q_R2,  mode);
						vec_mapNone_Q_R2.clear();
						vec_mapNone_R2.clear();
						vec_mapNone_Q.clear();
					  }
					  else  //fasta R2 and R1
					  {
						  t_fileName=_R2_fname+"_MapNone.fasta";
						WriteFasta(t_fileName, vec_mapNone,100, mode);
						t_fileName=_R1_fname+"_MapNone.fasta";
						WriteFasta(t_fileName, vec_mapNone_R2,100, mode);
						vec_mapNone_R2.clear();
					  }
						
				  }
				  else  //either no R2 or not Three primer, so R1 is R1 and R2 is R2
				  {
					  if(_vecSeq_Q.size()>0)
					  {
						t_fileName=_R1_fname+"_MapNone.fastq";
						WriteFastq(t_fileName, vec_mapNone, vec_mapNone_Q,  mode);
						vec_mapNone_Q.clear();
					  }
					  else
					  {
						  t_fileName=_R1_fname+"_MapNone.fasta";
						  WriteFasta(t_fileName, vec_mapNone,100, mode);
						  
					  }
					  
					  if(_vecSeq_R2.size()>0)
					  {
						  if(_vecSeq_Q_R2.size()>0)
						  {
								t_fileName=_R2_fname+"_MapNone.fastq";
								WriteFastq(t_fileName, vec_mapNone_R2, vec_mapNone_Q_R2,  mode);
								vec_mapNone_Q_R2.clear();
								vec_mapNone_R2.clear();
						  }
						  else
						  {
							  t_fileName=_R2_fname+"_MapNone.fasta";
								WriteFasta(t_fileName, vec_mapNone_R2,100, mode);
								vec_mapNone_R2.clear();
						  }
					  }
				  }
				  numOfWritesDone_unmapped++;
				  vec_mapNone.clear();
	      }
	    
	    //cout<<"done map both"<<endl;
#ifdef debug
		flush(cout);
#endif
	  }//end of each 1000 seqs read write
      
    }//end of for loop of sequence data vec

  //one last writing
 // cout<<"start writing last round of output........"<<endl;
    
  string t_fileName;
  //mapBoth
  vector<string> header;
  header.push_back("lengthAlign");
  header.push_back("matchRate");
  header.push_back("seq_start");
  header.push_back("seq_end");
  header.push_back("isotype_start");
  header.push_back("isotype_end");
  header.push_back("Seq_length");

  bool writeHeader=true;
  
  vector<vector<double> > vec_stats(7);
  for(unsigned int s=0;s<_vecIsotype.size();s++)
    {
      if(pt_vec_mapped[s].size()>0)
	{
	  //cout<<"\t------writing files at i-----:"<<s<<endl;
	  //cout<<"numOfWritesDone:"<<numOfWritesDone_mapped[s]<<endl;
	  if(numOfWritesDone_mapped[s]==0)
	    {
	      mode=ofstream::trunc;
	      writeHeader=true;
	    }
	  else
	    {
	      mode=ofstream::app;
	      writeHeader=false;
	    }
	  
	  if(s<_vecIsotype.size())
	    {
	      t_fileName=_R1_fname+ _vecIsotype.at(s).GetName()+".fasta";
	      vector<double> temp(pt_vec_map_overlap[s].begin(),pt_vec_map_overlap[s].end());
	      vec_stats[0]= temp;
	      temp.clear();

	      vec_stats[1]=pt_vec_map_match_rate[s];

	      vector<double> temp1(pt_vec_map_pos_start_seq[s].begin(),pt_vec_map_pos_start_seq[s].end());
	      vec_stats[2]=temp1;
	      temp1.clear();

	      vector<double> temp2(pt_vec_map_pos_end_seq[s].begin(),pt_vec_map_pos_end_seq[s].end());
	      vec_stats[3]=temp2;
	      temp2.clear();

	      vector<double> temp3(pt_vec_map_pos_start[s].begin(),pt_vec_map_pos_start[s].end());
	      vec_stats[4]=temp3;
	      temp3.clear();

	      vector<double> temp4(pt_vec_map_pos_end[s].begin(), pt_vec_map_pos_end[s].end());
	      vec_stats[5]=temp4;
	      temp4.clear();

	      vector<double> temp5(pt_vec_map_seq_len[s].begin(), pt_vec_map_seq_len[s].end());
	      vec_stats[6]=temp5;
	      temp5.clear();
	      WriteTextTableFile(t_fileName+"_ConstantStat.txt", vec_stats, ' ', writeHeader,mode, header);
			
	      pt_vec_map_overlap[s].clear();
	      pt_vec_map_match_rate[s].clear();
	      pt_vec_map_pos_start_seq[s].clear();
	      pt_vec_map_pos_end_seq[s].clear();
	      pt_vec_map_pos_start[s].clear();
	      pt_vec_map_pos_end[s].clear();
	      pt_vec_map_seq_len[s].clear();
	    }
	  else
	    {
	      t_fileName=_R1_fname+ "notFoundIsotype.fasta";
	    }
	  //cout<<"writing fasta"<<endl;
	  // t_fileName=_mapBoth_fname+ _vecForward.at(s).GetName();
	  WriteFasta(t_fileName, pt_vec_mapped[s],100, mode);
	  if(_demux)
	    {
			/*cout<<"===========Demuxing............."<<endl;
			cout<<"\tisotype index:"<<s<<endl;
			cout<<"\tvecSeq_Q.size():"<<_vecSeq_Q.size()<<endl;
			cout<<"\tvecSeq_R2.size():"<<_vecSeq_R2.size()<<endl;
			cout<<"\tvecSeq_Q_R2.size():"<<_vecSeq_Q_R2.size()<<endl;
			*/
	      if(_vecSeq_Q.size()<=0)
				  {
					  //cout<<"-----write R1 fasta"<<endl;
						t_fileName=_R1_fname+ _vecIsotype.at(s).GetName()+"_demux.fasta";
						WriteFasta(t_fileName, pt_vec_demux[s],100, mode);
						pt_vec_demux[s].clear();
				  }else
				  {
						//cout<<"-----write R1 fastq"<<endl;
						t_fileName=_R1_fname+ _vecIsotype.at(s).GetName()+"_demux.fastq";
						WriteFastq(t_fileName, pt_vec_demux[s],pt_vec_demux_Q[s], mode);
						pt_vec_demux[s].clear();
						pt_vec_demux_Q[s].clear();
				  }
				  //Read2 
				  if(_vecSeq_R2.size()>0)
				  {
						if(_vecSeq_Q_R2.size()<=0)
						  {
							  //cout<<"-----write R2 fasta"<<endl;
								t_fileName=_R2_fname+ _vecIsotype.at(s).GetName()+"_demux.fasta";
								WriteFasta(t_fileName, pt_vec_demux_R2[s],100, mode);
								pt_vec_demux[s].clear();
						  }else   //fastq
						  {
								//cout<<"-----write R2 fastq"<<endl;
								t_fileName=_R2_fname+ _vecIsotype.at(s).GetName()+"_demux.fastq";
								WriteFastq(t_fileName, pt_vec_demux_R2[s],pt_vec_demux_Q_R2[s], mode);
								pt_vec_demux_R2[s].clear();
								pt_vec_demux_Q_R2[s].clear();
						  }
				  }

		  //t_fileName=_map_fname+ _vecIsotype.at(s).GetName()+"_demux.fasta";
	      //WriteFasta(t_fileName, pt_vec_demux[s],100, mode);
	      //pt_vec_demux[s].clear();
	    }
	  //cout<<"writing info"<<endl;
	  //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
	  pt_vec_mapped[s].clear();
	  pt_vec_demux[s].clear();
	  //WriteTextFile(t_fileName+"_ConstantStat.txt", pt_vec_map_pos_end[s], ' ', 1,mode);
	  //pt_vec_map_pos_end[s].clear();//g_vec_len_mapBoth.at(s).clear();
	  //cout<<"vecBoth cleared:"<<pt_vec_mapped[s].size()<<endl;
	  //WriteTextFile(t_fileName+"_primerdimerStat.txt", pt_vec_isotype_primerDimer_pos_cross[s], ' ', 1,ofstream::app);
	  //pt_vec_isotype_primerDimer_pos_cross[s].clear();//g_vec_len_mapBoth.at(s).clear();
	  /*if(g_trim_flag)
	    {
	    WriteFasta(t_fileName+"_trim.fas", g_vec_mapBoth_trim.at(s), 100,ofstream::app);
	    g_vec_mapBoth_trim.at(s).clear();
	    }*/
	}
    }
  
  //none
  if(vec_mapNone.size()>0)
    {
      //cout<<"\t------writing files at unmapped last-----:"<<endl;
      //cout<<"vecBoth:"<<vec_mapBoth.size()<<endl;
      if(numOfWritesDone_unmapped==0)
	{
	  mode=ofstream::trunc;
	}
      else
	{
	  mode=ofstream::app;
	}
	if(type==ThreePrime&&_vecSeq_R2.size()>0)//R1 is R2 and R2 is R1
		  {
			  if(_vecSeq_Q_R2.size()>0) 
			  {
				t_fileName=_R2_fname+"_MapNone.fastq";
				WriteFastq(t_fileName, vec_mapNone, vec_mapNone_Q,  mode);
				
				t_fileName=_R1_fname+"_MapNone.fastq";
				WriteFastq(t_fileName, vec_mapNone_R2, vec_mapNone_Q_R2,  mode);
				vec_mapNone_Q_R2.clear();
				vec_mapNone_R2.clear();
				vec_mapNone_Q.clear();
			  }
			  else  //fasta R2 and R1
			  {
				  t_fileName=_R2_fname+"_MapNone.fasta";
				WriteFasta(t_fileName, vec_mapNone,100, mode);
				t_fileName=_R1_fname+"_MapNone.fasta";
				WriteFasta(t_fileName, vec_mapNone_R2,100, mode);
				vec_mapNone_R2.clear();
			  }
				
		  }
		  else  //either no R2 or not Three primer, so R1 is R1 and R2 is R2
		  {
			  if(_vecSeq_Q.size()>0)
			  {
				  //cout<<"Write fastq R1............"<<endl;
				t_fileName=_R1_fname+"_MapNone.fastq";
				WriteFastq(t_fileName, vec_mapNone, vec_mapNone_Q,  mode);
				vec_mapNone_Q.clear();
			  }
			  else
			  {
				  t_fileName=_R1_fname+"_MapNone.fasta";
				  WriteFasta(t_fileName, vec_mapNone,100, mode);
				  
			  }
			  
			  if(_vecSeq_R2.size()>0)
			  {
				  if(_vecSeq_Q_R2.size()>0)
				  {
						t_fileName=_R2_fname+"_MapNone.fastq";
						WriteFastq(t_fileName, vec_mapNone_R2, vec_mapNone_Q_R2,  mode);
						vec_mapNone_Q_R2.clear();
						vec_mapNone_R2.clear();
				  }
				  else
				  {
					  t_fileName=_R2_fname+"_MapNone.fasta";
						WriteFasta(t_fileName, vec_mapNone_R2,100, mode);
						vec_mapNone_R2.clear();
				  }
			  }
		  }
		  numOfWritesDone_unmapped++;
		  vec_mapNone.clear();
      //t_fileName=_unmap_fname;
      //WriteFasta(t_fileName, vec_mapNone,100, mode);
      //#fileCounter_mpBoth<-fileCounter_mpBoth+1;
      //vec_mapNone.clear();
    }
  //cout<<"done with one last write"<<endl;
  
  cout<<endl;
  //cleanup

  delete [] pt_vec_mapped;
  //delete [] pt_vec_mapForward;

  delete [] pt_vec_map_pos_end;

  delete [] pt_vec_map_pos_start;
  delete [] pt_vec_map_pos_end_seq;
  delete [] pt_vec_map_pos_start_seq;
  delete [] pt_vec_map_match_rate;
  delete [] pt_vec_map_overlap;
  //delete [] pt_vec_mapForward_pos_end;

  delete [] numOfWritesDone_mapped;
  //if(tempAS_arr!=NULL)
	//  delete[] tempAS_arr;
  
  //for sequence output and quality string if possible
  delete [] pt_vec_demux;
  if(pt_vec_demux_Q!=NULL)
	  delete[] pt_vec_demux_Q;
  if( pt_vec_demux_R2!=NULL)
	  delete [] pt_vec_demux_R2;
  
  if( pt_vec_demux_Q_R2!=NULL)
	  delete [] pt_vec_demux_Q_R2;
}//end of function of map isotype

