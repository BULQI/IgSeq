#include <stdlib.h>
#include "score.hpp"
#include <cmath>
#include <iostream>
#include <stdio.h>
using namespace std;
//all matrix data are from ncbi ftp://ftp.ncbi.nih.gov/blast/matrices/
//so far supports only BLOSUM50 and nuc4.4

//********************NUC44***********************
//the actually nuc44 matrix
//the scale information is from matlab. not sure where this is coming from, but
//the score are from ncbi.
//order is A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
static int nuc44Int_1d[] ={ 5,  -4,  -4,  -4,  -4,   1,   1,  -4,  -4,   1,  -4,  -1,  -1,  -1,  -2,
			    -4,  5,  -4,  -4,  -4,   1,  -4,   1,   1,  -4,  -1,  -4,  -1,  -1,  -2,
			    -4, -4,   5,  -4,   1,  -4,   1,  -4,   1,  -4,  -1,  -1,  -4,  -1,  -2,
			    -4, -4,  -4,   5,   1,  -4,  -4,   1,  -4,   1,  -1,  -1,  -1,  -4,  -2,
			    -4, -4,   1,   1,  -1,  -4,  -2,  -2,  -2,  -2,  -1,  -1,  -3,  -3,  -1,
			    1,   1,  -4,  -4,  -4,  -1,  -2,  -2,  -2,  -2,  -3,  -3,  -1,  -1,  -1,
			    1,  -4,   1,  -4,  -2,  -2,  -1,  -4,  -2,  -2,  -3,  -1,  -3,  -1,  -1,
			    -4,   1,  -4,   1,  -2,  -2,  -4,  -1,  -2,  -2,  -1,  -3,  -1,  -3,  -1,
			    -4,   1,   1,  -4,  -2,  -2,  -2,  -2,  -1,  -4,  -1,  -3,  -3,  -1,  -1,
			    1,  -4,  -4,   1,  -2,  -2,  -2,  -2,  -4,  -1,  -3,  -1,  -1,  -3,  -1,
			    -4,  -1,  -1,  -1,  -1,  -3,  -3,  -1,  -1,  -3,  -1,  -2,  -2,  -2,  -1,
			    -1,  -4,  -1,  -1,  -1,  -3,  -1,  -3,  -3,  -1,  -2,  -1,  -2,  -2,  -1,
			    -1,  -1,  -4,  -1,  -3,  -1,  -3,  -1,  -3,  -1,  -2,  -2,  -1,  -2,  -1 , 
			    -1,  -1,  -1,  -4,  -3,  -1,  -1,  -3,  -1,  -3,  -2,  -2,  -2,  -1,  -1,
			    -2,  -2,  -2,  -2,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1
		            };
//*****the following is used to turn a two D array into one data, so we can pass it to the function.
//nuc44Int[] is an array of pointers, pointing to the original nuc44Int_1d. with this, we can go found the scores
//in a one data format, pointing to each row of the above array. All other score matrices are implemented
//like this.
static const int* nuc44Int[] ={nuc44Int_1d,nuc44Int_1d+15,nuc44Int_1d+30,nuc44Int_1d+45,nuc44Int_1d+60,nuc44Int_1d+75,
			 nuc44Int_1d+90,nuc44Int_1d+105, nuc44Int_1d+120, nuc44Int_1d+135,nuc44Int_1d+150, 
			 nuc44Int_1d+165, nuc44Int_1d+15*12,nuc44Int_1d+15*13, nuc44Int_1d+15*14 
                };
const char ScoreMatrix::NucAlphabet[]={ 'A' ,  'T' ,  'G'  ,'C',   'S',   'W',   'R' ,  'Y',   'K',   'M',   'B' ,  'V' ,  'H' , 'D',   'N'};
//the reason we have to do it this way is because we can not do it simply by passing 
//a two D array to the function. it is not supported. you have to know the size of
//the dimension size except the first one.
//of course, another alternative is to define template and then do it.
ScoreMatrix nuc44(nuc44Int,0.277316, ScoreMatrix::NucAlphabet,15);//again not sure where this scale come from by Matlab.
//****************************************

//***************************************
/*nuc44_variant matrix, set up this in order for the D alignment according to Matlab code. make the highly similar sequence aligned
 *set up the mismatch to be higher cost. -14 everywhere.
 */
//order is A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
static int nuc44Int_1d_v[] ={ 5,  -14,  -14,  -14,  -14,   1,   1,  -14,  -14,   1,  -14,  -5,  -5,  -5,  -8,
			    -14,  5,  -14,  -14,  -14,   1,  -14,   1,   1,  -14,  -5,  -14,  -5,  -5,  -8,
			    -14, -14,   5,  -14,   1,  -14,   1,  -14,   1,  -14,  -5,  -5,  -14,  -5,  -8,
			    -14, -14,  -14,   5,   1,  -14,  -14,   1,  -14,   1,  -5,  -5,  -5,  -14,  -8,
			    -14, -14,   1,   1,  -5,  -14,  -8,  -8,  -8,  -8,  -5,  -5,  -12,  -12,  -5,
			    1,   1,  -14,  -14,  -14,  -5,  -8,  -8,  -8,  -8,  -12,  -12,  -5,  -5,  -5,
			    1,  -14,   1,  -14,  -8,  -8,  -5,  -14,  -8,  -8,  -12,  -5,  -12,  -5,  -5,
			    -14,   1,  -14,   1,  -8,  -8,  -14,  -5,  -8,  -8,  -5,  -12,  -5,  -12,  -5,
			    -14,   1,   1,  -14,  -8,  -8,  -8,  -8,  -5,  -14,  -5,  -12,  -12,  -5,  -5,
			    1,  -14,  -14,   1,  -8,  -8,  -8,  -8,  -14,  -5,  -12,  -5,  -5,  -12,  -5,
			    -14,  -5,  -5,  -5,  -5,  -12,  -12,  -5,  -5,  -12,  -5,  -8,  -8,  -8,  -5,
			    -5,  -14,  -5,  -5,  -5,  -12,  -5,  -12,  -12,  -5,  -8,  -5,  -8,  -8,  -5,
			    -5,  -5,  -14,  -5,  -12,  -5,  -12,  -5,  -12,  -5,  -8,  -8,  -5,  -8,  -1 , 
			    -5,  -5,  -5,  -14,  -12,  -5,  -5,  -12,  -5,  -12,  -8,  -8,  -8,  -5,  -5,
			    -8,  -8,  -8,  -8,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5,  -5
		            };
static const int* nuc44Int_v[] ={nuc44Int_1d_v,nuc44Int_1d_v+15,nuc44Int_1d_v+30,nuc44Int_1d_v+45,nuc44Int_1d_v+60,nuc44Int_1d_v+75,
			 nuc44Int_1d_v+90,nuc44Int_1d_v+105, nuc44Int_1d_v+120, nuc44Int_1d_v+135,nuc44Int_1d_v+150, 
			 nuc44Int_1d_v+165, nuc44Int_1d_v+15*12,nuc44Int_1d_v+15*13, nuc44Int_1d_v+15*14 
                };
ScoreMatrix nuc44_v(nuc44Int_v,0.277316, ScoreMatrix::NucAlphabet,15);//again not sure where this scale come from by Matlab.
//****************************************

//*********************************************
//BLOSUM50, aa
//the scale information is from ncbi. not sure where this is coming from, but
//the score are from ncbi.
//order is A,R,N,D,C,Q,E,G,H,I,L,K,M,F,P,S,T,W,Y,V,B,Z,X,*
static int blosum50_1d[] ={  5,-2,-1,-2,-1,-1,-1, 0,-2,-1,-2,-1,-1,-3,-1, 1, 0,-3,-2, 0,-2,-1,-1,-5,
			    -2, 7,-1,-2,-4, 1, 0,-3, 0,-4,-3, 3,-2,-3,-3,-1,-1,-3,-1,-3,-1, 0,-1,-5,
			    -1,-1, 7, 2,-2, 0, 0, 0, 1,-3,-4, 0,-2,-4,-2, 1, 0,-4,-2,-3, 4, 0,-1,-5,
			    -2,-2, 2, 8,-4, 0, 2,-1,-1,-4,-4,-1,-4,-5,-1, 0,-1,-5,-3,-4, 5, 1,-1,-5,
			    -1,-4,-2,-4,13,-3,-3,-3,-3,-2,-2,-3,-2,-2,-4,-1,-1,-5,-3,-1,-3,-3,-2,-5,
			     -1,1, 0, 0,-3, 7, 2,-2, 1,-3,-2, 2, 0,-4,-1, 0,-1,-1,-1,-3, 0, 4,-1,-5,
			     -1,0, 0, 2,-3, 2, 6,-3, 0,-4,-3, 1,-2,-3,-1,-1,-1,-3,-2,-3, 1, 5,-1,-5,
			      0,-3,0,-1,-3,-2,-3,8,-2,-4,-4,-2,-3,-4,-2,0,-2,-3,-3,-4,-1,-2,-2,-5,
			     -2,0,1,-1,-3,1,0,-2,10,-4,-3,0,-1,-1,-2,-1,-2,-3,2,-4,0,0,-1,-5,
			     -1,-4,-3,-4,-2,-3,-4,-4,-4,5,2,-3,2,0,-3,-3,-1,-3,-1,4,-4,-3,-1,-5,
			     -2,-3,-4,-4,-2,-2,-3,-4,-3,2,5,-3,3,1,-4,-3,-1,-2,-1,1,-4,-3,-1,-5,
			     -1,3,0,-1,-3,2,1,-2,0,-3,-3,6,-2,-4,-1,0,-1,-3,-2,-3,0,1,-1,-5,
			     -1,-2,-2,-4,-2,0,-2,-3,-1,2,3,-2,7,0,-3,-2,-1,-1,0,1,-3,-1,-1,-5,
			     -3,-3,-4,-5,-2,-4,-3,-4,-1,0,1,-4,0,8,-4,-3,-2,1,4,-1,-4,-4,-2,-5,
			     -1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10,-1,-1,-4,-3,-3,-2,-1,-2,-5,
			     1,-1,1,0,-1,0,-1,0,-1,-3,-3,0,-2,-3,-1,5,2,-4,-2,-2,0,0,-1,-5,
			     0,-1,0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1,2,5,-3,-2,0,0,-1,0,-5,
			     -3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1,1,-4,-4,-3,15,2,-3,-5,-2,-3,-5,
			     -2,-1,-2,-3,-3,-1,-2,-3,2,-1,-1,-2,0,4,-3,-2,-2,2,8,-1,-3,-2,-1,-5,
			     0,-3,-3,-4,-1,-3,-3,-4,-4,4,1,-3,1,-1,-3,-2,0,-3,-1,5,-4,-3,-1,-5,
			     -2,-1,4,5,-3,0,1,-1,0,-4,-4,0,-3,-4,-2,0,0,-5,-3,-4,5,2,-1,-5,
			     -1,0,0,1,-3,4,5,-2,0,-3,-3,1,-1,-4,-1,0,-1,-2,-2,-3,2,5,-1,-5,
			     -1,-1,-1,-1,-2,-1,-1,-2,-1,-1,-1,-1,-1,-2,-2,-1,0,-3,-1,-1,-1,-1,-1,-5,
			     -5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,1 };
static const int* blosum50Int[] ={blosum50_1d,blosum50_1d+24,blosum50_1d+24*2,blosum50_1d+24*3, blosum50_1d+24*4,
			    blosum50_1d+24*5,blosum50_1d+24*6,blosum50_1d+24*7,blosum50_1d+24*8,blosum50_1d+24*9,
			    blosum50_1d+24*10,blosum50_1d+24*11,blosum50_1d+24*12,blosum50_1d+24*13,blosum50_1d+24*14,
			    blosum50_1d+24*15,blosum50_1d+24*16,blosum50_1d+24*17,blosum50_1d+24*18,blosum50_1d+24*19,
			    blosum50_1d+24*20,blosum50_1d+24*21,blosum50_1d+24*22,blosum50_1d+24*23};



ScoreMatrix blosum50(blosum50Int,log(2)/3,ScoreMatrix::AaAlphabet,24);//again not sure where this scale come from, but I got it by google.
                                           //it is from NCBI C toolkit cross reference
                                           //http://www.ncbi.nlm.nih.gov/IEB/ToolBox/C_DOC/lxr/source/data/BLOSUM50
const char ScoreMatrix::AaAlphabet[]={'A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X','*'};
//****************************************************************


//constructor
ScoreMatrix::ScoreMatrix(const int*  matrix[ ] , const double& scale, const char* alphabet, const int& alphabet_array_size)
  :
  _matrix(matrix),_scale(scale),
  _alphabet(alphabet),_scaledScoreFlag(false), c_alphabet_array_size(alphabet_array_size)
{
  //do nothing 
}
const int ScoreMatrix::GetAlphabetLength()
{
  return this->c_alphabet_array_size;
}
const int** ScoreMatrix::GetScoreMatrix() const
{
  return this->_matrix;
}
  
double ScoreMatrix::GetScale() const
{
  return this->_scale;
}

ScoreMatrix::~ScoreMatrix()
{
  //here, we simply set the pointer to null, because we did not allocate the memory when we build it.
  //the memory holding the data is declared outside and we simply point to it when we build the object.
  this->_matrix=NULL;
}
const char* ScoreMatrix::GetAlphabet() const
{
  return this->_alphabet;
}
  
double ScoreMatrix::GetScore(const char& first, const char& second) const
{
  //look through
  
  int firstIndex=-1, secondIndex=-1;
  char firstUp=toupper(first);
  char secondUp=toupper(second);
  for(int i=0;i<c_alphabet_array_size;i++)
    {
      
      if(firstUp==_alphabet[i])
	{
	  firstIndex=i;
	}
      if(secondUp==_alphabet[i])
	{
	  secondIndex=i;
	}
      
      if(firstIndex!=-1&&secondIndex!=-1)
	{
	  break;
	}
    }
 
  if(firstIndex==-1)
	{
	  cout<<"********WARNING:can not find the letter \'"<<first<<"\', using \'"<<_alphabet[c_alphabet_array_size-1]<<"\' instead."<<endl; 
	  firstIndex=c_alphabet_array_size-1;
	}
  if(secondIndex==-1)
	{
	  cout<<"********WARNING:can not find the letter \'"<<first<<"\', using \'"<<_alphabet[c_alphabet_array_size-1]<<"\' instead."<<endl; 
	  secondIndex=c_alphabet_array_size-1;
	}

  if(_scaledScoreFlag)
    {
      //cout<<"in no scaled case:**"<<endl;
      return _matrix[firstIndex][secondIndex]/_scale;
    }
  else
    {
      //cout<<"in yes scaled case:"<<endl;
      return _matrix[firstIndex][secondIndex];
    }
}

void ScoreMatrix::SetScaledScoreFlag(const bool& flag)
{
  this->_scaledScoreFlag=flag;
}

//***************************************************************
//the tsm1 matrix
//order is A   T   C H K
static int tsm1Int_1d[] ={  10,  -4,  -4,  -4,  -4,
			     -4,   5,  -4,  -4,   5,
			     -4,  -4,   5,  -4,   6,
			     -4,  -4,  -4,  -4,  -4,
			     -4,   5,   6,  -4,   5,           };
static const int* tsm1Int[] ={tsm1Int_1d,tsm1Int_1d+5*1,tsm1Int_1d+5*2,tsm1Int_1d+5*3,tsm1Int_1d+5*4                };
static const char tsm1Alphabet[]={ 'A' , 'C',   'H',    'K', 'T'   };
//the reason we have to do it this way is because we can not do it simply by passing 
//a two D array to the function. it is not supported. you have to know the size of
//the dimension size except the first one.
//of course, another alternative is to define template and then do it.
ScoreMatrix tsm1(tsm1Int,1, tsm1Alphabet,5);//again not sure where this scale come from by Matlab.
//***************

//***************************************************
//the tsm2 matrix
//order is A   T   C  G
static int tsm2Int_1d[] ={  10,  -9,  -9,  -9,
			    -9,  10,   -9,   -9,
			    -9,  -9,   10,   -9,
			    -9,  -9,   -9,   10
};
static const int* tsm2Int[] ={tsm2Int_1d,tsm2Int_1d+4*1,tsm2Int_1d+4*2,tsm2Int_1d+4*3                };
static const char tsm2Alphabet[]={ 'A' , 'T', 'C',   'G'   };
//the reason we have to do it this way is because we can not do it simply by passing 
//a two D array to the function. it is not supported. you have to know the size of
//the dimension size except the first one.
//of course, another alternative is to define template and then do it.
ScoreMatrix tsm2(tsm2Int,1, tsm2Alphabet,4);
//this matrix is used to test the example in the reference:Barton,G. 1993. CABIO. An efficient Algorithm to locate all locally ptimal alignments between two sequences alowing for gaps. Vol 9. no.6. 1993 P729-34
//******************************************

//*************************NUC44 high penalty of MisMatch
//the actually nuc44 matrix modified with high penalty of mismatch
//the scale information is kept same as the original one from NUC44
//the score are from ncbi.
//order is A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
static int factor=4;
static int nuc44HPInt_1d[] ={ 5,  -4*factor,-4*factor,  -4*factor,  -4*factor, 1,   1,  -4*factor,  -4*factor,   1,    -4*factor,  -1*factor,  -1*factor,  -1*factor,  -2*factor,
			    -4*factor,  5,  -4*factor,  -4*factor,  -4*factor,   1,  -4*factor,   1,   1,  -4*factor,  -1*factor,  -4*factor,  -1*factor,  -1*factor,  -2*factor,
			    -4*factor, -4*factor,   5,  -4*factor,   1,  -4*factor,   1,  -4*factor,   1,  -4*factor,  -1*factor,  -1*factor,  -4*factor,  -1*factor,  -2*factor,
			    -4*factor, -4*factor,  -4*factor,   5,   1,  -4*factor,  -4*factor,   1,  -4*factor,   1,  -1*factor,  -1*factor,  -1*factor,  -4*factor,  -2*factor,
			    -4*factor, -4*factor,   1,   1,  -1*factor,  -4*factor,  -2*factor,  -2*factor, -2*factor,  -2*factor,  -1*factor,  -1*factor,  -3*factor,  -3*factor,  -1*factor,
			    1,   1,  -4*factor,  -4*factor,  -4*factor,  -1*factor,  -2*factor,  -2*factor,  -2*factor,  -2*factor,  -3*factor,  -3*factor,  -1*factor,  -1*factor,  -1*factor,
			    1,  -4*factor,   1,  -4*factor,  -2*factor,  -2*factor,  -1*factor,  -4*factor,  -2*factor,  -2*factor,  -3*factor,  -1*factor,  -3*factor,  -1*factor,  -1*factor,
			    -4*factor,   1,  -4*factor,   1,  -2*factor,  -2*factor,  -4*factor,  -1*factor,  -2*factor,  -2*factor,  -1*factor,  -3*factor,  -1*factor,  -3*factor,  -1*factor,
			    -4*factor,   1,   1,  -4*factor,  -2*factor,  -2*factor,  -2*factor,  -2*factor,  -1*factor,  -4*factor,  -1*factor,  -3*factor,  -3*factor,  -1*factor,  -1*factor,
			    1,  -4*factor,  -4*factor,   1,  -2*factor,  -2*factor,  -2*factor,  -2*factor,  -4*factor,  -1*factor,  -3*factor,  -1*factor,  -1*factor,  -3*factor,  -1*factor,
			    -4*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -3*factor,  -3*factor,  -1*factor,  -1*factor,  -3*factor,  -1*factor,  -2*factor,  -2*factor,  -2*factor,  -1*factor,
			    -1*factor,  -4*factor,  -1*factor,  -1*factor,  -1*factor,  -3*factor,  -1*factor,  -3*factor,  -3*factor,  -1*factor,  -2*factor,  -1*factor,  -2*factor,  -2*factor,  -1*factor,
			    -1*factor,  -1*factor,  -4*factor,  -1*factor,  -3*factor,  -1*factor,  -3*factor,  -1*factor,  -3*factor,  -1*factor,  -2*factor,  -2*factor,  -1*factor,  -2*factor,  -1*factor , 
			    -1*factor,  -1*factor,  -1*factor,  -4*factor,  -3*factor,  -1*factor,  -1*factor,  -3*factor,  -1*factor,  -3*factor,  -2*factor,  -2*factor,  -2*factor,  -1*factor,  -1*factor,
			    -2*factor,  -2*factor,  -2*factor,  -2*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor,  -1*factor
		            };
static const int* nuc44HPInt[] ={nuc44HPInt_1d,nuc44HPInt_1d+15,nuc44HPInt_1d+30,nuc44HPInt_1d+45,nuc44HPInt_1d+60,nuc44HPInt_1d+75,
			 nuc44HPInt_1d+90,nuc44HPInt_1d+105, nuc44HPInt_1d+120, nuc44HPInt_1d+135,nuc44HPInt_1d+150, 
			 nuc44HPInt_1d+165, nuc44HPInt_1d+15*12,nuc44HPInt_1d+15*13, nuc44HPInt_1d+15*14 
                };
const char NucHPAlphabet[]={ 'A' ,  'T' ,  'G'  ,'C',   'S',   'W',   'R' ,  'Y',   'K',   'M',   'B' ,  'V' ,  'H' , 'D',   'N'};
//the reason we have to do it this way is because we can not do it simply by passing 
//a two D array to the function. it is not supported. you have to know the size of
//the dimension size except the first one.
//of course, another alternative is to define template and then do it.
ScoreMatrix nuc44HP(nuc44HPInt,0.277316, NucHPAlphabet,15);//again not sure where this scale come from by Matlab.

//*************************NUC44DM1, nuc44 one way degenerate match.
//the actually nuc44 matrix modified in order to allow the degenerated match between sequences and isotypes/adaptors
//the special thing about this score matrix is that it is directional, meaning it is fine to have degeneracy on
//the isotype/adaptors, but not so on the sequences. For example, if there is a N in sequence, we treat it
//as a mismatch, no matter what is in the isotype/adaptor. But if there is degenerated nt on the adaptor/isotye,
//it could a match depending on what is in the sequence. score ('N', 't') is mismatch, but score ('t','N') is
//always good match, assuming first nt is from the sequence and second from the adaptor.
//we manually curate the matrix based on the nuc44 matrix.
//.
//order is A   T   G   C   S   W   R   Y   K   M   B   V   H   D   N
//Again, degeneracy on the isotype/adapter (second) is allowed and could be a good match;
//but degeneracy on the sequence (first input) is normally a bad thing, meaning
//no certainty for the reading. It can only be a match when we coincidently the
//same degeneracy on the isotype/adapter.
//we always assume the sequence is the first input and isotype/adapter the second.
//the rows are for score for each same nt of the sequence against different nts from isotype/adapter.
//each column is for the same nt from the isotype/adapter against the sequence (first)
                           /*  A    T    G    C    S    W    R    Y    K    M    B   V   H   D   N*/
static int nuc44DM1Int_1d[] ={ 5,  -4,  -4,  -4,  -4,   5,   5,  -4,  -4,   5,  -4,  5,  5,  5,  5,
	/*T*/		      -4,   5,  -4,  -4,  -4,   5,  -4,   5,   5,  -4,   5, -4,  5,  5,  5,
	/*G*/	              -4,  -4,   5,  -4,   5,  -4,   5,  -4,   5,  -4,  -5,  5, -4,  5,  5,
	/*C*/  	              -4,  -4,  -4,   5,   5,  -4,  -4,   5,  -4,   5,   5,  5,  5, -4,  5,
	/*S*/		      -4,  -4,   1,   1,   5,  -4,   1,   1,   1,   1,   5,  5,  1,  1,  5,
	/*W*/		       1,   1,  -4,  -4,  -4,   5,   1,   1,   1,   1,   1,  1,  5,  5,  5,
	/*R*/		       1,  -4,   1,  -4,   1,   1,   5,  -4,   1,   1,   1,  5,  1,  5,  5,
	/*Y*/		      -4,   1,  -4,   1,   1,   1,  -4,   5,   1,   1,   5,  1,  5,  1,  5,
	/*K*/		      -4,   1,   1,  -4,   1,   1,   1,   1,   5,  -4,   5,  1,  1,  5,  5,
	/*M*/		       1,  -4,  -4,   1,   1,   1,   1,   1,  -4,   5,   1,  5,  5,  1,  5,
	/*B*/		      -4,  -1,  -1,  -1,   2,  -1,  -1,   2,   2,  -1,   5,  2,  2,  2,  5,
	/*V*/		      -1,  -4,  -1,  -1,   2,  -1,   2,  -1,  -1,   2,   2,  5,  2,  2,  5,
	/*H*/		      -1,  -1,  -4,  -1,  -1,   2,  -1,   2,  -1,   2,   2,  2,  5,  2,  5,
	/*D*/		      -1,  -1,  -1,  -4,  -1,   2,   2,  -1,   2,  -1,   2,  2,  2,  5,  5,
	/*N*/		      -2,  -2,  -2,  -2,   1,   1,   1,   1,   1,   1,   3 , 3,  3,  3,  5
		             };
/* how to calculate the score for degnerated match, row first and column second.
 *eg. score('T','S'), since S -> G or C, so, this is mismatch anyway
 *  score('T', 'W'), since W->A or T, so this is a match. 
 *  score ('S', 'W'), since s->G or C and w->A or T, so this is 
 *     mismatch.
 *  score('S','R'), since s->G or C and R->A or G, when s is G, we have
 *     score('G', 'R')=5 (match); when s is C, score('C', 'R')=-4. If we
 *     assume it is equal likely for S to be G and C, we have 
 *     score=0.5*5+0.5*(-4)=0.5, we round it up to 1 for the integer format  
 */
static const int* nuc44DM1Int[] ={nuc44DM1Int_1d,nuc44DM1Int_1d+15,nuc44DM1Int_1d+30,nuc44DM1Int_1d+45,nuc44DM1Int_1d+60,nuc44DM1Int_1d+75,
			 nuc44DM1Int_1d+90,nuc44DM1Int_1d+105, nuc44DM1Int_1d+120, nuc44DM1Int_1d+135,nuc44DM1Int_1d+150, 
			 nuc44DM1Int_1d+165, nuc44DM1Int_1d+15*12,nuc44DM1Int_1d+15*13, nuc44DM1Int_1d+15*14 
                };
const char NucDM1Alphabet[]={ 'A' ,  'T' ,  'G'  ,'C',   'S',   'W',   'R' ,  'Y',   'K',   'M',   'B' ,  'V' ,  'H' , 'D',   'N'};
//IUPAC code: S -> G or C;   W -> A or T; R -> A or G; Y -> C or T; K -> G or T;
//   M -> A or C ; B -> C, G or T; V -> A or C or G; H -> A, C or T; D->A, G or T;
//N-> any base
//the reason we have to do it this way is because we can not do it simply by passing 
//a two D array to the function. it is not supported. you have to know the size of
//the dimension size except the first one.
//of course, another alternative is to define template and then do it.
ScoreMatrix nuc44DM1(nuc44DM1Int,0.277316, NucDM1Alphabet,15);//again not sure where this scale come from by Matlab.

/*-------------------------MatchMatrix----------------------------------------------*/

//constructor
MatchMatrix::MatchMatrix(const double*  matrix[ ] ,  const char* alphabet, const int& alphabet_array_size)
  :
  _matrix(matrix),
  _alphabet(alphabet), c_alphabet_array_size(alphabet_array_size)
{
  //do nothing 
}
const int MatchMatrix::GetAlphabetLength()
{
  return this->c_alphabet_array_size;
}
const double** MatchMatrix::GetMatchMatrix() const
{
  return this->_matrix;
}
  
MatchMatrix::~MatchMatrix()
{
  //here, we simply set the pointer to null, because we did not allocate the memory when we build it.
  //the memory holding the data is declared outside and we simply point to it when we build the object.
  this->_matrix=NULL;
}
const char* MatchMatrix::GetAlphabet() const
{
  return this->_alphabet;
}
//here, need to be care, we assume orders and asymmetry degeneracy
double MatchMatrix::GetScore(const char& first, const char& second) const
{
  //look through
  
  int firstIndex=-1, secondIndex=-1;
  char firstUp=toupper(first);
  char secondUp=toupper(second);
  for(int i=0;i<c_alphabet_array_size;i++)
    {
      
      if(firstUp==_alphabet[i])
	{
	  firstIndex=i;
	}
      if(secondUp==_alphabet[i])
	{
	  secondIndex=i;
	}
      
      if(firstIndex!=-1&&secondIndex!=-1)
	{
	  break;
	}
    }
 
  if(firstIndex==-1)
	{
	  cout<<"********WARNING:can not find the letter \'"<<first<<"\', using \'"<<_alphabet[c_alphabet_array_size-1]<<"\' instead."<<endl; 
	  firstIndex=c_alphabet_array_size-1;
	}
  if(secondIndex==-1)
	{
	  cout<<"********WARNING:can not find the letter \'"<<first<<"\', using \'"<<_alphabet[c_alphabet_array_size-1]<<"\' instead."<<endl; 
	  secondIndex=c_alphabet_array_size-1;
	}


  //return _matrix[firstIndex][secondIndex]/_scale;
  return _matrix[firstIndex][secondIndex];
}



//now define a match matrix mmf, match matrix feng
//strict score matrix, NOT USED YET. need to work on it more.
//the idea is that we strictly disallow the degeneracy. <--- implement this later!!!
               /*      A      T      G      C       S       W       R       Y       K      M      B     V     H      D    N*/
static double mmf_strict_1d[] ={ 
			   /*A*/   0,     1,     1,     1,      1,      0,      0,      1,      1,     0,     1,    0,    0,     0,   0,
       		   /*T*/   1,     0,     1,     1,      1,      0,      1,      0,      0,     1,     0,    1,    0,     0,   0,
		       /*G*/   1,     1,     0,     1,      0,      1,      0,      1,      0,     1,     0,    0,    1,     0,   0,
		       /*C*/   1,     1,     1,     0,      0,      1,      1,      0,      1,     0,     0,    0,    0,     1,   0,
		       /*S*/   1,     1,     0.5,   0.5,    0,      1,      0.5,    0.5,    0.5,   0.5,   0,    0,    0.5,   0.5, 0,
		       /*W*/   0.5,   0.5,   1,     1,      1,      0,      0.5,    0.5,    0.5,   0.5,   0.5,  0.5,  0,     0,   0,
		       /*R*/   0.5,   1,     0.5,   1,      0.5,    0.5,    0,      1,      0.5,   0.5,   0.5,  0,    0.5,   0,   0,
		       /*Y*/   1,     0.5,   1,     0.5,    0.5,    0.5,    1,      0,      0.5,   0.5,   0,    0.5,  0,     0.5, 0,
		       /*K*/   1,     0.5,   0.5,   1,      0.5,    0.5,    0.5,    0.5,    0,     1,     0,    0.5,  0.5,   0,   0,
		       /*M*/   0.5,   1,     1,     0.5,    0.5,    0.5,    0.5,    0.5 ,   1,     0,     0.5,  0,    0,     0.5, 0,
		       /*B*/   1,     0.666, 0.666, 0.666,  0.333,  0.666,  0.666,  0.333,  0.333, 0.666, 0,    0.333,0.333, 0.333,0,
		       /*V*/   0.666, 1,     0.666, 0.666,  0.333,  0.666,  0.333,  0.666,  0.666, 0.333, 0.333,0,    0.333, 0.333,0,
		       /*H*/   0.666, 0.666, 1,     0.666,  0.666,  0.333,  0.666,  0.333 , 0.666, 0.333, 0.333,0.333,0,     0.333,0,
		       /*D*/   0.666, 0.666, 0.666, 1,      0.666,  0.333,  0.333,  0.666,  0.333, 0.666, 0.333,0.333,0.333, 0,    0,
		       /*N*/   0.75,  0.75,  0.75,  0.75,   0.5,    0.5,    0.5,    0.5,    0.5,   0.5,   0.25, 0.25, 0.25,  0.25, 0
		             };
static const double* mmf_strict_score[] ={mmf_strict_1d,mmf_strict_1d+15,mmf_strict_1d+30,mmf_strict_1d+45,mmf_strict_1d+60,mmf_strict_1d+75,
			 mmf_strict_1d+90,mmf_strict_1d+105, mmf_strict_1d+120, mmf_strict_1d+135,mmf_strict_1d+150, 
			 mmf_strict_1d+165, mmf_strict_1d+15*12,mmf_strict_1d+15*13, mmf_strict_1d+15*14 
                };

               /*      A      T      G      C       S       W       R       Y       K      M      B     V     H      D    N*/
static double mmf_1d[] ={ 
			   /*A*/   0,     1,     1,     1,      1,      0,      0,      1,      1,     0,     1,    0,    0,     0,   0,
       		   /*T*/   1,     0,     1,     1,      1,      0,      1,      0,      0,     1,     0,    1,    0,     0,   0,
		       /*G*/   1,     1,     0,     1,      0,      1,      0,      1,      0,     1,     0,    0,    1,     0,   0,
		       /*C*/   1,     1,     1,     0,      0,      1,      1,      0,      1,     0,     0,    0,    0,     1,   0,
		       /*S*/   1,     1,     0.5,   0.5,    0,      1,      0.5,    0.5,    0.5,   0.5,   0,    0,    0.5,   0.5, 0,
		       /*W*/   0.5,   0.5,   1,     1,      1,      0,      0.5,    0.5,    0.5,   0.5,   0.5,  0.5,  0,     0,   0,
		       /*R*/   0.5,   1,     0.5,   1,      0.5,    0.5,    0,      1,      0.5,   0.5,   0.5,  0,    0.5,   0,   0,
		       /*Y*/   1,     0.5,   1,     0.5,    0.5,    0.5,    1,      0,      0.5,   0.5,   0,    0.5,  0,     0.5, 0,
		       /*K*/   1,     0.5,   0.5,   1,      0.5,    0.5,    0.5,    0.5,    0,     1,     0,    0.5,  0.5,   0,   0,
		       /*M*/   0.5,   1,     1,     0.5,    0.5,    0.5,    0.5,    0.5 ,   1,     0,     0.5,  0,    0,     0.5, 0,
		       /*B*/   1,     0.666, 0.666, 0.666,  0.333,  0.666,  0.666,  0.333,  0.333, 0.666, 0,    0.333,0.333, 0.333,0,
		       /*V*/   0.666, 1,     0.666, 0.666,  0.333,  0.666,  0.333,  0.666,  0.666, 0.333, 0.333,0,    0.333, 0.333,0,
		       /*H*/   0.666, 0.666, 1,     0.666,  0.666,  0.333,  0.666,  0.333 , 0.666, 0.333, 0.333,0.333,0,     0.333,0,
		       /*D*/   0.666, 0.666, 0.666, 1,      0.666,  0.333,  0.333,  0.666,  0.333, 0.666, 0.333,0.333,0.333, 0,    0,
		       /*N*/   0.75,  0.75,  0.75,  0.75,   0.5,    0.5,    0.5,    0.5,    0.5,   0.5,   0.25, 0.25, 0.25,  0.25, 0
		             };
//IUPAC code: S -> G or C;   W -> A or T; R -> A or G; Y -> C or T; K -> G or T;
//   M -> A or C ; B -> C, G or T; V -> A or C or G; H -> A, C or T; D->A, G or T;

/* how to calculate the score for degnerated match, 
 *eg. score('T','S'), since S -> G or C, so, this is mismatch anyway, 1
 *  score('T', 'W'), since W->A or T, so this is a match,0. 
 *  score ('S', 'W'), since s->G or C and w->A or T, so this is 
 *     mismatch, score of 1.
 *  score('S','R'), since s->G or C and R->A or G, when s is G, we have
 *     score('G', 'R')=0 (match); when s is C, score('C', 'R')=1, If we
 *     assume it is equal likely for S to be G and C, we have 
 *     score=0*0.5+0.5*(1)=0.5, we have 0.5. so all between 0 and 1.
 *in the end, we add all scores to get a match score of the sequence.
 *that is basically how many mismatches for the matching.  
 */
static const double* mmf_score[] ={mmf_1d,mmf_1d+15,mmf_1d+30,mmf_1d+45,mmf_1d+60,mmf_1d+75,
			 mmf_1d+90,mmf_1d+105, mmf_1d+120, mmf_1d+135,mmf_1d+150, 
			 mmf_1d+165, mmf_1d+15*12,mmf_1d+15*13, mmf_1d+15*14 
                };
const char mmfAlphabet[]={ 'A' ,  'T' ,  'G'  ,'C',   'S',   'W',   'R' ,  'Y',   'K',   'M',   'B' ,  'V' ,  'H' , 'D',   'N'};
//IUPAC code: S -> G or C;   W -> A or T; R -> A or G; Y -> C or T; K -> G or T;
//   M -> A or C ; B -> C, G or T; V -> A or C or G; H -> A, C or T; D->A, G or T;
//N-> any base
//the reason we have to do it this way is because we can not do it simply by passing 
//a two D array to the function. it is not supported. you have to know the size of
//the dimension size except the first one.
//of course, another alternative is to define template and then do it.
MatchMatrix mmf(mmf_score, mmfAlphabet,15);//again not sure where this scale come from by Matlab.
MatchMatrix mmf_strict(mmf_strict_score, mmfAlphabet, 15);