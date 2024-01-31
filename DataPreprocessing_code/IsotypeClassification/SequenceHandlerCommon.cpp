#include "SequenceHandlerCommon.hpp"
#include <string>
using namespace std;


unsigned int CompareStrings(const string& str1, const string& str2)
{
  unsigned int ret=0;
  unsigned int len1, len2, len;
  len1=str1.length();len2=str2.length();
  len=len1;

  
  if(len>len2)
    {
      len=len2;
      ret=len1-len2;
    }
  else
    {
      ret=len2-len1;
    }
  for(unsigned int i=0;i<len;i++)
    {
      
      if(str1.at(i)!=str2.at(i))
	{
	  ret++;
	}
    }
  
  return ret;
}


SequenceString ReverseComplement(SequenceString& seq)
{
  string tempStr=seq.GetSequence();
  SequenceString temp(seq.GetName(), "");
  string tempStrReturn("");
  //cout<<"temStr (seq get sequence):"<<tempStr<<endl;
  //cout<<"length:"<<tempStr.length()<<endl;
  for(unsigned int i=tempStr.length();i>0;i--)
    {
      //cout<<"\tloop i="<<i<<endl;
      switch(tempStr.at(i-1))
	{
	case 'A':
	case 'a':
	  tempStrReturn.push_back('T');
	  break;
	case 'T':
	case 't':
	  tempStrReturn.push_back('A');
	  break;
	case 'U':
	case 'u':
	  tempStrReturn.push_back('A');
	  break;
	case 'G':
	case 'g':
	  tempStrReturn.push_back('C');
	  break;
	case 'C':
	case 'c':
	  tempStrReturn.push_back('G');
	  break;
	case 'Y':
	case 'y':
	  tempStrReturn.push_back('R');
	  break;
	case 'R':
	case 'r':
	  tempStrReturn.push_back('Y');
	  break;
	case 'S':
	case 's':
	  tempStrReturn.push_back('S');
	  break;
	case 'W':
	case 'w':
	  tempStrReturn.push_back('W');
	  break;
	case 'K':
	case 'k':
	  tempStrReturn.push_back('M');
	  break;
	case 'M':
	case 'm':
	  tempStrReturn.push_back('K');
	  break;
	case 'B':
	case 'b':
	  tempStrReturn.push_back('V');
	  break;
	case 'D':
	case 'd':
	  tempStrReturn.push_back('H');
	  break;
	case 'H':
	case 'h':
	  tempStrReturn.push_back('D');
	  break;
	case 'V':
	case 'v':
	  tempStrReturn.push_back('B');
	  break;
	case 'N':
	case 'n':
	default:
	  tempStrReturn.push_back('N');
	  break;
	}
      
    }

  temp.SetSequence(tempStrReturn);

  return temp;
}

SequenceString Reverse(SequenceString& seq)
{
  string tempStr=seq.GetSequence();
  SequenceString temp(seq.GetName(), "");
  //string tempStrReturn(tempStr);
	const char * cstr=tempStr.c_str();
	char * cstr_ret=new char[tempStr.length()+1];
unsigned int index1=0;
unsigned int index2=tempStr.length()-1;
  //cout<<"temStr (seq get sequence):"<<tempStr<<endl;
  //cout<<"length:"<<tempStr.length()<<endl;
  while(index2>=index1)
    {
	//char tempChar=tempStrReturn[index1];
	cstr_ret[index1]=cstr[index2];
	cstr_ret[index2]=cstr[index1];
	index1++;
	index2--;
    }

	cstr_ret[tempStr.length()]='\0';
	string tempStrReturn(cstr_ret);
  temp.SetSequence(tempStrReturn);
	
	delete [] cstr_ret;
  return temp;
}

double MatchBarcodes(const SequenceString& seq, const SequenceString& barcode, const MatchMatrix* mm)
{
  double score=0;
  unsigned seq_len=seq.GetSequence().length();
  string seqStr=seq.GetSequence();
  
  unsigned barcode_len=barcode.GetSequence().length();
  if(seq_len<barcode_len)
    {
      seqStr.append(barcode_len-seq_len,'N');
    }
  string barcodeStr=barcode.GetSequence();
  for(unsigned i=0;i<barcode_len;i++)
    {
      //compare each char for score
      score+=mm->GetScore(seqStr.at(i),barcodeStr.at(i));
    }
  return score;
}
