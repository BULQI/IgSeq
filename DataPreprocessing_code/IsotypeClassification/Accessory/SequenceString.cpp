
#include "SequenceString.hpp"
#include<iostream>
#include <string>
#include <sstream>
#include <stdlib.h>
using namespace std;
SequenceString::SequenceString():c_name(""),c_seq("")
{
  //empty
}
SequenceString::~SequenceString()
{
  //empty
}
SequenceString::SequenceString(const string& _name, const string& _seq) : c_name(_name),c_seq(_seq)
{
  //empty
}
void SequenceString::SetName(const string& _name)
{
  c_name=_name;
}
void SequenceString::SetSequence(const string& _seq)
{
  c_seq=_seq;
  //c_len=_seq.length()
}

const string SequenceString::GetName() const
{
  return c_name;
}

const string SequenceString::GetSequence() const
{
  return c_seq;
}

const unsigned int SequenceString::GetLength() const
{
  return c_seq.length();
}

string SequenceString::toString(bool _fasta) const
{
  ostringstream ss("");
  if(!_fasta)
    {
       ss<<c_seq.length() <<"-character long SequenceString instance\n"<<c_name<<":"<<c_seq<<"\n";
    }
  else
    {
      ss<< ">"<<c_name<<"\n"<<c_seq<<"\n";
    }
  return ss.str();
}

bool SequenceString::operator < (const SequenceString& other) const
{
	int ret=c_seq.compare(other.c_seq);
	//cout<<"compare in the function:"<<ret<<endl;
	if(ret<=0)
		return true;
	else
		return false;
}

void SequenceString::Serialize(ofstream& _ofs) const
{
//do necessary checking
  if(!_ofs.is_open())
    {
      cout<<"**ERROR**: closed file buffer. quit..."<<endl;
      exit(-1);
    }
  //serializing......
  const char* p_char; //pointer used to direct the writing to the file stream
  //save name 
  unsigned l=c_name.length()+1;
  p_char=(char*)&l;
  _ofs.write(p_char, sizeof(unsigned));
  p_char=c_name.c_str();
  _ofs.write(p_char, c_name.length()+1);

  //seq
  l=c_seq.length()+1;
  p_char=(char*)&l;
  _ofs.write(p_char, sizeof(unsigned));
  p_char=c_seq.c_str();
  _ofs.write(p_char, c_seq.length()+1);
}
void SequenceString::Deserialize(ifstream& _ifs)
{
//do necessary checking
  if(!_ifs.is_open())
    {
      cout<<"**ERROR**: closed file buffer. quit..."<<endl;
      exit(-1);
    }
  //serializing......
  char* p_char; //pointer used to direct the writing to the file stream
  //save name
  unsigned l;//=c_name.length()+1;
  p_char=(char*)&l;
  _ifs.read(p_char, sizeof(unsigned));
  p_char=new char[l];
  _ifs.read(p_char,l);
  c_name=string(p_char);
  delete[] p_char;

  //
  //l=c_seq.length()+1;
  p_char=(char*)&l;
  _ifs.read(p_char, sizeof(unsigned));
  p_char=new char[l];
  _ifs.read(p_char, l);
  c_seq=string(p_char);
  delete[] p_char;
}

SequenceString SequenceString::Sub(const unsigned& _start, const unsigned & _end) const
{
  SequenceString ret;

  //get the sub string
  ret.c_name=c_name;
  unsigned temp_start, temp_end;
  temp_start=_start;
  temp_end=_end;
  if(temp_start>c_seq.size())
    {
      temp_start=c_seq.size();
    }
  if(temp_end==(unsigned)-1)
    {
      temp_end=c_seq.size();
    }
  ret.c_seq=c_seq.substr(_start, _end-_start+1);
  return ret;
}

unsigned SequenceString::GetLetterCount(const char& _c) const
{
  //we need to itera through the string to get the count
  unsigned count=0;
  for(unsigned i=0;i<c_seq.size();i++)
    {
      if(c_seq.at(i)==_c)
	{
	  count++;
	}
    }

  return count;
}
