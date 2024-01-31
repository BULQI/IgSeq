#include "FASTQ.hpp"
#include <sstream>

Fastq::Fastq(const string& _seqName, const SequenceString& _ss, const string& _q ):
	c_name(_seqName), c_ss(_ss), c_q(_q)
	{
		//empty
	}
Fastq::~Fastq()//destructor empty one
{
	//empty 
}


void Fastq::SetName(const string& _name)
{
	c_name.assign(_name);
}
void Fastq::SetSequenceString(const SequenceString& _seq)
  {
	  c_ss=_seq;
  }
void Fastq::SetQualityString(const string& _q)
{
	c_q.assign(_q);
}
  
  const string Fastq::GetName()const
  {
	  return c_name;
  }
  const string Fastq::GetSequence()const
  {
	  return c_ss.GetSequence();
  }
  const SequenceString Fastq::GetSequenceString()const
  {
	  return c_ss;
  }
  const string Fastq::GetQualityString() const
  {
	  return c_q;
  }

  string Fastq::toString(bool _fasta) const
  {
	  ostringstream oss("");
	  if(!_fasta)
		{
			oss<<"@"<<c_name<<"\n";
			oss<<c_ss.GetSequence()<<"\n";
			oss<<"+\n";
			oss<<c_q<<"\n";
		   //oss<<c_ss.toString(_fasta)<<"\n";
		   //oss<<
		}
	  else
		{
		  oss<< ">"<<c_name<<"\n"<<c_ss.GetSequence()<<"\n";//<<"+\n"<<c_q<<"\n";
		  //oss<<"\tname:"<<c_name<<endl;
		  //oss<<"\tq:"<<c_q<<endl;
		}
	  return oss.str();
  }
  string Fastq::toFasta() const
  {
	  return toString(true);
  }

  //bool operator < (const SequenceString& other) const;

  void Fastq::Serialize(ofstream& _ofs)const
  {
	  //do necessary checking
	  if(!_ofs.is_open())
		{
		  cout<<"**ERROR**: closed file buffer. quit..."<<endl;
		  exit(-1);
		}
	 
	  //serialize the sequence string, seq
	  c_ss.Serialize(_ofs);
	  //serializing the rest......
	  const char* p_char; //pointer used to direct the writing to the file stream
	  //save name 
	  unsigned l=c_q.length()+1;
	  p_char=(char*)&l;
	  _ofs.write(p_char, sizeof(unsigned));
	  p_char=c_q.c_str();
	  _ofs.write(p_char, c_q.length()+1);

  }
  void Fastq::Deserialize(ifstream& _ifs)
  {
	  //do necessary checking
	  if(!_ifs.is_open())
		{
		  cout<<"**ERROR**: closed file buffer. quit..."<<endl;
		  exit(-1);
		}
	   //read sequence string first
	  c_ss.Deserialize(_ifs);
	   
	  //deserializing quality......
	  char* p_char; //pointer used to direct the writing to the file stream
	  //save name
	  unsigned l;//=c_name.length()+1;
	  p_char=(char*)&l;
	  _ifs.read(p_char, sizeof(unsigned));
	  p_char=new char[l];
	  _ifs.read(p_char,l);
	  c_q=string(p_char);
	  delete[] p_char;

	  //set name
	  c_name=c_ss.GetName();
  }