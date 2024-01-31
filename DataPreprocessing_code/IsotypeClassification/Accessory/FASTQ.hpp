#ifndef FASTQ_HPP
#define FASTQ_HPP

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "SequenceString.hpp"

using namespace std;

//in this module , we define the class of data structure of FASTQ file format

//define each individual fastq sequence 

class Fastq
{
public:
	//fastq();
	Fastq(const string& _seqName="", const SequenceString& _ss=SequenceString(), const string& _q="" );
	~Fastq();//destructor empty one

  //no copy constructor so far

  void SetName(const string& _name);
  void SetSequenceString(const SequenceString& _seq);
  void SetQualityString(const string& _q);
  
  const string GetName()const;
  const string GetSequence()const;
  const SequenceString GetSequenceString()const;

  const string GetQualityString() const;

  //const unsigned GetLength()const;

  //unsigned GetLetterCount(const char& _c) const; 

  string toString(bool _fasta=false) const;
  string toFasta() const;
  //bool operator < (const SequenceString& other) const;

  void Serialize(ofstream& _ofs)const;
  void Deserialize(ifstream& _ifs);
protected:
  string c_name; //name, the first line 
  SequenceString c_ss;//sequencestring holding the sequence 
  string c_q;//qualtiy score
};	

#endif