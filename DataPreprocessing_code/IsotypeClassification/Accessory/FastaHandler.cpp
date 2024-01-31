#include "FastaHandler.hpp"
#include <stdlib.h>
#include <string.h>
//#include <stdio.h>
#include <fstream>
#include <zlib.h>
#include "string_ext.hpp"
#include "GzTools.hpp"
#include "../SequenceHandlerCommon.hpp"
//#include ".hpp"

#define MAX 6000

using namespace std;
size_t ReadFasta(const string& _fname, vector<SequenceString>& _seqStrVec, bool toUpper)
{
	//need to know whether this is compressed file;
	//we are using the ".gz" suffix as a criterion for the compressed file.
	//we only do gzip compression.not others.
	bool gzflag=_fname.substr(_fname.size()-3, 3)==".gz";
	
  ifstream ifs_p;//(_fname.c_str());
  FILE* fb=NULL;
  GZ_UTILITY gu{0};
  if(gzflag) //gziped
  {
  //cout<<"....1"<<endl;
	fb=gzOpen_B(_fname, gu, "rb");
	if(!fb)
	{
		cout<<">>>>>>ERROR:input file\""<<_fname<<"\" can not be opened, quit....\n"<<endl;
		return string::npos;
	}
  } 
  else
  {
	 ifs_p.open(_fname);
	if(!ifs_p.is_open())
    {
      cout<<">>>>>>ERROR:the input file \""<<_fname<<"\" can not be opened, quit....\n";
      return string::npos;
    }  
  }
  //so far, we have open files, now we need to read the file first.
  string line; //out_line;

  string gene_info;
  string gene_sequence;
  
  int gene_number=0;
  
  int line_count=0;
  string temp_seq;

  cout<<"read the fasta file........."<<endl;
  cout<<"line#: ";
  //  double d_storage, d_current;
  string temp_string;
  bool ok=true;
  while(ok)//(!ifs_p.eof())
    {
      line_count++;
      if(line_count/5000*5000==line_count)
		{
		  cout<<"..... "<<line_count;
		  cout.flush();
		}
      //found_one=false;
	  //cout<<"start the loop......."<<endl;
	  if(gzflag) //compressed file 
	  {
		//cout<<"\tready to read........"<<endl;
		ok=getline_B(fb, line,gu);
	  }
	  else
		getline(ifs_p, line);
      
	  chomp_ext(line);

      if(line.compare("\n")==0||line.length()==0||line.compare("\t")==0||line.compare("\t\t")==0)
		{
		  cout<<"...read an empty row (skip!)"<<endl;
		   //last, check the file end 
		  if(!gzflag)
			ok=!ifs_p.eof();
		
		  continue;
		}
            //cout<<"here......"<<endl;
      if(line[0]=='>')
		{
		  //this is a new refSeq header line,
		  temp_string=line.substr( 1,line.length());
				  
		  if(gene_number>0) //we already read something, so we have to store it to the vector
			{

			  SequenceString ss(gene_info,to_upper_str(temp_seq));
			  _seqStrVec.push_back(ss);
			  
			}//otherwise, we are doing the first one, so do not bother
		  gene_info=temp_string;
		  temp_seq.clear();
		  gene_number++;
		  
		}
      else
		{//append the seqence to it
		  temp_seq.append(line);
		  //temp_seq.append("\n");
		}

		//at the end we need to check for end of file for regular reading (not gz'ed file)
		if(!gzflag)
			ok=!ifs_p.eof();	
    }
  
  //here we have to update/add the last gene
  //read_gene_sequence(temp_seq, gene_sequence)==-1)
  SequenceString ss(gene_info,to_upper_str(temp_seq));
  _seqStrVec.push_back(ss);
  

  //we are done with reading the promoter sequence file.
  cout<<"\nfinish reading the file........\n"
      <<"\tsummary: total "<<line_count<<" line read in and \n"
      <<"\t\t"<<gene_number<<" sequences store in the vector...."<<endl;
  
  if(gzflag)
  {
	gzClose_B(fb,gu);
  }
  else
  {  
	ifs_p.close();
  }
  
  return gene_number;
}

//in here we read ONE record of fasta 
//the caller needs to open/close  the file
/*input: 
 *	fb FILE stream pointer, assuming it is ok. but we need to check for the end of file status.
 *	ss SequenceString reference, used to return the read record.
 *	hangover, string buffer used to hold the hangover from the previous read. 
 *		will be passed between calls.
 *  gu, GZ_UTILITY, used by the gz tools for reading the compressed files. 
 *		will be passed between calls
 *	gz bool, used to indicating whether this is a compressed stream. 
 *output:
 *	return the bool indcating whether the reading COULD BE CONTINUE.
 *		true: meaning we can go on to read and current record is good. 
 *		false when we reach the end or there is an error (either stream reading error or parsing error (no '>', etc).
 *		the outer caller need to check whether it is a eof or error. 
 *      Most importantly, it could happen that the reading can not go on (false), but the current 
 *		read is still valid. Just check to see whether it is a good read. if the read is no good
 *		due to the error, the rescord is set to be empty.!!!!
 */
bool get1FastaSeq(FILE* fb, SequenceString& ss, string& hangover, GZ_UTILITY& gu, const bool& gz)
{
	bool ok=true;
	string line;
	string sname;
	string sseq;
	bool gotHeader=false;
	char buf[MAX];
	//GZ_UTILITY gu;
	
	//now the first thing to do is to check the hangover buffer 
	//to get the reads from the previous read.
	if(hangover.length()!=0)
	{//there are something in there,
		if(hangover[0]!='>')
		{
			//this is impossible. it should not happen.
			cout<<"Corrupted file reading, reading stopped!!"<<endl;
			return false;
		}
		sname=move(hangover.substr( 1));
		gotHeader=true;
	}
	
	//assume the file is opened. start reading.
	while (ok)
	{
		/*
		//need to check whether we are reaching the end of file
		if((gz&&!ok)||((!gz)&&feof(fb)))
		{
			//cout<<"***\t first line reaching end...."<<endl;
			//reaching the end, done. 
			return false;
		}
		*/
		if(gz)
		{
			//cout<<"\tready to read........"<<endl;
			ok=getline_B(fb, line,gu);
			//if(!ok)
			//	cout<<"\t----no good"<<endl;
		}
		else
		{
			//cout<<"read the fline"<<endl;
			char* st=fgets(buf,(size_t)MAX, fb);
			//ok=(feof(fb)&&st!=nullptr);
			ok=(st!=nullptr);
		}  

		if(!ok)
		{//could be error or reaching the end.
			 if(sseq.length()==0 ||sname.length()==0)
			{
				//cout<<"sname:"<<sname<<"\n sseq:"<<sseq<<endl;
				ss.SetName("");
				ss.SetSequence("");
				ok= false;
			}
			else{
				ss.SetName(move(sname));
				ss.SetSequence(move(sseq));
				ok=false;
				
			}
				
			  //return false;
			  //cout<<"\t\t 0--0-breaking the loop inside....."<<endl;
			  break; 
		}
		if(!gz)
			line=buf;
		
		chomp_ext(line);
		
		if(line[0]=='>') //first line or header line, need to check for whether to continue or finish this one.
		{
			if(!gotHeader)
			{
				sname=move(line.substr( 1,line.length()));
				gotHeader=true;
			}
			else //need to finish this	
			{
				ss.SetName(move(sname));
				ss.SetSequence(move(sseq));
	
				hangover=move(line);
				break;
			}
				
		}
		else 
		{
			sseq.append(line);
		}
		/*//read again 
		if((gz&&!ok)||((!gz)&&feof(fb)))
		{
			//cout<<"\t****reaching end"<<endl; 
			//reaching the end, done. 
			return false;
		}
		
		if(gz)
		{
			//cout<<"\tready to read........"<<endl;
			ok=getline_B(fb, line,gu);
		}
		else
		{
			char* st=fgets(buf,MAX, fb);
			ok=(st!=nullptr);
		}  

		if(!ok||(gz&&line[0]=='>')||((!gz)&&buf[0]=='>'))
		{
			  //cout<<"...rerror occured , !!!"<<endl;
			  
			  return false;
		}
		if(gz)
			sseq=move(line);
		else
			sseq=buf;*/
	}
	
	return ok;
}	

//now we want to read and concate the fasta files 
//
size_t concatenateFasta(const string& fn_r1, const string& fn_r2, 
			FileType ft, string outFile_name, 
			ios_base::openmode mode,
			bool rc)
{
		//need to know whether this is compressed file;
	//we are using the ".gz" suffix as a criterion for the compressed file.
	//we only do gzip compression.not others.
	bool gzflag=(ft==GZ_FASTA);
	
	//opent file first
	
	FILE* fb1=NULL; FILE* fb2=NULL;
	GZ_UTILITY gu1{0};
	GZ_UTILITY gu2{0};
  if(gzflag) //gziped
  {
  //cout<<"....1"<<endl;
	fb1=gzOpen_B(fn_r1, gu1, "rb");
	fb2=gzOpen_B(fn_r2, gu2, "rb");
	if(!fb1||!fb2)
	{
		cout<<">>>>>>ERROR:input file\""<<fn_r1<<" or "<<fn_r2<<"\" can not be opened, quit....\n"<<endl;
		return string::npos;
	}
  } 
  else
  {
	 fb1=fopen(fn_r1.c_str(), "r");
	 fb2=fopen(fn_r2.c_str(), "r");
	if(fb1==nullptr||fb2==nullptr)
    {
      cout<<">>>>>>ERROR:the input file \""<<fn_r1<<"\" or \""<<fn_r2<<"\" can not be opened, quit....\n";
      return string::npos;
    }  
  }
  //call to read
	bool ok1=true, ok2=true;
	string hangover1, hangover2;
	SequenceString ss1;
	SequenceString ss2;
	SequenceString ssc;
	unsigned count=0;
	
	ofstream ofs((outFile_name).c_str(), mode);
  
  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<outFile_name<<"\" can not be opened, quit....\n";
      return 0;
    }
  while(ok1&&ok2)
  {
	  ok1=get1FastaSeq(fb1, ss1,hangover1,gu1, gzflag);
	  
	  ok2=get1FastaSeq(fb2, ss2,hangover2, gu2, gzflag);
	  
	  
	  //now we should do concate
	  if(ss1.GetName()!=""&&ss2.GetName()!="")
	  {
		count++;
		/*cout<<"-----seq:"<<count<<endl;
		  cout<<"record 1:"<<ss1.toString()<<endl;
		  cout<<"record 2:"<<ss2.toString()<<endl;
		*/
		//
		ssc.SetName(ss1.GetName());
		string tempSeq=ss1.GetSequence();
		if(rc)
		{
			//only reverse complement the read2
			
			tempSeq.append(ReverseComplement(ss2).GetSequence());
			
		}
		else
		{
			tempSeq.append(ss2.GetSequence());
		}
		ssc.SetSequence(move(tempSeq));
		writ1FastaSeq(ofs, ssc);
	  }
	  
  }
  
  if(gzflag)
  {
	gzClose_B(fb1,gu1);
	gzClose_B(fb2,gu2);
  }
  else
  {  
	fclose(fb1);
	fclose(fb2);
  }
  ofs.close();
  return count;
}

//write 
void WriteFasta(const string& _fname, vector<SequenceString>& _seqStrVec, const unsigned int& _width,  ios_base::openmode mode)
{
  ofstream ofs((_fname).c_str(), mode);
  
  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<_fname<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  
  for(unsigned int i=0;i<_seqStrVec.size();i++)
    {
      ofs<<">"<<_seqStrVec.at(i).GetName()<<"\n";
      for(unsigned int j=0;j*_width<_seqStrVec.at(i).GetLength();j++)
	{
	  ofs<<_seqStrVec.at(i).GetSequence().substr(j*_width, _width)<<"\n";
	}
    }

  ofs.close();
}

void WriteTextFile(const string& _fname, vector<unsigned int>& _seqStrVec, const char& c, const unsigned int& _width, ios_base::openmode mode)
{
  ofstream ofs((_fname).c_str(), mode);
  
  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<_fname<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  
  for(unsigned int i=0;i<_seqStrVec.size();i++)
    {
      //ofs<<">"<<_seqStrVec.at(i).GetName()<<"\n";
      ofs<<_seqStrVec.at(i);
      if((i+1)%_width==0||i==_seqStrVec.size()-1)
	{
	  ofs<<"\n";
	}
      else
	{
	  ofs<<c;
	}
    }

  ofs.close();
  
}

//Writing a text table file, with header or not
//input: _fname, file name
//       _seqStrVec, vector of vector of  numbers to be written
//       _c, a char to delimite the columns
//       _header, whether to write the header
//       mode, to open the files, truc (new) or append
//       _headerStr, the header names to be written.
void WriteTextTableFile(const string& _fname, vector<vector<double> >& _seqStrVec, const char& c, const bool& _header, ios_base::openmode mode, vector<string> _headerStr)
{
  ofstream ofs((_fname).c_str(), mode);
  
  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<_fname<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  unsigned colNum=_seqStrVec.size();
  if(colNum==0)
    {
      ofs.close();
      return ;
    }
  //now check for the headers
  if(_header)//write header, so go head
    {
      //cout<<"Wring header line......."<<endl;
      if(_headerStr.size()==0)
	{
	  //cout<<"size is zero"<<endl;
	  //_headerStr=new vector<string>;
	
	  for(unsigned int i=0;i<colNum;i++)
	    {
	      _headerStr.push_back("NA");
	    }
	}
      if(_headerStr.size()<colNum)
	{
	  //cout<<"heaer line is smaller"<<endl;
	  for(unsigned int i=_headerStr.size();i<colNum;i++)
	    {
	      _headerStr.push_back("NA");
	    }
	}
      for(unsigned int i=0;i<colNum;i++)
	{
	  //cout<<"writing "<<i<<"...";
	  ofs<<_headerStr.at(i);
	  //cout<<"elsment :"<<_headerStr.at(i)<<"xxx";
	  if(i!=colNum-1)
	    ofs<<c;
	  else
	    ofs<<"\n";
	}
      //cout<<endl;
    }
  

  unsigned maxSize=_seqStrVec.at(0).size();
  //cout<<"maxsize is :"<<maxSize<<endl;
  for(unsigned int i=1;i<colNum;i++)
    {
      
      if(maxSize<_seqStrVec.at(i).size())
	{
	  maxSize=_seqStrVec.at(i).size();
	}
	
    }
  //cout<<"maxsize is (after):"<<maxSize<<endl;

  vector<double> tempStr;
  for(unsigned int j=0;j<maxSize;j++)
    {
      //ofs<<">"<<_seqStrVec.at(i).GetName()<<"\n";
      
      for(unsigned int i=0;i<_seqStrVec.size();i++)
	{
	  tempStr=_seqStrVec.at(i);
	  //in case the vectors are of unequal length
	  if(j<tempStr.size())
	    {
	      ofs<<tempStr.at(j);
	    }
	  else
	    {
	      ofs<<"";
	    }
	  if(i!=colNum-1)
	    ofs<<c;
	  else
	    ofs<<"\n";
	}
      
	
	/*if((i+1)%_width==0||i==_seqStrVec.size()-1)
	{
	  ofs<<"\n";
	}
      else
	{
	  ofs<<c;
	  }*/
    }

  ofs.close();
}//end of function



//Writing a text table file, with header or not
//input: _fname, file name
//       _seqStrVec, vector of vector of  numbers to be written
//       _c, a char to delimite the columns
//       _header, whether to write the header
//       mode, to open the files, truc (new) or append
//       _headerStr, the header names to be written.
void WriteTextTableFile(const string& _fname, vector<vector<string> >& _seqStrVec, const char& c, const bool& _header, ios_base::openmode mode, vector<string> _headerStr)
{
  ofstream ofs((_fname).c_str(), mode);
  
  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<_fname<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  unsigned colNum=_seqStrVec.size();
  if(colNum==0)
    {
      ofs.close();
      return ;
    }
  //now check for the headers
  if(_header)//write header, so go head
    {
      //cout<<"Wring header line......."<<endl;
      if(_headerStr.size()==0)
	{
	  //cout<<"size is zero"<<endl;
	  //_headerStr=new vector<string>;
	
	  for(unsigned int i=0;i<colNum;i++)
	    {
	      _headerStr.push_back("NA");
	    }
	}
      if(_headerStr.size()<colNum)
		{
		  //cout<<"heaer line is smaller"<<endl;
		  for(unsigned int i=_headerStr.size();i<colNum;i++)
			{
			  _headerStr.push_back("NA");
			}
		}
      for(unsigned int i=0;i<colNum;i++)
		{
		  //cout<<"writing "<<i<<"...";
		  ofs<<_headerStr.at(i);
		  //cout<<"elsment :"<<_headerStr.at(i)<<"xxx";
		  if(i!=colNum-1)
			ofs<<c;
		  else
			ofs<<"\n";
		}
      //cout<<endl;
    }
  

  unsigned maxSize=_seqStrVec.at(0).size();
  //cout<<"maxsize is :"<<maxSize<<endl;
  for(unsigned int i=1;i<colNum;i++)
    {
      
      if(maxSize<_seqStrVec.at(i).size())
		{
		  maxSize=_seqStrVec.at(i).size();
		}
		
    }
  //cout<<"maxsize is (after):"<<maxSize<<endl;

  //vector<string> tempStr;
  string sw;
  for(unsigned int j=0;j<maxSize;j++)
    {
      //ofs<<">"<<_seqStrVec.at(i).GetName()<<"\n";
      //sw.assign("");
      for(unsigned int i=0;i<_seqStrVec.size();i++)
		{
		  //tempStr=_seqStrVec.at(i);
		  //in case the vectors are of unequal length
		  if(j<_seqStrVec.at(i).size())
			{
			  //ofs<<tempStr.at(j);
			  sw.append(_seqStrVec.at(i).at(j));
			}
		  else
			{
			  //ofs<<"";
			  sw.append("");
			}
		  if(i!=colNum-1)
		  {
			//ofs<<c;
			sw.append(1, c);
		  }
		  else
		  {
			//ofs<<"\n";
			sw.append("\n");
		  }
		}
	  if((j%10000==0&&j>0)||j==maxSize-1)
	  {
		  ofs<<sw;
		  sw.assign("");
	  }
		
	/*if((i+1)%_width==0||i==_seqStrVec.size()-1)
	{
	  ofs<<"\n";
	}
      else
	{
	  ofs<<c;
	  }*/
    }

  ofs.close();
}//end of function

bool writ1FastaSeq(ofstream& ofs, const SequenceString& ss, const unsigned& _width)
{
	ofs<<">"<<ss.GetName()<<"\n";
    for(unsigned int j=0;j*_width<ss.GetLength();j++)
	{
	  ofs<<ss.GetSequence().substr(j*_width, _width)<<"\n";
	}
	return true;
}

