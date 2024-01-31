#include "FastqHandler.hpp"
#include <stdlib.h>
//#include <stdio.h>
#include <fstream>
#include <zlib.h>
#include "string_ext.hpp"
#define MAX 5000

using namespace std;

//not very efficient. might need to rewrite. but how to avoid the code duplication?
size_t ReadFastq(const string& _fname, vector<SequenceString>& _seqStrVec, vector<string>& _vecQ, bool toUpper)
{
	
	cout<<"=============reading 333++++++++++"<<endl;
		vector<Fastq> v_fastq;
		size_t n=ReadFastq(_fname, v_fastq, toUpper);
		_seqStrVec.clear();
		_vecQ.clear();
		if(v_fastq.size()>0)
		{
			for(unsigned i =0;i<v_fastq.size();i++)
			{
				_seqStrVec.push_back(v_fastq.at(i).GetSequenceString());
				_vecQ.push_back(v_fastq.at(i).GetQualityString());
			}
		}
		return n;
}
//return string::npos upon error 
size_t ReadFastq(const string& _fname, vector<Fastq>& _seqStrVec, bool toUpper)
{
	//cout<<"=============reading 2++++++++++"<<endl;
	//need to know whether this is compressed file;
	//we are using the ".gz" suffix as a criterion for the compressed file.
	//we only do gzip compression.not others.
	bool gzflag=_fname.substr(_fname.size()-3, 3)==".gz";
	ifstream ifs_p;
	
	FILE* fb=nullptr;
	
	GZ_UTILITY gu{0};
	gu.z_stream_init=false;
	
	if(gzflag)	
	{	
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
		//cout<<"......2"<<endl;
		if(!ifs_p.is_open())
		{
			cout<<">>>>>>ERROR:the input file \""<<_fname<<"\" can not be opened, quit....\n"<<endl;
			return string::npos;
		}
	}
	
	//cout<<"....4 : ifs good :"<<ifs_p.good()<<endl;
	//cout<<".....3:ifs_p_open:"<<ifs_p.is_open()<<endl;
  //so far, we have open files, now we need to read the file first.
  string line; //out_line;

  string gene_info;
  string gene_sequence;
  string quality;
  int gene_number=0;
  
  int line_count=0;
  string temp_seq;

  //cout<<"read the fastq file........."<<line<<line_count<<gene_number<<endl;
  //cout<<"line#: ";
  //  double d_storage, d_current;
  string temp_string;
  bool ok=true;
  while(ok)
    {
      
	  //start reading lines and parsing
	  //fq files line 1: name
	  //	line 2: seqence
	  //	line 3: +
	  //	line 4: quality
	  //
	  //cout<<"start the loop......."<<endl;
	  //break;
	  if(gzflag) //compressed file 
	  {
		//cout<<"\tready to read........"<<endl;
		ok=getline_B(fb, line, gu);
	  }
	  else
		getline(ifs_p, line);  
      //break;
	  chomp_ext(line);
	  //break;
	  //cout<<"i:"<<line_count<<"--line:"<<line<<endl;

      if(line.compare("\n")==0||line.length()==0||line.compare("\t")==0||line.compare("\t\t")==0)
		{
		  cout<<"...read an empty row (skip!)"<<endl;
		  //last, check the file end 
		  if(!gzflag)
			ok=!ifs_p.eof();
		  continue;
		}
        //first line has to be "@xxxxxxx    
		//cout<<"here......"<<endl;
      if(line[0]=='@')
		{
		  //this is a new refSeq header line,
		  gene_info=line.substr( 1,line.length());
		  gene_number++;
		}
      else
		{//
			cout<<"ERROR: corrupted file, no lines with '@...' please check! Reading is stopped"<<endl;
			/*if(gzflag)
			  {
				gzClose_B(fb);
			  }
			  else
			  {  
				ifs_p.close();
			  }*/
			break;
			
		}
		
	//last, check the file end
		if(!gzflag)	
			ok=!ifs_p.eof();
		
	  //now doing the line 2
	  //sequence
	  if(gzflag) //compressed file 
	  {
		ok=getline_B(fb, line,gu);
	  }
	  else
		getline(ifs_p, line);  
      
	  chomp_ext(line);
	  if(toUpper)
		gene_sequence.assign(to_upper_str(line));
	  else
		gene_sequence.assign(line);

	//cout<<"\t"<<gene_sequence<<endl;
	//last, check the file end 
		if(!gzflag)
			ok=!ifs_p.eof();
		
	  //reading the third line
	  if(gzflag) //compressed file 
	  {
		ok=getline_B(fb, line,gu);
	  }
	  else
		getline(ifs_p, line);  
	
	  chomp_ext(line);
	  
	  if(line[0]!='+'||line.size()!=1)
		{
		  //there is something wrong. we need to quit
		  cout<<"ERROR: something wrong. the file is corrupted!! Reading is stopped!!"<<endl;
		  break;
    	}
	//cout<<line<<endl;
	//last, check the file end 
		if(!gzflag)
			ok=!ifs_p.eof();	
	
	//reading the 4th line
		//reading the third line
	  if(gzflag) //compressed file 
	  {
		ok=getline_B(fb, line,gu);
	  }
	  else
	  {
		getline(ifs_p, line);  
	  }
	  chomp_ext(line);
	  
	  quality.assign(line);
		
	 //cout<<line<<endl; 
	
	  //done, now we need to push to the vector 
	  Fastq fq(gene_info,SequenceString(gene_info, gene_sequence),quality);
	  _seqStrVec.push_back(fq);
	  line_count++;
      if(line_count%10000==0)
		{
		  cout<<"..... "<<line_count;//<<endl;
		  cout.flush();
		}
	  //last, check the file end 
		if(!gzflag)
			ok=!ifs_p.eof();   //we don't have to check for gzipped case, since it will return false anyway.
    }

  //we are done with reading the promoter sequence file.
  cout<<"\nfinish reading the file........\n"
      <<"\tsummary: total "<<line_count<<" records read in and \n"
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

void WriteFastq(const string& _fname, vector<Fastq>& _fqVec,  ios_base::openmode mode)
{
  ofstream ofs((_fname).c_str(), mode);
  
  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<_fname<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  
  for(unsigned int i=0;i<_fqVec.size();i++)
    {
      ofs<<"@"<<_fqVec.at(i).GetName()<<"\n";
      
	  ofs<<_fqVec.at(i).GetSequence()<<"\n";
	  
	  ofs<<"+\n"<<_fqVec.at(i).GetQualityString()<<"\n";
	
    }

  ofs.close();
}

void WriteFastq(const string& _fname, vector<SequenceString>& _seqStrVec, 
	vector<string>& _qVec,
	ios_base::openmode mode)
{
  ofstream ofs((_fname).c_str(), mode);
  
  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<_fname<<"\" can not be opened, quit....\n";
      exit(-1);
    }
  if(_qVec.size()!=_seqStrVec.size())
  {
	  cout<<"ERROR: the sizes of the sequence string and qualtiy string don't match. please check........"<<endl;
  }
  for(unsigned int i=0;i<_qVec.size();i++)
    {
      ofs<<"@"<<_seqStrVec.at(i).GetName()<<"\n";
      
	  ofs<<_seqStrVec.at(i).GetSequence()<<"\n";
	  
	  ofs<<"+\n"<<_qVec.at(i)<<"\n";
	
    }

  ofs.close();	
		
}
/*
 *NOTE: IMPORT need to check the return bool, false means can not CONTINUE, but the record might be good.
 */
bool get1FastqSeq(FILE* fb, SequenceString& ss, string& qs, GZ_UTILITY& gu, const bool& gz)
{
	bool ok=true;
	string line;
	string sname;
	string sseq;
	
	char buf[MAX];
	
	//assume the file is opened. start reading.
	//cout<<"read 1"<<endl;
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
			//cout<<"in here ,error"<<endl;
			 if(sseq.length()==0 ||sname.length()==0||qs.length()==0)
			{
				//cout<<"sname:"<<sname<<"\n sseq:"<<sseq<<endl;
				ss.SetName("");
				ss.SetSequence("");
				qs.clear();
				ok= false;
			}
			else
			{
				ss.SetName(move(sname));
				ss.SetSequence(move(sseq));
				
				ok=false;
				
			}			
			  //return false;
			  //cout<<"\t\t 0--0-breaking the loop inside....."<<endl;
			  return ok; 
		}
		if(!gz)
			line=buf;
		//cout<<"read 2"<<endl;
		chomp_ext(line);
		//cout<<"line read:"<<endl;
		if(line[0]=='@') //first line or header line, need to check for whether to continue or finish this one.
		{
			ss.SetName(move(line.substr( 1,line.length())));
				//gotHeader=true;	
		}
		else 
		{
			//cout<<"in here error"<<endl;
			ok=false; 
			return ok;
			//ss.SetSequence(line);
		}
		//cout<<"read again......."<<endl;

		//read again , second line for sequence
		if(!ok)
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
			line=buf;
		}  
		chomp_ext(line);
		if(!ok)//||((!gz)&&buf[0]=='>'))
		{
			  //cout<<"...rerror occured , !!!"<<endl;
			  
			  return false;
		}
		ss.SetSequence(line);
		//cout<<"read 3"<<endl;
		//read again again
		//cout<<"read again and again ......."<<endl;
		//read again , second line for sequence
		if(!ok)//||((!gz)&&feof(fb)))
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
			line=buf;
		}  
		chomp_ext(line);
		if(!ok||line.length()!=1||line[0]!='+')//||((!gz)&&buf[0]=='>'))
		{
			  cout<<"...rerror occured , !!!"<<endl;
			  
			  return false;
		}
		//don't do anything since this is '+' line
		//cout<<"read 4"<<endl;
		//read again and again and again for quality string 
		if(gz)
		{
			//cout<<"\tready to read........"<<endl;
			ok=getline_B(fb, line,gu);
		}
		else
		{
			char* st=fgets(buf,MAX, fb);
			ok=(st!=nullptr);
			line=buf;
		}
		if(!ok)//||((!gz)&&feof(fb)))
		{
			//cout<<"\t****reaching end"<<endl; 
			//reaching the end, done. 
			return false;
		}
		chomp_ext(line);
		qs.assign(line);
		//cout<<"---header:"<<ss.toString()<<endl;
	
	
	return ok;
}

//now we want to read and concate the fasta files 
//
size_t concatenateFastq(const string& fn_r1, const string& fn_r2, 
		FileType ft,string outFile_name, 
		ios_base::openmode mode,
		bool rc )
{
	//need to know whether this is compressed file;
	//we are using the ".gz" suffix as a criterion for the compressed file.
	//we only do gzip compression.not others.
	bool gzflag=(ft==GZ_FASTQ);
	
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
	string qs1, qs2;
	SequenceString ss1;
	SequenceString ss2;
	SequenceString ssc;
	string qsc;
	unsigned count=0;
	
	ofstream ofs((outFile_name).c_str(), mode);
  
  if(!ofs.is_open())
    {
      cout<<">>>>>>ERROR:the output file \""<<outFile_name<<"\" can not be opened, quit....\n";
      return 0;
    }
	//cout<<"start doing it ...."<<endl;
  while(ok1&&ok2)
  {
	  ok1=get1FastqSeq(fb1, ss1,qs1,gu1, gzflag);
	  
	  ok2=get1FastqSeq(fb2, ss2,qs2, gu2, gzflag);
	  
	  string tempStr;
	  //cout<<"ss1:"<<ss1.GetName()<<endl;
	  //now we should do concate
	  if(ss1.GetName()!=""&&ss2.GetName()!="")
	  {
		  //cout<<"Get one !!!"<<endl;
		  ssc.SetName(ss1.GetName());
		  qsc.assign(qs1);
		  
		  tempStr.assign(ss1.GetSequence());
		  count++;
		  if(rc)
		  {
			  tempStr.append(ReverseComplement(ss2).GetSequence());
			  qsc.append(move(flipStr(qs2)));
		  }
		  else 
		  {
			  tempStr.append(ss2.GetSequence());
			  qsc.append(move(qs2));
		  }
		  ssc.SetSequence(move(tempStr));
		  writ1FastqSeq(ofs, ssc, qsc);
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

bool writ1FastqSeq(ofstream& ofs, const SequenceString& ss, const string& qs)
{
	ofs<<"@"<<ss.GetName()<<"\n";
    ofs<<ss.GetSequence()<<"\n";
	ofs<<"+\n";
	ofs<<qs<<"\n";
	
	return true;
}