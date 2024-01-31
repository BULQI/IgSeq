#include "FileHandler.hpp"
#include "string_ext.hpp"
#include "FASTQ.hpp"
#include "FastaHandler.hpp"
#include "FastqHandler.hpp"
#include <string>
#include <string.h>

#define MAX 5000

bool is_file(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISREG(buf.st_mode);
}

string basename (const std::string& str)
{
  //std::cout << "Splitting: " << str << '\n';
  std::size_t found = str.find_last_of("/\\");
  //std::cout << " path: " << str.substr(0,found) << '\n';
  //std::cout << " file: " << str.substr(found+1) << '\n';
  return str.substr(found+1);
}


bool is_dir(const char* path) {
    struct stat buf;
    stat(path, &buf);
    return S_ISDIR(buf.st_mode);
}

bool exist(const char* path) {
    struct stat buf;
	
    int x= stat(path, &buf);
	//cout<<"inside the exist:x="<<x<<endl;
    return x==0;
}
//use an open file stream to check whther theis is a text file
bool is_text(const string& fname)
{
	FILE* f=fopen(fname.c_str(),"r");
	if(f==nullptr)
	{
		cout<<"ERROR:can not open the file........."<<endl;
		return false;
	}
	bool flag=true;
	char c1=0;
	int count=0;
	while(count <=20) //keep going for ten of them .....
	{
		if(feof(f)) //reaching the end  
		{
			fclose(f);
			return true; 
		}
		//treat this as a
		int numRead=fread(&c1, sizeof(char), (size_t)1,f);
		if(numRead!=1)
		{
			cout<<"ERROR:  an error/corruption occured......"<<endl;
			//fclose(f);
			flag=false;
		}
		if(((int)c1>127||(int)c1<32)
						&&(int)c1!=9
						&&(int)c1!=10
						&&(int)c1!=13)
				{
					flag=false;
					break;
				}
		count++;
	}
	fclose(f);
	return flag;
}
bool is_fasta(const string& fname)
{
	//assuming this is a text file 
	//first reset the stream 
	FILE* f=fopen(fname.c_str(),"r");
	if(f==nullptr)
	{
		cout<<"ERROR:can not open the file........."<<endl;
		return false;
	}
	
	//now read first 6 lines if they are availale
	//#int MAX=5000;
	char buf[MAX];
	int lcount=0;
	unsigned rcount=0;
	bool flag=true;
	string hangover;
	bool gotHeader=false, gotSeq=false;
	//read the first line to test for '>'
	char* r=fgets(buf, MAX, f);
	if(r ==nullptr)
	{
		flag=false;
		return flag;
	}
	if(feof(f))
	{
		//fclose(f);
		flag=true;
		if(lcount==0)
			flag=false ;//an empty file.
		return flag;
	}
	//this is the first line, we need a 
	if(buf[0]!='>')
	{
		flag=false;
		return flag;
	}
	else
	{
		gotHeader=true;
	}
	
//start reading the rest 
	while(rcount<4)
	{
		if(feof(f))
		{
			//fclose(f);
			flag=(gotHeader&&gotSeq);
			break;
		}
		
		r=fgets(buf, MAX, f);
		if(r ==nullptr)
		{
			flag=false;
			break;
		}
		//this is the first line, we need a 
		if(buf[0]!='>')
		{
			gotSeq=true;
			if(lcount>100)//too long now.
			{
				flag=false;
				break;
			}
		}
		else
		{
			rcount++;
			lcount=0;//reset
			gotHeader=true;
		}
	}//end of while
	fclose(f);
	return flag;
}

bool is_fastq(const string& fname)
{
	//assuming this is a text file 
	//first reset the stream 
	FILE* f=fopen(fname.c_str(),"r");
	if(f==nullptr)
	{
		cout<<"ERROR:can not open the file........."<<endl;
		return false;
	}
	
	//now read first 6 lines if they are availale
	//#int MAX=5000;
	char buf[MAX];
	int lcount=0;
	bool flag=true;
	while(lcount<16)
	{
		//cout<<"--------line count:"<<lcount<<endl;
		//cout<<"\treading 1sth line:"<<feof(f)<<endl; 
		if(feof(f))
		{
			flag=true;
			if(lcount==0)
				flag=false ;//an empty file.
			break;
		}
		char* r=fgets(buf, MAX, f);
		if(r ==nullptr)
		{
			flag=false;
			break;
		}
		//this is the first line, we need a 
		if(buf[0]!='@')
		{
			flag=false;
			break;
		}
		//read a second line 
		//cout<<"\treading 2nd line:"<<endl;
		if(feof(f))
		{
			flag=false;
			break;

		}
		r=fgets(buf, MAX, f);
		if(r ==nullptr)
		{
			flag=false;
			break;
		}
		//this is not the first line, we need a sequence line  
		if(buf[0]=='@')
		{
			flag=false;
			break;
		}
		
		//read a third line 
		//cout<<"\treading 3rd line:"<<endl;
		if(feof(f))
		{
			flag=false;
			break;
		}
		r=fgets(buf, MAX, f);
		if(r ==nullptr)
		{
			flag=false;
			break;
		}
		//this is the third line, we need a '+' 
		if(buf[0]!='+')
		{
			flag=false;
			break;
		}
		
		//read a forth line 
		//cout<<"\treading 4th line:"<<endl;
		if(feof(f))
		{
			flag=false;
			break;
		}
		r=fgets(buf, MAX, f);
		if(r ==nullptr)
		{
			flag=false;
			break;
		}
		//this is the fourth line, we need a 
		/*if(buf[0]!='+')
		{
			flag=false;
			break;
		}*/
		lcount+=4;
	}//end of while
	
	fclose(f);
	return flag;
}

FileType check_gzFileType(const string& fname)
{
	GZ_UTILITY gu{0};
	//assuming this is a 
	FILE* fb=gzOpen_B(fname, gu, "rb");
	
	if(!fb)
	{
		cout<<">>>>>>ERROR:input file\""<<fname<<"\" can not be opened, quit....\n"<<endl;
		return UNKNOWN;
	}
	
	//
	string line;
	bool ok=true;
	FileType ft=UNKNOWN;
	int lcount=0;
	FileType ft_previous=UNKNOWN;
	while(ok)
	{
		//cout<<"---line:"<<lcount<<endl;
		if(lcount>16)
			break;
		
		ok=getline_B(fb, line,gu);
		//cout<<"\tread line 1sth:"<<line<<endl;
		if(!ok)
		{
			ft=UNKNOWN;
			break; 
		}
		if(line.length()<1) //don't allow empty line 
		{
			ft=GZ;
			break;
		}
		if(line.length()>1&&line[0]=='@')
		{
			if(ft_previous!=GZ_FASTQ&&ft_previous!=UNKNOWN)
			{
				ft=GZ;
				break;
			}
			ft=GZ_FASTQ;
		}
		if(line.length()>1&&line[0]=='>')
		{
			if(ft_previous!=GZ_FASTA&&ft_previous!=UNKNOWN)
			{
				ft=GZ;
				break;
			}
			ft=GZ_FASTA;
		}			
		if(line.length()>1&&line[0]!='@'&&line[0]!='>')
		{
			ft=GZ;
			break;
		}
		
		//reading the second line
		ok=getline_B(fb, line,gu);
		//cout<<"\tread line 2nd:"<<line<<endl;
		if(!ok)
		{
			ft=GZ;
			break; 
		}
		if(line.length()<=0)
		{
			ft=GZ;
			break;
		}
		if(ft==GZ_FASTA)
		{
			lcount+=2;
			
			ft_previous=GZ_FASTA;
			//don't read any further. just return assuming this is a good fasta file.
			break;
		}			
		//can only be FASTQ
		//read a third line 
		ok=getline_B(fb, line,gu);
		//cout<<"\tread line 3rd:"<<line<<endl;
		if(!ok)
		{
			ft=GZ;
			break; 
		}
		if(line.length()!=1||line[0]!='+')
		{
			ft=GZ;
			break;
		}
		
		//read the fourth line
		ok=getline_B(fb, line,gu);
		//cout<<"\tread line 4th:"<<line<<endl;
		if(!ok)
		{
			ft=GZ;
			break; 
		}
		if(line.length()<=0)
		{
			ft=GZ;
			break;
		}
		lcount+=4;
		ft=GZ_FASTQ;
		ft_previous=GZ_FASTQ;
	}
	gzClose_B(fb,gu);
	return ft;
}
//in here we need to read the file and check the first two bytes
//to identify what type this is.
//assuming the file exist
FileType getFileType_deep(const string& fname)
{
	
	const char* fn=fname.c_str();
	
	FILE* f=fopen(fn,"rb"); //binary read mode 
	if(f==nullptr)
	{
		cout<<"ERROR: can not open the file ....."<<endl;
		return UNKNOWN;
	}
	//here it must be a regular file 
	unsigned char c1, c2;
	size_t numRead;
	
	//bool gz_flag=false;
	//bool text_flag=false;
	FileType ft {UNKNOWN};
	//cout<<"*****deteting....."<<endl;
	//do{
		//read the first two bytes 
		numRead=fread(&c1, sizeof(char), (size_t)1,f);
		if(numRead!=1)
		{
			cout<<"ERROR: either we reach the end of file too early or an error occured......"<<endl;
			fclose(f);
			return UNKNOWN;
		}
		numRead=fread(&c2, sizeof(char), (size_t)1,f);
		if(numRead!=1)
		{
			cout<<"ERROR: either we reach the end of file too early or an error occured......"<<endl;
			fclose(f);
			return UNKNOWN;
		}
		fclose(f);
		//check the chars
		switch((int)c1)
		{
			case 31: //
				if((int)c2==139)
				{
					//gz_flag=true;
					ft=check_gzFileType(fname);
				}
				else
				{					
					ft=UNKNOWN;
				}				
				break;
			case 62:
				if(is_fasta(fname))
					ft= FASTA;
				else
				{
					if(is_text(fname))
					{
						ft= TXT;
					}
					else
					{
						ft= UNKNOWN;
					}
				}
				break;
			case 64:
				if(is_fastq(fname))
					ft= FASTQ;
				else
				{
					if(is_text(fname))
					{
						ft= TXT;
					}
					else
					{
						ft= UNKNOWN;
					}
				}
				break;
			default: //otherwise, could be a text, check for regular ascii chars 
				ft=UNKNOWN;
				if(((int)c1<=127&&(int)c1>=32)
						||(int)c1==9
						||(int)c1==10
						||(int)c1==13)
				{
					if(((int)c2<=127&&(int)c2>=32)
						||(int)c2==9
						||(int)c2==10
						||(int)c2==13)
						{
							//check for text format
							if(is_text(fname))
							{
								ft= TXT;
							}
							else
							{
								ft= UNKNOWN;
							}
						}
				}
				
				break;
		}

	//fclose(f);
	return ft;
}

FileType getFileType_byName(const string& fname)
{
	const char* f=fname.c_str();
	if(!exist(f))
		return NOT_EXIST;
	//
	if(is_dir(f))
	{
		return DIR;
	}
	//
	if(!is_file(f))
	{
		return UNKNOWN;
	}
	//here it must be a regular file 
	size_t dot_pos;
	string suffix;
	//bool tryAgain=false;
	size_t size=fname.size()-1;
	bool gz_flag=false;
	
	//cout<<"*****deteting....."<<endl;
	do{
		//cout<<"\ttrying......."<<endl;
		//tryAgain=false;
		//Check for the suffix.
		dot_pos=fname.find_last_of('.',size);
		//cout<<"dot-Pos:"<<dot_pos<<";string:npos"<<string::npos<<endl;
		if(dot_pos==string::npos)
		{
			//can not find it, then 
			if(gz_flag)
				return GZ;
			return UNKNOWN;
		}
		
		
		//suffix without '.'
		suffix=fname.substr(dot_pos+1, size-dot_pos);
		
		suffix=to_lower_str(suffix);
		//cout<<"suffix:"<<suffix<<endl;
		
		size=dot_pos-1;
		//cout<<"size:"<<size<<endl;
		if(suffix.compare("txt")==0)
		{
			//cout<<"here"<<endl;
			if(gz_flag)
				return GZ_TXT;
			return TXT;
		}
		if(suffix.compare("fasta")==0||suffix.compare("fa")==0||suffix.compare("fas")==0)
		{
			if(gz_flag)
				return GZ_FASTA;
			return FASTA;
		}
		if(suffix.compare("fq")==0||suffix.compare("fastq")==0)
		{
			if(gz_flag)
				return GZ_FASTQ;
			return FASTQ;
		}
		
		if(suffix.compare("gz")==0||suffix.compare("gzip")==0)
		{
			
			if(gz_flag)
			{//multi-gz'ed
				cout<<"WARNING: In detecting file type, multiple compression found!!"<<endl;
				return UNKNOWN;
			}
			gz_flag=true;
			//tryAgain=true;
		}
	}	
	while (gz_flag);
	
	return UNKNOWN;
}

//get file type, note FileType is user-defined enum.
FileType getFileType(const string& fname, bool deep)
{
	if(!deep)
		return getFileType_byName(fname);
	else 
		return getFileType_deep(fname);
}


//a file handler to read file into a fasta vector vector
//we try to detect the following thing in this function 
// 1, file exist?
// 2, is it a file or directory
// 3, is it a fasta, fastq or gziped fasta, gziped fastq. No other type supported so far
// 
// We will return a vector holding the squenences of the file (SquenceStrings)
// we will return the total number of sequences read in. If none or error, we will return 
// string::npos.
size_t readFile2SeqStrVector(const string& _fname, vector<SequenceString>& _vec, vector<string>* _vec_Q)
{
	//check for files
	if(!exist(_fname.c_str()))
	{
		cout<<"ERROR: the specified file ("<<_fname<<") does not exist!!"<<endl;
		return string::npos;
	}
	//check for files
	if(is_dir(_fname.c_str()))
	{
		cout<<"ERROR: the specified file ("<<_fname<<") is a directory!!"<<endl;
		return string::npos;
	}
	//check file type 
	FileType ft=getFileType(_fname);
	size_t numReads;
	//cout<<"before loop"<<endl;
	switch(ft)
	{
		case FASTA:
	
		case GZ_FASTA:
			numReads=ReadFasta(_fname, _vec);
			break;
			
		case FASTQ:
		case GZ_FASTQ:
		{
			/*if(_vec_Q==NULL)
			{
				cout<<"ERROR: Quality vector(s) for fastq data has not been initalizated correctly. check !!"<<endl;
				exit(-1);
			}*/
			//cout<<"reading.FASTQ"<<endl;
			vector<Fastq> _vec_fq;
			//cout<<"we are here........."<<endl;
			numReads=ReadFastq(_fname, _vec_fq);
			//cout<<"done reading fastq:"<<numReads<<endl;
			if(numReads==string::npos)
				return numReads;
			//rewrite it into fasta vector
			for(unsigned i=0;i<numReads;i++)
			{
				_vec.push_back(_vec_fq.at(i).GetSequenceString());
				if(_vec_Q!=nullptr)
					_vec_Q->push_back(_vec_fq.at(i).GetQualityString());
			}
			numReads=_vec.size();
			break;
		}	
		case GZ:
		case GZ_TXT:
		case TXT:
		case UNKNOWN:
		case DIR:
		default:
			cout<<"ERROR: unsupported file format for this action"<<endl;
			numReads= string::npos;
			break;
	}
	cout<<"Done inside file read vect"<<endl;
	return _vec.size();
}	


size_t concatnateSeqFiles(const string& fn_r1, const string& fn_r2)
{
	size_t totalNum=string::npos;
	//now start 
	if(!exist(fn_r1.c_str())||!exist(fn_r2.c_str()))
	{
		cout<<"ERROR: file(s) doesn't exist. please check..."<<endl;
		return totalNum;
	}
	//check filetype
	FileType ft1=getFileType(fn_r1,true); //deep mode, open file to check.
	FileType ft2=getFileType(fn_r2,true); //deep mode, open file to check.
		
	if(ft1!=ft2)
	{
		cout<<"ERROR: the two read files are not of identical format. please check........"<<endl;
	}		
	switch(ft1)
	{
		case FASTA:
		case GZ_FASTA:
		//
			break;
		case FASTQ:
		case GZ_FASTQ:
			break;
		default:
			cout<<"ERROR:only fasta/fastq files are supported currently, please check the file format...."<<endl;
			break;
	}
	
	
	
	return 0;
}


