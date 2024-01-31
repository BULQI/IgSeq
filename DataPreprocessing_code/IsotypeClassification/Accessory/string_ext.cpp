#include "string_ext.hpp"
using namespace std;
#include <iostream>
//get rid of space in front of a string
void chomp_front_ext(string& s)
{
  for( unsigned int i=0;i<s.length();)
    {
      if(s[i]==' '||s[i]=='\t'||s[i]=='\n'||s[i]=='\r')
		{
		  s.erase(i,1);
		}
      else
		break;
    }
  //return 0;
}

//we work in-situ. means we get rid of white space in the front and
//shift forward the whole string
//return the new length of the string. or the old lenth if no chomping
unsigned chomp_front_ext(char* s, const unsigned& len /*remember that the real length is one more
							*since there is the null in the end '\0'*/
						)
{
	
	if(len==0)
	{
		return 0;
	}
	//assuming len is set up correctly, non-negative 
	if(s[len]!='\0')
	{
		//not correctly ended. so set it to be correctly terminated
		cout<<"WARNING: the C-style string is not set correctly terminated"<<endl;
		s[len]='\0';
		
	}
	unsigned real_start=0;
	for( ;real_start<len;)
    {
      if(s[real_start]==' '||s[real_start]=='\t'||s[real_start]=='\n'||s[real_start]=='\r')
		{
		  real_start++; //point to the next one
		}
	  else
		break;
    }
	
	if(real_start==0) //we are done here
		return len;
	
	unsigned real_len=len-real_start;
	for (unsigned int i=0; i<=real_len;i++)
	{
		//We need to shift towards front
		s[i]=s[real_start];
		real_start++;
	}
	
   return real_len;
}

void chomp_end_ext(string& s)
{
  for(int i=s.length()-1;i>=0;)
    {
      if(s[i]==' '||s[i]=='\t'||s[i]=='\n'||s[i]=='\r')
	{
	  s.erase(i,1);
	  i=s.length()-1;
	}
      else
	break;
    }
  //return 0;
}

unsigned chomp_end_ext(char* s, const unsigned& len)
{
	if(len==0)
	{
		return 0;
	}
	//assuming len is set up correctly, non-negative 
	if(s[len]!='\0')
	{
		//not correctly ended. so set it to be correctly terminated
		cout<<"WARNING: the C-style string is not set correctly terminated"<<endl;
		s[len]='\0';
		
	}
	unsigned real_end=len; //here is one more, pointing to the null char \0
	for(unsigned int i=len-1;i>=0;i--)
    {
      if(s[i]==' '||s[i]=='\t'||s[i]=='\n'||s[i]=='\r')
		{
			real_end =i;
		}
      else
		break;
    }
	return real_end;
}
void chomp_ext(string & s)
{
  chomp_front_ext(s);
  chomp_end_ext(s);
  //return 0;
}

unsigned chomp_ext(char* s, const unsigned& len)
{
		unsigned real_len;
		real_len=chomp_end_ext(s, len);
		real_len=chomp_front_ext(s,len);
		return real_len;
}

//return the number of elements parsed.5
int split_ext(const string& s, string* buf, const char& delim)
{
  int len=0;
  
  len++;
  int start_p=0;
  //first run to get the number of substrings
  for( unsigned int i=0;i<s.length();i++)
    {
      if(s[i]==delim)
	{
	  len++;  
	}
    } 
  //cout<<"first run, lenght is "<<len<<endl;
  //second run
  //*buf=new string[len]();
  len=0; 
  for( unsigned int i=0;i<s.length();i++)
    {
      if(s[i]==delim)
	{
	  //len++;
	  buf[len]=s.substr(start_p,i-start_p);
	  start_p=i+1;
	  //cout<<"len is "<<len<<";i is "<<i<<"buf now is :"<<(*buf)[len]<<endl;
	  len++;
	}
    }
  //after the run, we have to take care of the rest of string
  buf[len]=s.substr(start_p, s.length()-start_p);
  return len+1;
}

int is_number(const char* str)
{
  int index=0;

  bool whiteSpaceFlag=false; //used to get rid of white space ' ','\t', '\n','\r'
  //white space has to be on the front end, not in between. (for end white space in
  //the end, check for trailingWhiteSpace
  
  int flag=0;//this is the flag to indicate whether the  
  //bit could be '+','-' in addition to numbers and dot '.',
  //this the first bit and it could be '+','-', flag==1;
  //this is not the first bit, it has to be numbers or dot '.' flag==0;
  //no scientific format allowed
  bool trailingWhiteSpace=false;
  int dotflag=0;//used to count dots
  int trueFloat=0;
  bool signFlag=false;//used to indicate whether the sing is in the end.
  
  while(str[index]!='\0')
    {
      if(index==0||whiteSpaceFlag)
	{
	  whiteSpaceFlag=true;
	  flag=1;
	}
      else
	{
	  flag=0;
	  whiteSpaceFlag=false;
	}

      //cout<<"Index is "<<index<<endl;
      switch (str[index])
	{
	case ' ':
	  ;
	case '\r':
	  ;
	case '\n':
	  ;
	case '\t':
	  if(!whiteSpaceFlag)
	    {
	      trailingWhiteSpace=true;
	      break;
	    }
	  break;  
	case '+':
	  ;
	case '-':
	  if(!flag)
	    {
	      return 0;
	    }
	  whiteSpaceFlag=false;
	  trailingWhiteSpace=false;
	  signFlag=true;
	  break;
	case '.':
	  if(dotflag>=1)//there is more than one 'dot',so it is not ok
	    return 0;
	  else
	    dotflag++;
	  break;
	  
	case '1':
	  ;
	case '2':
	  ;
	case '3':
	  ;
	case '4':
	  ;
	case '5':
	  ;
	case '6':
	  ;
	case '7':
	  ;
	case '8':
	  ;
	case '9':
	  if(dotflag==1)//a float one, 
	    //so there is no-zero number trailing the 'dot', 
	    //this is a true float number
	    trueFloat++;
	  ;
	case '0':
	  signFlag=false;
	  whiteSpaceFlag=false;
	  if(trailingWhiteSpace) //following trailing space with a number, not allowed
	    return 0;
	  break;
	default:
	  return 0;

	}
      index++;
    }
  if(index==0)
    {
      //this is an empty string. 
      return 0;
    }
  if(signFlag)//means we have '+' or '-' in the end, we don't allow this at this momnet
    return 0;
  if(dotflag==1&&trueFloat>0)
    {
      return 2;//for flat number
    }
  else
    return 1;//for integer
  return 0;
}


//check whether it is one of the following
//1,2,3,4,5,6,7,8,9,0
bool is_number_char( char c)
{
  bool temp=true;
  switch (c)
    {  
    case '1':
      ;
    case '2':
      ;
    case '3':
      ;
    case '4':
      ;
    case '5':
      ;
    case '6':
      ;
    case '7':
      ;
    case '8':
      ;
    case '9':
      ;
    case '0':
      break;
    default:
      temp=false;
      
    }
  
  
  return temp;
  
}


//Definition:
//input --- string s, contains string of format #-#:# for parsing
//          char1 and char2, the char1 is for consecutive char such as for now using '-'
//                           char2 is for including char ':'
//output --- vector vec, the vector contains the numbers parsed from string
//           int, the whether there is error for the string
//        return 0 of good, return non-zero for error
int parseNumberString(const string& s, vector<int>& vec, char stopping_char1, char stopping_char2)
{
  
  int len=s.size();
  const char* s_char=s.c_str();
  char sub[1000];
  int sub_start=0;//the index for current char to be stored in the sub buffer
  bool consecutive_mode=true;
  bool previous_consecutive_mode=true;
  bool start_mode=true;
  int start_number=-1;
  int end_number=-1;
  bool populate_vec=false;
  bool stop_mode=false;
  bool char_round=false;

  //go through the string to parse
  for(int i=0;i<len;i++)  //the char string is one longer than the original string, with '\0' in the end
    {
      //cout<<"Round i"<<i<<endl;
      //cout<<"count is "<<i<<endl;
      
      if(!(s_char[i]=='-'||s_char[i]==':'||is_number_char(s_char[i])||s_char[i]=='\0'))
	{
	  cout<<"******ERROR in parsing the number strings:\n\t\tthe numeric string \""<<s<<"\"is of incorrect format,quit......\n"; 
	  return -1;
	}
      populate_vec=false;//don't insert numbers in the vector yet
      if(s_char[i]=='-')
	{
	  //cout<<"\t***Con"<<endl;
	  consecutive_mode=true;
	  populate_vec=true;

	}
      if(s_char[i]==':')
	{
	  //cout<<"\t***in Con"<<endl;
	  consecutive_mode=false;
	  populate_vec=true;
	}

      if(i==len-1)//the last one, we have to stop, just like finding '-' or ':', the conditions are the same!!!
	{
	  stop_mode=true;
	  //populate_vec=true; 
	}

      //if we are here then we have to check to see whether we want to insert number/populate
      if(populate_vec||stop_mode)
	{
	  if(char_round&&populate_vec)
	    {
	      cout<<"******ERROR:unknown format for input \n\t"<<s<<"\n";
	      cout<<"\t";
	      for(int k=0;k<i;k++)
		{
		  cout<<" ";
		}
	      cout<<"^\n\n";
	      return -1;
	    }
	  char_round=true;
	  if(stop_mode)
	    {

	      if(!populate_vec)//we stopping at a finish char '-' or ':'
		{
		  sub[sub_start]=s_char[i];
		  sub[sub_start+1]='\0';
		}

	      if(populate_vec)//this means we stop when encountering '-' or ':', 
		             //which has not set up end number yet,
		             //for now, we think this is not appropriate, 
		             //but just issuing a WARNINg, and then leave it alone
		{
		  cout<<"WARNING: the string ends at '"<<stopping_char1<<"' or '"<<stopping_char2<<"', which are ignored...\n";
		}
	    }
	  
	  if(start_mode) //we are doing nothig, since it is a starting number, we have
	                 //to check further to do anything
	    {
	      start_number=atoi(sub);
	      start_mode=false;
	      
	    }
	  else  //not the starting mode, so we have to populate 
	        //the fields based on the previous consecutive mode
	    {
	      end_number=atoi(sub);

	      if(previous_consecutive_mode)
		{
		  if(end_number<start_number)
		    {
		      cout<<"******ERROR in parsing the number strings:\n\t\t please using increasing order for field index\n"<<endl;
		      return -1;
		    }
		  for(int j=start_number;j<end_number;j++)
		    {
		      vec.push_back(j);
		    }
		}
	      else
		{
		  vec.push_back(start_number);
		}
	      
	      //reset for next round
	      start_number=end_number;
	    }
	  
	  if(stop_mode)
	    {
	      vec.push_back(start_number);
	    }
	  
	  previous_consecutive_mode=consecutive_mode;
	  populate_vec=false;

	  //reset the buffer pointer for remembering the numbers
	  sub_start=0;
	  sub[0]='\0';
	    

	}
      else //not populate_vec, then it must a number for this round, so put it in a temp_char
 	{
	  sub[sub_start]=s_char[i];
	  sub_start++;
	  sub[sub_start]='\0';
	  char_round=false;
	}
    }

  return 0;
}

string to_upper_str(const string& s)
{
  string ret_s=s;
  //char temp_c;
  for(unsigned i=0;i<ret_s.size();i++)
    {
      if(ret_s[i]<=122&&ret_s[i]>=97) //lower case [a-z]
	{
	  ret_s[i]=ret_s[i]-32;
	}
    }
  return ret_s;
}
string to_lower_str(const string& s)
{
  string ret_s=s;
  //char temp_c;
  for(unsigned i=0;i<ret_s.size();i++)
    {
      if(ret_s[i]<=90&&ret_s[i]>=65) //upper case [A-Z]
	{
	  ret_s[i]=ret_s[i]+32;
	}
      
    }
  return ret_s;
}

bool stringCompare_ext(const string& s1, const string& s2)
{
  int ret=s1.compare(s2);
  if(ret<=0)
    {
      return true;
    }
  else
    {
      return false;
    }
}
string flipStr(const string& _s1)
{
  char* temp=new char[_s1.size()+1];
  unsigned count=0;
  //char c;
  //unsigned indexFromEnd=_s1.size();
  for(unsigned i=_s1.size();i>0;i--)
    {
      //count=0;
      temp[count]=_s1.at(i-1);
      count++;
    }
  temp[count]='\0';
  string newStr(temp);
  delete [] temp;
  return newStr;
}

char DnaComplement(const char& c)
{
 
  char temp;
  //SequenceString temp(seq.GetName(), "");
  //string tempStrReturn("");
  //cout<<"temStr (seq get sequence):"<<tempStr<<endl;
  //cout<<"length:"<<tempStr.length()<<endl;
  
      switch(c)
	{
	case 'A':
	case 'a':
	  temp='T';
	  break;
	case 'T':
	case 't':
	  temp='A';
	  break;
	case 'U':
	case 'u':
	  temp='A';
	  break;
	case 'G':
	case 'g':
	  temp='C';
	  break;
	case 'C':
	case 'c':
	  temp='G';
	  break;
	case 'Y':
	case 'y':
	  temp='R';
	  break;
	case 'R':
	case 'r':
	  temp='Y';
	  break;
	case 'S':
	case 's':
	  temp='S';
	  break;
	case 'W':
	case 'w':
	  temp='W';
	  break;
	case 'K':
	case 'k':
	  temp='M';
	  break;
	case 'M':
	case 'm':
	  temp='K';
	  break;
	case 'B':
	case 'b':
	  temp='V';
	  break;
	case 'D':
	case 'd':
	  temp='H';
	  break;
	case 'H':
	case 'h':
	  temp='D';
	  break;
	case 'V':
	case 'v':
	  temp='B';
	  break;
	case 'N':
	case 'n':
	default:
	  temp='N';
	  break;
	}

      return temp;
}
bool is_nucleotide(const char &c)
{ 
      switch(c)
	{
	case 'A':
	case 'a':
	case 'T':
	case 't':
	case 'U':
	case 'u':
	case 'G':
	case 'g':
	case 'C':
	case 'c':
	case 'Y':
	case 'y':
	case 'R':
	case 'r':
	case 'S':
	case 's':
	case 'W':
	case 'w':
	case 'K':
	case 'k':
	case 'M':
	case 'm':
	case 'B':
	case 'b':
	case 'D':
	case 'd':
	case 'H':
	case 'h':
	case 'V':
	case 'v':
	case 'N':
	case 'n':
		return true;
		break;
	default:
		return false;
		break;
	}
    return false;
}
