#include "GzTools.hpp"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "FileHandler.hpp"

//for printing out the content of the char array.
//static void debugChars(const char* s, size_t t, const string& n="");
//not implemented yet.
//to read one line from an opened gzFile stream 
gzFile getline(gzFile _if, string& _line, const char& _delimit)
{
	return 0;
}

//In there, just for fun, we will create two different ways of get line
//
/*typeA, we are doing the simple one with the gzgets
* NOTE: THIS ONE HAS BEEN TESTED. COULD BE WRONG. NEED TO BE THOROULY TESTED 
*Input: _if, an opened gzFile pointer, 
*			_line, a line buf to read things in 
*Output: bool TRUE if the reading is good. 			 
*/
bool getline_A(gzFile _if, string& _line )
							//, const char& _delimit)
{
	unsigned v_size=BUF_SIZE;
	char* v=new char [ v_size ];
    unsigned pos = 0;
    for ( ;; ) {
        if ( gzgets( _if, &v[ pos ], v_size - pos ) == 0 ) {
            // end-of-file or error
            int err;
            const char *msg = gzerror( _if, &err );
            if ( err != Z_OK ) {
                // handle error
				cout<<"ERROR: "<<msg<<endl;
            }
			//clean up to avoid buffer leaking
			delete [] v;
            return false;
        }
		//gzgets read in len -1 chars or when a newline char is read in
        unsigned read = strlen( &v[ pos ] );
		
		//we first check for the case where read in something >=1
		//instead of zero read 
        if ( read >= 1 && v[ pos + read - 1 ] == '\n' ) {
            if ( read >= 2 && v[ pos + read - 2 ] == '\r' ) {
                pos = pos + read - 2;
            } else {
                pos = pos + read - 1;
            }
            break;
        }
		//we did not see new line so far, but we do stop reading in more 
		//we are reaching the end of the buffer
		// case 1: read=0??? read in zero thing, ? how?
		//this is only possible when this is the last 
		//line and there is no new line in the end. an empty ending line 
		// case 2: pos+read <v_size -1, we read in something, but
		//not fill the whole buf and we did not see a new line.
		//we are done with reading, but no new line.  we are the line
		// finishing but has something.
        if ( read == 0 || pos + read < v_size - 1 ) {
            pos = read + pos;
            break;
        }
			
		//if we are here, we are reading new line and did not break,
		//so 
		//this is the case, where we need keep reading  
        pos = v_size - 1;
        v_size=v_size * 2;
		char* w=new char[v_size];
		strcpy(w,v);
		delete [] v;
		v=w;
    }
    v[ pos ]='\0';
    _line.assign(v);
	return true;
}

///-----------the following section contains code for getlinne B function
#define windowBits 15
#define ENABLE_ZLIB_GZIP 16
//define a default gzip untility variable to be used
//GZ_UTILITY _gu {0};


/*
//this is a stream for gzreading and decompression
static z_stream strm = {0}; //initialize all members to be zero
static unsigned char gzip_in[CHUNK];
static unsigned char gzip_out[OUT_CHUNK];
*/
//static 
int init_gzip_stream(GZ_UTILITY& gu, bool full/*FILE* file,char* out*/ 
				 
			)
{// unsigned     
        
		//need to check for errors
        int x=Z_OK;
		if(full)
		{//The application must initialize zalloc, zfree and opaque before calling the init function. 
			gu.strm.zalloc = Z_NULL;
			gu.strm.zfree = Z_NULL;
			gu.strm.opaque = Z_NULL;
			gu.strm.next_in = Z_NULL; //bytes holding the input data 
			gu.strm.avail_in = 0;
			x=inflateInit2 (& gu.strm,  windowBits| ENABLE_ZLIB_GZIP);
			//x=inflateInit (& strm);//,  windowBits | ENABLE_ZLIB_GZIP);
			if(x!=Z_OK)
			{
				cout<<"******ERROR: in initialize the strm buffer..........."<<endl;
				//exit(-1); we don't have to quit here, we return the signal to the caller and let the caller to decide and handle.
			}
			else
			{
				//cout<<"****good at initializing buffer............."<<endl;
			}
			gu.strm.next_in = gu.gzip_in; //bytes holding the input data 
			//strm.avail_in = 0;
			gu.strm.next_out = gu.gzip_out; //bytes holding the output data
			gu.strm.avail_out=OUT_CHUNK;
		}
		else
		{ //The application must update next_in and avail_in when avail_in has dropped to zero. 
		  //It must update next_out and avail_out when avail_out has dropped to zero. 
			if(gu.strm.avail_in==0)
			{
				gu.strm.next_in=gu.gzip_in;
			}
			if(gu.strm.avail_out==0)
			{
				gu.strm.next_out=gu.gzip_out;
				gu.strm.avail_out=OUT_CHUNK;
			}
		}
		//All other fields are set by the compression library and must not be updated by the application.
		
    return x; //strm;
}
/*
void debugZstream(const string& s)
{
	
		cout<<"\n-*&Debugging strm*--:"<<s<<endl;
		cout<<"\t^next_in:"<<strm.next_in<<endl;
		cout<<"\t^avail_in:"<<strm.avail_in<<endl;
		cout<<"\t*next_out:"<<strm.next_out<<endl;
		cout<<"\t*avail_out:"<<strm.avail_out<<endl;
		
		cout<<"\t&gzip_out initialized lenth:"<<OUT_CHUNK<<endl;
		cout<<"\t^gzip_out detected lenth:"<<strlen((char*)gzip_out)<<endl;
		cout<<"\t^gzip_out string: "<<gzip_out<<endl;
		cout<<"\t~gzip_in first char: ID1--"<<(int)(gzip_in[0])<<";ID2--"<<(int)(gzip_in[1])
			<<";CM--"<<(int)(gzip_in[2])<<endl;
		cout<<"\t*pointer addresses---gzip_in:strm.next_in ("<< static_cast<const void *>(gzip_in)<<":"<<static_cast<const void *>(strm.next_in)<<")"<<endl;
		cout<<"\t* gzip_out:strm.next_out ("<<static_cast<const void *>(gzip_out)<<":"<<static_cast<const void *>(strm.next_out)<<")"<<endl;		
		//cout<<"\tsizeof(char):"<<sizeof(char)<<"sizeof(unsigned char):"<<sizeof(unsigned char)<<endl;
		cout<<"\t*length---next_out: "<<strlen((char*)strm.next_out)<<";gzip_in:"<<strlen((char*)gzip_out)<<endl;
		cout<<"\t*end at null????gzip_in of "<<strlen((char*)gzip_out )<<":"<<(unsigned)gzip_out[strlen((char*)gzip_out)]<<endl<<endl;
		
		//debugging the bugger
		/// *cout<<"\t**Debugging:gzip_out buf diplay:"<<endl;
		cout<<"\t**";
		for(unsigned i=0;i<OUT_CHUNK;i++)
		{
			if(gzip_out[i]!='\0')
				cout<<"["<<i<<"]:"<<(char)gzip_out[i]<<";";
			else 
				cout<<"["<<i<<"]:"<<"\\0"<<";";
		}
		cout<<endl;
		/// * /
		debugChars((char*) gzip_out, OUT_CHUNK, "gzip_out (strm output)");
}
*/
/* void debugChars(const char* s, size_t t, const string& n)
{
	cout<<"**Debugging:\""<<n<<"\" buf diplay:"<<endl;
		cout<<"\t**";
		
		for(unsigned i=0;i<t;i++)
		{
			cout<<"["<<i<<"]:";
			switch (s[i])
			{
			case '\0':
				cout<<"\\0"<<";";
				
				break;
			case '\n':
				cout<<"\\n"<<";";
				break;
			default: 
				cout<<(char)s[i]<<";";
				break;
			}
		}
		cout<<endl;

}
*/
//static
 bool inflate_gzip( GZ_UTILITY& gu, const size_t& bytes_read, const size_t& bytes_avail, int& code )
 {
			z_stream& strm =gu.strm;
			strm.avail_in = (int)bytes_read;
            
            strm.avail_out = (int)bytes_avail;
			unsigned avail_out_before=strm.avail_out;
            int rt=inflate (& strm, Z_NO_FLUSH);//with Z_NO_FLUSH, we mean to continue if there are more input/more space for output
  //            printf ("string %s, length %i",gzip_out, strm.avail_out);
			unsigned avail_out_after=strm.avail_out;
			bool ok=true;
			int y=0;
			switch(rt)
			{
			case Z_OK: //if some progress has been made (more input processed or more output produced), 
				ok=true;
#ifdef DEBUG
				cout<<"z_stream OK:"<<endl;
				//cout<<"\t"<<strm.msg<<endl;
#endif
				break;
			case Z_STREAM_END:  //if the end of the compressed data has been reached and all uncompressed output has been produced, 
				ok=true;
				//strm.avail_in=0;
#ifdef DEBUG
				cout<<"z_stream reaching end:"<<endl;
				//cout<<"\t"<<strm.msg<<endl;
#endif
				y=inflateReset(&strm);
				if(y!=Z_OK)
				{
					cout<<"ERRRROR::::::::::after z_stream_end and reset"<<endl;
				}
				break;
			case Z_NEED_DICT: // if a preset dictionary is needed at this point, 
				ok=false;   //not what does it mean!!
				cout<<"z_stream (need) dictionary error upon inflating:"<<endl;
				//cout<<"\t"<<strm.msg<<endl;
				break;
			case Z_DATA_ERROR: // if the input data was corrupted (input stream not conforming to the zlib format or incorrect check value, in which case strm->msg points to a string with a more specific error), If Z_DATA_ERROR is returned, the application may then call inflateSync() to look for a good compression block if a partial recovery of the data is to be attempted.
				ok=false;
				cout<<"z_stream data corruption error upon inflating:"<<endl;
				cout<<"\t"<<strm.msg<<endl;
				break;
			case Z_STREAM_ERROR: // if the stream structure was inconsistent (for example next_in or next_out was Z_NULL, or the state was inadvertently written over by the application), 
				ok=false;
				cout<<"z_stream stream error upon inflating (possibly parameter setting error):"<<endl;
				cout<<"\t"<<strm.msg<<endl;
				break;
			case Z_MEM_ERROR: // if there was not enough memory, 
				ok=false;
				cout<<"z_stream (not enough) memory error upon inflating:"<<endl;
				cout<<"\t"<<strm.msg<<endl;
				break;
			case Z_BUF_ERROR: //if no progress was possible or if there was not enough room in the output buffer when Z_FINISH is used. Note that Z_BUF_ERROR is not fatal, and inflate() can be called again with more input and more output space to continue decompressing. 
				ok=true;
				cout<<"z_stream buffer error upon inflating:"<<endl;
				//cout<<"\t"<<strm.msg<<endl;
				break;
			default:
				ok=false;
				cout<<"Unknown error, something is wrong"<<endl;
				cout<<"\t"<<strm.msg<<endl;
				break;
			}
					//<--------------a bug in the original code ??? be careful
			int x=Z_OK;
#ifdef DEBUG
		debugZstream("inside the inflate() before reset");
#endif
			//checking for another specical case when we reading the headers we report Z_OK, but no data compressing, 
			//this is a special case. In this case, Z_OK is reported and, but strm.avail did not change, keep at OUT_CHUNK;
			//and this is not like Z_BUF_ERROR, in which strm.avail_out did not and kept at OUT_CHUNK (although it is possilbe
			//to be other values too
			
			code=rt;
			if(rt==Z_OK&&avail_out_before==avail_out_after&&avail_out_after==OUT_CHUNK)
			{
				code=100; //<---this is my own defined INFLATE_HEADER_ONLY 
			}
			if(ok)		
			 x=init_gzip_stream(gu,false);  //<--we are resetting the parameters of the strm (checking for avail_in and avail_out for conditions)
			
#ifdef DEBUG
		debugZstream("inside the inflate() at the end");
#endif
		if(x!=Z_OK)
			return false;
    return ok;// all OK
}
/*
static char* first_line=(char*)&gzip_out[0];
static char* current_line=first_line;
static char* next_line=first_line;

static unsigned current_line_size=0;
static unsigned next_line_size=0;

static char hangover[OUT_CHUNK+2];
static bool z_stream_init=false;
*/
//In there, just for fun, we will create two different ways of get line
//but will only be called (used) by getline for better typedef, these 
//two will be moved to cpp, so not be called directly by other files

//assuming the file has been opened 
bool getline_B(FILE* _f, string& l, GZ_UTILITY& gu)
			//, const char& _delimit)
{
	//assume at this point, the file has been opened correctly.	
	bool ok=true;
	unsigned size_line=BUF_SIZE;
	gu.current_line=gu.next_line;
	gu.current_line_size=gu.next_line_size;
	char* line=new char [size_line];
	line[0]=0; //set it to empty
	unsigned line_len=0;
	unsigned hangover_len=0;
	size_t bytes_read=0;
	//size_t bytes_avail=0;
	//case 1: current line is not initialized 
	//case 2: current line is empty
	//case 3: reading outside of the buffer now.
	//in all these cases, we need to read the file buffer and inflate
	do{	
		//cout<<"reading......."<<endl;
#ifdef DEBUG	
		cout<<"reading......."<<endl;
#endif
		if(!gu.current_line || gu.current_line_size==0 || ((unsigned char*) gu.next_line) >= gu.gzip_out + OUT_CHUNK) //this is necessary, since in here we check for whether we have done with reading a line in the available buffer. the enclosed section check whether we have the read in buffer processed/inflated and whether to read in more input from the file buffer.
		{   //if we are doing a new call of "getline", we are checking actually the next_line in the previous call.
			// the third condition in above will not used if we have done a correct job previously, just keep it here for safety.
			//if we are in the current call, but come here because we are looping from previous unfound "new line" case, we are 
			//be jumping using current_line_size==0 condition (2nd condition), we will not use the third condition either.
			
			//if we are there we need either to read file or inflate more. in any case we will need to rest to here ??using th estrm.next_out 
			gu.current_line= (char*)(gu.strm.next_out);
			gu.next_line=gu.current_line; //reset this in case next_line is out of boundary of gizp_out and no new line found, this might create a new read.(this is not necessary, since we set the current line to zero, so we should be fine)
			bytes_read=gu.strm.avail_in;//<---
			//bytes_avail=strm.avail_out; //<---
			
			if(gu.strm.avail_in==0&&gu.strm.avail_out!=0) //strm.avail_in !=0 means we have more input to decompress and strm.avail_out means we have more output to write 
			{//we are here means we need to do more reading to get data 
			//we here check whether we have more file to read before reading the file .
				if (feof (_f)) { //we are done 
					//cout<<"reaching the end of the input file"<<endl;
					int err=inflateEnd (& (gu.strm));
					//check for error 
					if(err ==Z_STREAM_ERROR)
						cout<<"ERROR: when closing the stream, Z_STREAM_ERROR is returned!!"<<endl;
					l.assign(line);//copy over everything (in case there is something left) and return
					return false;
				}
#ifdef DEBUG
				cout<<"reading the file ....."<<endl;
				//before reading we need to reset everything, we assume in here we have copy over every thing to correct buffer
				//mainly hangover 
				//we only need to reset the pointers and size
				//int x=init_gzip_stream(false) ; //not a full one 
				//if(x!=Z_OK)
				//	return false;
//#ifdef DEBUG
		debugZstream("before reading the bytes");
				cout<<"-----bytes_read:"<<bytes_read<<endl;
#endif
				
				bytes_read = fread (gu.gzip_in, sizeof (char), CHUNK, _f);
#ifdef DEBUG
		cout<<"-----bytes_read:"<<bytes_read<<endl;
		debugZstream("after reading the bytes");
#endif
			
			}			
			//in any case, either we have read more or not, we need to inflate the buffer 
//			cout<<"before inflating....."<<endl;
			//here now we have other things to do, figure out the length of current 
			gu.current_line_size=gu.strm.avail_out; //remember it before inflating, as holder 
			int inflate_code =Z_OK;
			ok=inflate_gzip(gu, bytes_read, gu.strm.avail_out, inflate_code);
			//if reading and inflation are OK.
			if(!ok){//so far, ok means everything is fine. otherwise (ok==false) we got an error, need to stop.
				/*//here we are done with reading the file by reaching the end, 
				if((unsigned char* )current_line>=gzip_out+bytes_read) //means, we are out of uncompressed bytes and ok, reach end of file
					return false ;*/
				cout<<"someting is wrong upon inflating. we are jumping out......."<<endl;
				l.assign(line);
				return false;
			}
			//calculate the size of current, it is a little complicated, note now the current_line_size holding 
			//the strm.avail_out before reading, which is the possible len to before inflating meaning the max-value for current_line_size
			if(gu.current_line_size>=gu.strm.avail_out)  //this is no wrapping up (because resetting) and we still have some free bytes
			{//could be no writing of the output due to some reason (decompressing and not more buffer)
				if(gu.strm.avail_out==OUT_CHUNK&&gu.current_line_size==OUT_CHUNK)//this is a very special case
				{
					//we need to check for z_stream comppressing eror , if not then we are good
					if(inflate_code !=Z_BUF_ERROR&&inflate_code!=100) //the only possible cases are Z_BUF_ERROR, Z_STREAM_END and Z_OK  if ok is true 
					{
						//note 100 is my own defined code, which is reading the header only
						gu.current_line_size=OUT_CHUNK;
					}
					else  //no read, so 
					{
						gu.current_line_size=0;
					}
				}
				else  //case where we read something or did not read anything.
				{
					gu.current_line_size-=gu.strm.avail_out;//we got the bytes
				}
			}
			else  //in here means we wrapping out due to the fact that we used up the buffer in inflating.
			{//and in this case, the only possible value is OUT_CHUNK.
				if(gu.strm.avail_out!=OUT_CHUNK)
				{	cout<<"ERROR: something is wrong, strm avail is not what we expected, please check!!!"<<endl;
					return false;
				}
									
				gu.current_line_size=gu.current_line_size-0;//0 is the value of strm.avail_out before updating the params in strm
			}
#ifdef DEBUG
			cout<<"&&&&size debug, current_line_size:"<<current_line_size<<endl;
			cout<<"after inflate .."<<endl;
			
			cout<<"Now start doing the line parsing........."<<endl;
#endif
			//copy over anything in hangover, because we might need to using it"
			if(hangover_len>0)
			{
				//check to see whether we have enough othwise to increase
				while((line_len+hangover_len)>size_line)
				{
					size_line*=2;
					char* temp=line;
					
					line=new char[size_line];
					strncpy(line,temp, size_line/2);
					delete [] temp;
				}
			//strcpy(line+line_len,hangover);
#ifdef DEBUG
				debugChars(line, size_line, "line (buf display) before copy handover in case ******");
				cout<<"line_len:"<<line_len<<"---hangover_len:"<<hangover_len<<endl;
#endif
				strncpy(line+line_len,gu.hangover, hangover_len); //<- this is just in case there are leftovers in the hangover region. this could because we did not find new line in the last round of reading and inflating 
				line_len += hangover_len;//<-------- size_line and line_len??????
#ifdef DEBUG				
				debugChars(line, size_line, "line (buf display) after copy handover in case ******");
#endif
				gu.hangover[0]=0;
				hangover_len=0;
				
			}
			
		}
		//cout<<"inflating....."<<endl;
#ifdef DEBUG
		cout<<"Before searching for \\n, current_line_size:"<<gu.current_line_size<<endl;
		cout<<"\t*******current:"<<gu.current_line<<endl;
#endif 		
		gu.next_line=strcnstr(gu.current_line,gu.current_line_size,'\n');//current line pointing to the inflated z_stream output region , and now we start parsing (meaning trying to find a newline)
		if(gu.next_line){//find a newline, then we need to point the next_line to the correct region 
#ifdef DEBUG
			cout<<"\tfound a new line"<<endl;
#endif			
			gu.next_line[0]=0; //finish/replace the end of  current line with '\0'
			gu.next_line++;  //point to the next line 
			while(gu.next_line-gu.current_line+line_len>size_line)
			{
				//size_line*=2;
				//delete [] line;
				//line=new char[size_line];
				size_line*=2;
				char* temp=line;
					
				line=new char[size_line];
				strncpy(line,temp, size_line/2);
				delete [] temp;
			}
			gu.next_line_size=gu.current_line_size-(gu.next_line-gu.current_line);
			gu.current_line_size=gu.next_line-gu.current_line;
			//next_line_size=(char*)(strm._)-next_line;
			
#ifdef DEBUG			
			debugChars((char*)gzip_out, OUT_CHUNK, "gzip_out(zstrm output)");
			cout<<"--before copying, current_line:"<<current_line<<"****current_line_size:"<<current_line_size<<endl;
			debugChars(line, size_line, "line (final output");
			cout<<"pointer address line:"<<static_cast<const void *>(line)<<"; plus line size("<<line_len<<"):"<<static_cast<const void *>(line+line_len)<<endl;
			cout<<"next line size :"<<next_line_size<<endl;
#endif			
			strncpy(line+line_len,gu.current_line, gu.current_line_size); //copy to the hangover/buffer zone, increase it in case this is too big 
#ifdef DEBUG			
			cout<<"--after copying, current_line:"<<gu.current_line<<"****current_line_size:"<<gu.current_line_size<<endl;
			debugChars(line, size_line, "line (final output)");

			cout<<"\t**the line is "<<line<<endl;
	
			//cout<<"\tgzout:"<<gzip_out<<endl;
#endif
			//cleanup and ready to jump out 
			gu.hangover[0]=0;//set the buffer to empty
			hangover_len=0;
			//cout<<"line_len before:"<<line_len;
			line_len+=(gu.next_line-gu.current_line);
#ifdef DEBUG
			cout<<";line_len after:"<<line_len<<endl;
			cout<<"done for the all loop rounds........."<<gu.current_line<<endl;
#endif
			break;
		}else{ //could not find a new line in this case we need to read more bytes or inflate more
			//here, we don't need to check for the length, we set up previously the current_line buffer, gzip.out to be the same size as the hangover. we are safe
#ifdef DEBUG
			//cout<<"\t current_line_size:"<<current_line_size<<endl;
			//cout<<"\t hangover_len:"<<hangover_len<<endl;
#endif			
			//please do make sure you set the hangover(dest buffer is bigger enough to hold the copying.) otherwise
			//anything could happen.!!!!! Yes, we do , in this case, hangover is 1 more than gzip_out(OUT_CHUNK)
			strncpy(gu.hangover,gu.current_line, gu.current_line_size);//before start reading, we need to copy over to hangover buffer
			//line[0]=0;// skip that one!!
#ifdef DEBUG
			cout<<"\t hangover_len:"<<hangover_len<<endl;
			cout<<"\t current_line_size:"<<gu.current_line_size<<endl;
#endif
			hangover_len = gu.current_line_size;
				//strm.avail_out/sizeof(char)-(current_line-(char*)strm.next_out);
#ifdef DEBUG
			cout<<"##within the hangover region, hangover_len:"<<hangover_len<<"; hangover[len]:"<<
					(unsigned)gu.hangover[hangover_len]<<endl;
			cout<<"\t current_line_size:"<<gu.current_line_size<<endl;
#endif			
			//dump all the data in the current_line, since we need to either read more or decompress more 
			//current_line[0]=0; //<--???by doing this we "destroyed" the gzip_out buffer holding the usedcurrentline data.
						//the one before strm.next_in
					
			gu.current_line_size=0;//this is necessary
			
			//in this case, we will not use next line, since there is no '\n' found. we will go back the loop read 
			//more data and will not go out of loop for a new call.
		}
#ifdef DEBUG		
		cout<<"done for one round......"<<hangover<<endl;
#endif
	}while(true);
	l.assign(line); 
	delete [] line ;
	return true;
}
FILE* gzOpen_B(const string& _fname,  GZ_UTILITY& gu, const string & _mode)
{
	//here we now simly call the zlib open 
	if(_fname.length()==0)
	{
		cout<<"ERROR: can not open gz file. the file name has not been specified correctly"<<endl;
		return NULL;
	}
	//cout<<"start opeing....."<<_fname<<endl;
	//open the file 
	FILE* f=fopen(_fname.c_str(), _mode.c_str());
	//cout<<"here....."<<endl;
	if(!f)
	{
		////cout<<"inside........"<<endl;
		cout<<"ERROR: file open error !!"<<endl;
		perror("\terr msg:");
		//clearerr(f);
		//fclose(f);
		return 0;
	}
	init_gu(gu);
	gu.z_stream_init=false;
	//initiate/start the zstrem
	if(!gu.z_stream_init)
	{	
		gu.z_stream_init=true;
		init_gzip_stream( gu,true/*_f,&line[0]*/);
	}
	
	
	return f;
}
void init_gu(GZ_UTILITY& gu)
{
	//cout<<"inside inti gu ...."<<endl;
	gu.strm = {0}; //initialize all members to be zero
	//unsigned char gzip_in[CHUNK];
	//unsigned char gzip_out[OUT_CHUNK];
	
	//
	gu.first_line=(char*)(&(gu.gzip_out[0]));
	gu.current_line=gu.first_line;
	gu.next_line=gu.first_line;

	gu.current_line_size=0;
	gu.next_line_size=0;

	//char hangover[OUT_CHUNK+2];
	gu.z_stream_init=false;
}
/*void init_gu2(GZ_UTILITY2& gu)
{
	cout<<"inside inti gu ...."<<endl;
	gu.strm = {0}; //initialize all members to be zero
	//unsigned char gzip_in[CHUNK];
	//unsigned char gzip_out[OUT_CHUNK];
	
	//
	gu.first_line=(char*)(&(gu.gzip_out[0]));
	gu.current_line=gu.first_line;
	gu.next_line=gu.first_line;

	gu.current_line_size=0;
	gu.next_line_size=0;

	//char hangover[OUT_CHUNK+2];
	gu.z_stream_init=false;
}*/

void gzClose_B(FILE* _f, GZ_UTILITY& gu)
{
	//reset everything.
	gu.z_stream_init=false;
	gu.strm={0};
	gu.first_line=(char*)&(gu.gzip_out[0]);
	gu.current_line=gu.first_line;
	gu.next_line=gu.first_line;

	gu.current_line_size=0;
	gu.next_line_size=0;

	//cout<<"close the file stream.........."<<endl;
	//we leave the buffers alone, as long as we set the flags correctly 
	//do I need to close or de-initialize the z_stream??
	fclose(_f);
}
//----------------end of the getline_B function section----------


//for the following two, we simply calling the correct gzopen and gzclose
//function, we have them here to make holes for the future usage if we
//we want to do more

/*Input:   _fname, input file name string
*		_mode, a string indicating the mode to open the file.
*					Note, we are doing the binary file, not text ifle. so 
*					the mode can only be rb, rw..... check fopen c++ function for detail.
*Output: return a point gzifle is a pointer. 
*/
gzFile gzOpen(const string& _fname, const string & _mode)
{
	//here we now simly call the zlib open 
	if(_fname.length()==0)
	{
		cout<<"ERROR: can not open gz file. the file name has not been specified correctly"<<endl;
		return NULL;
	}
	//check for file existence.
	if(!exist(_fname.c_str()))
	{
		return NULL;
	}
	//now call zlib
	gzFile f= gzopen(_fname.c_str(), _mode.c_str() );
	if(f==NULL)
	{
		int eno;
		const char* err=gzerror(f,&eno); 
		cout<<"ERROR: gz file open error(error code:"<<eno<<")!!!"<<endl;
		cout<<"\terror msg:"<<err<<endl;
	}
	return f;
}

//return null upon return, in case the gzclose were called twice on
//the same file handler.
void  gzClose(gzFile _if)
{
	int state =gzclose( _if);
	//check
	switch (state) 
	{
	case Z_OK:
		cout<<"gz file closed successfully"<<endl;
		break;
	case Z_STREAM_ERROR:
		cout<<"gz file closed unsuccessfully: file not valid"<<endl;
		break;
	case Z_ERRNO:
		cout<<"gz file closed unsuccessfully: file operation error"<<endl;
		break;
	case Z_MEM_ERROR:
		cout<<"gz file closed unsuccessfully: out of memory"<<endl;
		break;
	case Z_BUF_ERROR:
		cout<<"gz file closed unsuccessfully: ended in the middle of a gzip stream"<<endl;
		break;
	default:
		cout<<"gz file closed unsuccessfully: unknown error"<<endl;
		return ;
	}
	
	return ;
}

//a helper function to find a char in a char buffer, note here we don't assume the char buffer is 
//a null-terminated char c string. so we have to require a size of char buffer.
//therefore we stop either null or size whichever comes first.
//so this is why we don't use the build in strstr, since it assume null terminator.

char* strcnstr(char* source, size_t size, const char& c)
{
	char* pos=NULL;
	
	//for loop
	for(unsigned i=0;i<size;i++)
	{
		if(c==source[i])
		{
			pos=&(source[i]);
			break;
		}
		if(source[i]=='\0')
		{
			pos=NULL; //not necessary, should be NULL anyway.
			break;
		}
	}

	return pos;
}	


