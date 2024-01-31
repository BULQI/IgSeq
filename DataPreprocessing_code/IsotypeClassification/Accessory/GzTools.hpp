#ifndef GZTOOLS_HPP
#define GZTOOLS_HPP

//this is a C++ accessory library to take care of gzip file input and output
//we rely on the zlib to do uncompression

#include <zlib.h>
#include <string>
#include <vector> 

using namespace std;

//in this library, do we need to do objects??? not sure
//not now!!! at this point (12/9/2018) we simply define a 
//helper functions to take care of reading from gzipped files
//just reading. not doing writtings for now.

//one thing to note is that gzFile is a pointer 

//to read one line from an opened gzFile stream 
//--->NOT IMPLEMENTED YET!!!<------
gzFile getline(gzFile _if, string& _line, const char& _delimit='\n');

//-------------------------------------------
//NOTE: THIS ONE HAS NOT BEEN TESTED. COULD BE WRONG. NEED TO BE THOROULY TESTED 
//This one use the simply the zgets for reading.
//gzgets read in len -1 chars or when a newline char is read in
bool  getline_A(gzFile _if, string& _line); //, const char& _delimit='\n');
//gzFile getline_B(gzFile _if, string& _line, const char& _delimit='\n');

//for the following two, we simply calling the correct gzopen and gzclose
//function, we have them here to make holes for the future usage if we
//we want to do more
gzFile gzOpen(const string& _fname, const string& _mode="rb");

//return null upon return, in case the gzclose were called twice on
//the same file handler.
void gzClose(gzFile _if);
//---------------end of above section-----------

//----------------for getline_B-------------

/*this is the class/struct for doing the decompression.
 *it contains buffers/flag/line etc in order to compress and read
 *we make this a struct in order to read/use mutiple copies so that 
 *we could read a few files at the same time.
 *
 */
//one thing to note is that gzFile is a pointer 
#define BUF_SIZE 0x1000
//#define BUF_SIZE 0x100
//#define BUF_SIZE 0x50
//#define DEBUG

#define CHUNK 0x50000
//#define CHUNK 0x200
//#define CHUNK 0x20
#define OUT_CHUNK CHUNK*1
 /*
typedef struct gz_utility2{
	//this is a stream for gzreading and decompression
	//z_stream strm ;//= {0}; //initialize all members to be zero
	unsigned char gzip_in[CHUNK];
	//unsigned char gzip_out[OUT_CHUNK];
	
	//
	char* first_line;//=(char*)&gzip_out[0];
	char* current_line;//=first_line;
	char* next_line;//=first_line;

	unsigned current_line_size;//=0;
	unsigned next_line_size;//=0;

	char hangover[OUT_CHUNK+2];
	bool z_stream_init;//=false;
} GZ_UTILITY2; */
 
typedef struct gz_utility{
	//this is a stream for gzreading and decompression
	z_stream strm ;//= {0}; //initialize all members to be zero
	unsigned char gzip_in[CHUNK];
	unsigned char gzip_out[OUT_CHUNK];
	
	//
	char* first_line;//=(char*)&gzip_out[0];
	char* current_line;//=first_line;
	char* next_line;//=first_line;

	unsigned current_line_size;//=0;
	unsigned next_line_size;//=0;

	char hangover[OUT_CHUNK+2];
	bool z_stream_init;//=false;
} GZ_UTILITY;
//extern GZ_UTILITY _gu;//this is a default one that can be used 
//assuming the file has been opened 
bool getline_B(FILE* _f, string& l, GZ_UTILITY& gu);

FILE* gzOpen_B(const string& _fname,  GZ_UTILITY& gu,const string & _mode="rb");

void gzClose_B(FILE* _f, GZ_UTILITY& gu);

int init_gzip_stream( GZ_UTILITY& gu,bool full=false/*FILE* file,char* out*/);
bool inflate_gzip( GZ_UTILITY& gu, const size_t& bytes_read, const size_t& bytes_avail, int& code);

void init_gu(GZ_UTILITY& gu);
//void init_gu2(GZ_UTILITY2& gu);
//this is the helper function that is called to print the stream contents
//should be called printZstream.
//void debugZstream(const string& s="");


//a helper function to find a char in a char buffer, note here we don't assume the char buffer is 
//a null-terminated char c string. so we have to require a size of char buffer.
//therefore we stop either null or size whichever comes first.
//so this is why we don't use the build in strstr, since it assume null terminator.
//input:
//	source char* a pointer to the char array that needs to be searched for a char c_str
//	size size_t the size of the source array. or the size to be searched.
//	c char the character to be search for in the array.
//output:
//	return the pointer to the found char in the array. nullptr if otherwise.
//
char* strcnstr(char* source, size_t size, const char& c);
#endif

