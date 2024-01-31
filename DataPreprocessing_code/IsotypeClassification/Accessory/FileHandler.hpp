#ifndef FILEHANDLER_HPP
#define FILEHANDLER_HPP
//in here we define some general functions to take care of files
#include <vector>
#include <string>
#include <iostream>
#include <sys/stat.h>
#include "SequenceString.hpp"


using namespace std;

enum FileType {GZ, FASTA, FASTQ, TXT, GZ_FASTA, GZ_FASTQ, GZ_TXT, DIR, UNKNOWN, NOT_EXIST};

//get file type, note FileType is user-defined enum.
//we will do two ways, either by file name, simply check the file name 
//	or deep mode, to read the file and the check the first two bytes/chars 
//	gzip file first byte is 31 and second 139
//	fasta '>'
//	fastq '@'
//input :
//	fname string name of the file 
//	deep bool to read first 2 byte to tell what file this is.
//return: 
//	FILETYPE. this is the userdefine and file type.
FileType getFileType(const string& fname, bool deep=false);

bool is_file(const char* path) ;

bool is_dir(const char* path) ;

bool exist(const char* path);

//a file handler to read file into a fasta vector vector
//we try to detect the following thing in this function 
// 1, file exist?
// 2, is it a file or directory
// 3, is it a fasta, fastq or gziped fasta, gziped fastq. No other type supported so far
// 
// We will return a vector holding the squenences of the file (SquenceStrings)
// we will return the total number of sequences read in. If none or error, we will return 
// string::npos.
//----updated 8/21/19
// expanded this to reading possibly fasta, fastq, fastq gzipped.
// we will return separately the sequence string and quality string in order to be compatible 
// with the previous version of the function.
size_t readFile2SeqStrVector(const string& _fname, vector<SequenceString>& _vec, vector<string>* _vec_Q=NULL);

/* read two files (pair end reads) and then cancatenate them 
 * assuming they are r1 and r2, we also need to revcomp the second 
 * reads.
 * input:
 * 		r1 string file name of read 1
 *		r2 string file name of read 2
 * output: 
 * 		return the total number of reads processed.
*/
size_t concatnateSeqFiles(const string& fn_r1, const string& fn_r2);

FileType check_gzFileType(const string& fname);

bool is_fastq(const string& fname);
bool is_fasta(const string& fname);
bool is_text(const string& fname);

FileType getFileType_deep(const string& fname);
FileType getFileType_byName(const string& fname);

//return the basename of a file path.
string basename (const std::string& str);

#endif