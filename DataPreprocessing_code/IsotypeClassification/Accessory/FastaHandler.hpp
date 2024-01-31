#ifndef FASTAHANDLER_HPP
#define FASTAHANDLER_HPP
#include <string>
#include <vector>
#include <iostream>
#include "FileHandler.hpp"
#include "SequenceString.hpp"
#include "GzTools.hpp"

//here, we won't define a class, instead we simply write up a few functions to
//read, write, etc
using namespace std;

//return total number of sequeces read in.
//toUpper =true : convert the character to upper case. false: keep the original letter
//now this can handle automatically either fasta or gzip'ed fasta 
//
size_t ReadFasta(const string& _fname, vector<SequenceString>& _seqStrVec, bool toUpper=false);


void WriteFasta(const string& _fname, vector<SequenceString>& _seqStrVec, const unsigned int& _width=50, ios_base::openmode mode=ios_base::out);

//writing a text file, mainly for one vector.
//input: _fname, file name
//       _seqStrVec, vector of numbers to be written
//       _c, a char to delimite the columns
//       _width, how many columns to write
//       mode, to open the files, truc (new) or append
void WriteTextFile(const string& _fname, vector<unsigned int>& _seqStrVec, const char& c='\t', const unsigned int& _width=1, ios_base::openmode mode=ios_base::out);

//this is not for fasta files. for "regular" table format files
//Writing a text table file, with header or not
//input: _fname, file name
//       _seqStrVec, vector of vector of  numbers to be written
//       _c, a char to delimite the columns
//       _header, whether to write the header
//       mode, to open the files, truc (new) or append
//       _headerStr, the header names to be written.
void WriteTextTableFile(const string& _fname, vector<vector<double> >& _seqStrVec, const char& c='\t', const bool& _header=true, ios_base::openmode mode=ios_base::out,vector<string> _headerStr=vector<string>());

void WriteTextTableFile(const string& _fname, vector<vector<string> >& _seqStrVec, const char& c='\t', const bool& _header=true, ios_base::openmode mode=ios_base::out, vector<string> _headerStr=vector<string>());

//concatenate two fasta files,could be gzipped fasta!!
//The caller need to be repsponsible for making sure that the input files could only be
// FASTA or GZ_FASTA!!!
/*input :
 *	fn_r1 and fn_r2, string, name of the input read 1 and read 2 file.
 *	ft, FileType, cutomized enum file type. In this case, can only be FASTA or GZ_FASTA
 *	rc, bool to indicate whether to reverse complement the second read r2 sequences.
 *		true by default.
 *output: 
 *	return total number of records processed. if there are error, will return string::npos
 *
 */
size_t concatenateFasta(const string& fn_r1, const string& fn_r2, 
		FileType ft, string outFile_name, 
		ios_base::openmode mode=ios_base::out,
		bool rc =true);


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
bool get1FastaSeq(FILE* fb, SequenceString& ss, string& hangover, GZ_UTILITY& gu, const bool& gz=false );

/*Assuming the file has been open and writable. 
 * return bool to indicate wheter writting is good.
 */
bool writ1FastaSeq(ofstream& ofs, const SequenceString& ss,const unsigned& _width=100);

#endif
