#ifndef FASTQHANDLER_HPP
#define FASTQHANDLER_HPP
#include <string>
#include <vector>
#include <iostream>
#include <string.h>
//#include <stdlib>
//#include <stdio>
#include "SequenceString.hpp"
#include "GzTools.hpp"
#include "FASTQ.hpp"
#include "FileHandler.hpp"
#include "../SequenceHandlerCommon.hpp"

//here, we won't define a class, instead we simply write up a few functions to
//read, write, etc
using namespace std;

//return total number of sequeces read in.
//toUpper =true : convert the character to upper case. false: keep the original letter
//
//the file could be gz'ed or regular fastq. then the sequence data are read into a vector of 
//fastq objects. 
//return string::npos upon error 
size_t ReadFastq(const string& _fname, vector<Fastq>& _seqStrVec, bool toUpper=false);

//read fastq files into two vectors, one holding the sequences, the other holding the
//quality string.
//
size_t ReadFastq(const string& _fname, vector<SequenceString>& _seqStrVec, vector<string>& _vecQ, bool toUpper=false);


void WriteFastq(const string& _fname, vector<Fastq>& _seqStrVec,  ios_base::openmode mode=ios_base::out);
//updated 8/21/2019, add this to write the fastq based on the separated sstring and quality.
void WriteFastq(const string& _fname, vector<SequenceString>& _seqStrVec, 
	vector<string>& _qVec,
	ios_base::openmode mode=ios_base::out);

//read a sequence from the file stream. 
//in here we read ONE record of fasta 
//the caller needs to open/close  the file
/*input: 
 *	fb FILE stream pointer, assuming it is ok. but we need to check for the end of file status.
 *	ss SequenceString reference, used to return the read record.
 *	qs, string quality string as a return  
 *		
 *  gu, GZ_UTILITY, used by the gz tools for reading the compressed files. 
 *		will be passed between calls. might not be needed. 
 *	gz bool, used to indicating whether this is a compressed stream. 
 *output:
 *	return the bool indcating whether the reading COULD BE CONTINUE.
 *		true: meaning we can go on to read and current record is good. 
 *		false when we reach the end or there is an error (either stream reading error or parsing error (no '>', etc).
 *		the outer caller need to check whether it is a eof or error. 
 *  NOTE:!!!    Most importantly, it could happen that the reading can not go on (false), but the current 
 *		read is still valid. Just check to see whether it is a good read. if the read is no good
 *		due to the error, the rescord is set to be empty.!!!!
 */
bool get1FastqSeq(FILE* fb, SequenceString& ss, string& qs, GZ_UTILITY& gu, const bool& gz=false );

//concatenate two fasta files,could be gzipped fasta!!
//The caller need to be repsponsible for making sure that the input files could only be
// FASTQ or GZ_FASTQ!!!
/*input :
 *	fn_r1 and fn_r2, string, name of the input read 1 and read 2 file.
 *	ft, FileType, cutomized enum file type. In this case, can only be FASTA or GZ_FASTA
 *output: 
 *	return total number of records processed. if there are error, will return string::npos
 *
 */
size_t concatenateFastq(const string& fn_r1, const string& fn_r2, 
		FileType ft,string outFile_name,ios_base::openmode mode=ios_base::out,
		bool rc =true);
		
bool writ1FastqSeq(ofstream& ofs, const SequenceString& ss, const string& qs);

#endif
