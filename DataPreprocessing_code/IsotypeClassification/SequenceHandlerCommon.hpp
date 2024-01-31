#ifndef SEQUENCEHANDLERCOMMON_HPP
#define SEQUENCEHANDLERCOMMON_HPP
#include "Accessory/SequenceString.hpp"
#include "score.hpp"

//we don't do in-place modification, but instead make a new 
//copy and return it. 
//
SequenceString ReverseComplement(SequenceString& seq);

//compare two strings character by character, return # of chars that are different. if two strings are of different size, 
//the ones longer are also counted as different chars
//for example: "abc" vs. "acc" return 1;
//"abc" vs "acb" return 2
//"abc" vs "ab" reurn 1
//"abc" vs "cab" return 3
unsigned int CompareStrings(const string& str1, const string& str2);
double MatchBarcodes(const SequenceString& seq, const SequenceString& barcode, const MatchMatrix* mm);
SequenceString Reverse(SequenceString& seq);

#endif
