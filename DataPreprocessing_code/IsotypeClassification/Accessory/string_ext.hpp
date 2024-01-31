#ifndef STRING_EXT_HPP
#define STRING_EXT_HPP
#include <string>
#include <vector>

#include <stdio.h>
#include <stdlib.h>

using std::vector;
using std::string;

void chomp_front_ext(string & s);
void chomp_end_ext(string & s);
void chomp_ext(string & s);

unsigned chomp_front_ext(char* s, const unsigned& len);
unsigned chomp_end_ext(char* s, const unsigned& len);
unsigned chomp_ext(char* s, const unsigned& len);


//the caller has to deallocate the memory!!!
//return the number of elements parsed
int split_ext(const string& s, string* buf, const char& delim);

//check whether the input str is a number format
//return 1 for a good integer number or a float format with trailing zeros follow '.'
//return 2 for a float integer number
//return 0 for a non-number format,
//no scientific format allowed
//+, - are consider as good format
int is_number(const char* str);

//check whether it is one of the following
//1,2,3,4,5,6,7,8,9,0
bool is_number_char( char c);

//change the string to upper/lower case, but only for alphabetical characters
string to_upper_str(const string& s);
string to_lower_str(const string& s);
//Definition:
//input --- string s, contains string of format #-#:# for parsing
//          char1 and char2, the char1 is for consecutive char such as for now using '-'
//                           char2 is for including char ':'
//output --- vector vec, the vector contains the numbers parsed from string
//           int, the whether there is error for the string
int parseNumberString(const string& s, vector<int>& vec, char stopping_char1='-', char stopping_char2=':');


//this is simply return the string class compare function, we need
//this because we want to use it in the sort function
bool stringCompare_ext(const string& s1, const string& s2);

//this is the function to return a string with a reverse order of the char in the string
string flipStr(const string&);

char DnaComplement(const char& c);

//return whether a char is a valid IUPAC nucleotide code. 
//Note: it doesn't include gap '-' or '.' as a valid code.
bool is_nucleotide(const char &c);

#endif
