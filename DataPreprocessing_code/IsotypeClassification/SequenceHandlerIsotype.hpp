#ifndef SEQUENCEHANDLER_ISOTYPE_HPP
#define SEQUENCEHANDLER_ISOTYPE_HPP

#include <vector>
#include <string>
#include "Accessory/SequenceString.hpp"
#include "score.hpp"

using namespace std;

enum mapType{ FivePrime, ThreePrime};
//IN this function, we take care to do the isotype identification
void MappingIsotypes(const vector<SequenceString>& _vecSeq, /*this is the sequence data that we want to find isotypes, is R1 if we do pair end.*/
		     vector<SequenceString>& _vecIsotype, /*this is the isotype sequences*/
		     const mapType& type, /*indicating whether it is 5'prime or 3' prime*/
		     ScoreMatrix* _sm, const double& _gapOpen, const double& _gapExtension,
		     const double& _mismatchRateThreshold, const unsigned _minimumOverlapLength,
		     const unsigned int& _offset, 
		     const string& _R1_fname,
		     const string& _R2_fname,
		     const bool& _demux=false,
			 const vector<string>& _vecSeq_Q=vector<string>(), /*read 1 sequence quality*/
			 const vector<SequenceString>& _vecSeq_R2=vector<SequenceString>(), /*R2 sequence data that we want to find isotypes*/
			 const vector<string>& _vecSeq_Q_R2=vector<string>()/*r2 quality.*/
		     //const string& _mapReverse_fname, const string& _mapNone_fname,
		     //const string& _mapCrossOver_fname, const string& _mapBreakOutsideCon_fname
		      ); 


/*void SetUpTrimFlag(const bool& _f);*/
/*void SetUpByIsotypeOutputFlag(const bool& _f);*/
#endif

