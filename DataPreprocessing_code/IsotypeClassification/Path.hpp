#ifndef PATH_HPP
#define PATH_HPP

#include "pairwiseAlignment.hpp"
#include <vector>

using namespace std;

//this entry class is used to define each local alignment entry/path
class Path
{
public:
  Path();  //~Path() in this case we don't need a destructor since we don't have dynamic allocation of memeory.
  Path( const unsigned int _OptimalIndex[2], const unsigned int _startIndex [2], const double& _optimalValue );
  void SetOptimalIndex(unsigned int _OptimalIndex[2]);
  void SetStartIndex(unsigned int _startIndex [2]);
  void SetOptimalValue(const double& _optimalValue);
  const unsigned int* GetOptimalIndex() const ;
  const unsigned int* GetStartIndex() const;
  double GetOptimalValue() const;
  bool isTraced() const;
  void setTraced() ;
private:
  unsigned int c_optimalIndex[2];
  unsigned int c_startIndex[2];
  double c_optimalValue;
  bool c_traced;
};

class PathElementEntry
{
	public:
	PathElementEntry();
	~PathElementEntry();
	
	//copy constructor
	PathElementEntry(const PathElementEntry& pe);
	
	unsigned GetNumberOfPathes() const;
	unsigned GetPathID(unsigned index) const;
	double GetPathOptimalScore(unsigned index) const;
	unsigned GetPathOptimalIndexPattern(unsigned index) const;
	unsigned GetPathOptimalIndexSubject(unsigned index) const;
	//in here, we have to set up thing at the same time together
	// input index: which path this is for the entry in the pathElementEntry vector, like first, second, third....
	//			pathID: the path ID in the path vector
	//		optimalScore: the m
	void AddPathInfoNoCheckingDuplication(const unsigned& index, const unsigned& pathID, 
			const double& optimalScore, 
			const unsigned& optimalIndexPattern, const unsigned& optimalIndexSubject);
	//combining the extra one with this current one.
	void AddPathInfo(const PathElementEntry& pe);
	
	//remove a path (indicated by "index") from the records.
	void RemovePath(const unsigned index);
	//based on the pathid, look up the index of the path in the records.
	unsigned LookUpPathIndex(const unsigned& pathID) const ;

	string toString() const;
	private:
	PathElementEntry& operator = (const PathElementEntry& p) { return *this;}; //disable it  
	//pathID and c_optimalScore a pair to "rember" the information about
	//pathes the current entry engaged. each entry can be joining multiple
	//path, thanking to the cases where mismatch score and gap penalty are 
	//identical. it can be more one, more than 2 (might be rare)
	//unsigned c_pathNumber; //this is used to indicate how many pathes the current entry engaged.
	vector<unsigned> c_pathID;
	vector<double> c_optimalScore;//used to indicate up to this current entry the optimal score for the specific path
	vector<unsigned> c_optimalIndexPattern;//used to indicate where the optimal score on pattern is.
	vector<unsigned> c_optimalIndexSubject;//same thing, but on subject.
};
 
 #endif
