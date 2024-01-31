#ifndef SCORE_HPP
#define SCORE_HPP

//here in this file we define the score matrix for the most commonly 
//used ones



class ScoreMatrix
{
  //empty constructor
public :
  //ScoreMatrix(); this is not allowed in case someone call it and we end up with a uninitialized score matrix
  ScoreMatrix(const int* matrix[] , const double& scale, const char* alphabet, const int &alphabet_array_size);
  
  //we did not define the copy constructor, since the default one would be good.

  //destructor
  ~ScoreMatrix();

  //getter
  const int** GetScoreMatrix() const;  
  double GetScale() const ;
  const char* GetAlphabet() const;
  void SetScaledScoreFlag(const bool& Flag);
  
  double GetScore(const char& first, const char& second) const;
  const int GetAlphabetLength();
  //static members
  //this is the defualt static alphabet array for usering to use
  static const char AaAlphabet[]; //default amino acid alphabetic array
  static const char NucAlphabet[];//default Nucletide array.

  //member
private:		       
  const int** _matrix;
  double _scale;
  const char* _alphabet;
  bool _scaledScoreFlag;//this one is used to indicated wether we return
  //scaled score or not. default using unscaled score
  const int c_alphabet_array_size;
};

extern ScoreMatrix nuc44 ;
extern ScoreMatrix nuc44_v;//define as a variant of nuc44 for specific usage of D alignment with large cost for mismatch
extern ScoreMatrix blosum50;

extern ScoreMatrix tsm1;
extern ScoreMatrix tsm2;
extern ScoreMatrix nuc44HP;
extern ScoreMatrix nuc44DM1;//this is newly added to do isotype mapping with degenerate nts for isotype/adapter. following cutadapt format

//this following class is for sequence matching, NOT alignment.
//this is used for barcode demux
//we assume the barcode match occurs without gaps, at the fix position
//(begining and/or end of sequence). We simply match, instead of doing
//alignment. the only thing complicated is the degeneracy. It occurs 
//asymmetry. we always force some order. align sequence index by the 
//expected barcodes. allow degeneracy on barcodes,not much on sequence index
//   sequence index
//     ||||||
//     barcode
//a score in this matrix is between 0 and 1, with 0 being match, 1 mismatch. Could 
//   be a fraction meaning partial match due to degeneracy.
class MatchMatrix
{
  //empty constructor
public :
  //ScoreMatrix(); this is not allowed in case someone call it and we end up with a uninitialized score matrix
  MatchMatrix(const double* matrix[] , const char* alphabet, const int &alphabet_array_size);
  
  //we did not define the copy constructor, since the default one would be good.

  //destructor
  ~MatchMatrix();

  //getter
  const double** GetMatchMatrix() const;  
  //double GetScale() const ;
  const char* GetAlphabet() const;
  //void SetScaledScoreFlag(const bool& Flag);
  
  double GetScore(const char& first, const char& second) const;
  const int GetAlphabetLength();
  //static members
  //this is the defualt static alphabet array for usering to use
  //static const char AaAlphabet[]; //default amino acid alphabetic array <--not implemented yet
  static const char NucAlphabet[];//default Nucletide array.

  //member
private:		       
  const double** _matrix;
  //double _scale;
  const char* _alphabet;
  //bool _scaledScoreFlag;//this one is used to indicated wether we return
  //scaled score or not. default using unscaled score
  const int c_alphabet_array_size;
};

extern MatchMatrix mmf ; //matching matrix feng

#endif
