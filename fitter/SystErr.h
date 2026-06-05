#include "SimpleCardReader.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

#include "TH1.h"


#define nRCTypes 8
#define nPIDTypes 2
#define nMMETypes 2
#define nPi0Types 4
#define nWallTypes 2

enum rc_types   { sgel4=0, sgeg4, sgml4, sgmg4, sge, sgm, mge, mgm };
enum pid_types  { sg=0, mg };
enum mme_types  { mme=0, nue };
enum pi0_types  { rc1st=0, rc2nd, pid1st, pid2nd };
enum wall_types { fv200=0, fvexp };

const std::string rc_names[nRCTypes]      = { "sgel4", "sgeg4", "sgml4", "sgmg4", "sge", "sgm", "mge", "mgm" };
const std::string pid1r_names[nPIDTypes]  = { "sg", "mg" };
const std::string pidmr_names[nPIDTypes]  = { "sg", "mg" };
const std::string mme_names[nMMETypes]    = { "mme", "nue" };
const std::string pi0_names[nPi0Types]    = { "rc1st", "rc2nd", "pid1st", "pid2nd" };
const std::string wall_names[nWallTypes]  = { "fv200", "fvexp" };


// global switch of samples
bool rc_fit_global = false;
bool pid1r_fit_global = false;
bool pidmr_fit_global = false;
bool mme_fit_global = false;
bool pi0_fit_global = false;

// switch of sub-samples
bool rc_fit[nRCTypes][nWallTypes];
bool pid1r_fit[nPIDTypes][nWallTypes];
bool pidmr_fit[nPIDTypes][nWallTypes];
bool mme_fit[nMMETypes][nWallTypes];
bool pi0_fit[nPi0Types][nWallTypes];

// input filenames
std::string rc_filename;
std::string pid1r_filename;
std::string pidmr_filename;
std::string mme_filename;
std::string pi0_filename;

// sample prefix
std::string rc_prefix;
std::string pid1r_prefix;
std::string pidmr_prefix;
std::string mme_prefix;
std::string pi0_prefix;

// sample output filename
std::string rc_output;
std::string pid1r_output;
std::string pidmr_output;
std::string mme_output;
std::string pi0_output;

// parameter intervals for dChi2=30
double rc_lo[nRCTypes][nWallTypes][4];
double rc_hi[nRCTypes][nWallTypes][4];
double pid1r_lo[nPIDTypes][nWallTypes][4];
double pid1r_hi[nPIDTypes][nWallTypes][4];
double pidmr_lo[nPIDTypes][nWallTypes][4];
double pidmr_hi[nPIDTypes][nWallTypes][4];
double mme_lo[nMMETypes][nWallTypes][4];
double mme_hi[nMMETypes][nWallTypes][4];
double pi0_lo[nPi0Types][nWallTypes][4];
double pi0_hi[nPi0Types][nWallTypes][4];

// cutoff values
double rc_cutoff[nRCTypes][nWallTypes];
double pid1r_cutoff[nPIDTypes][nWallTypes];
double pidmr_cutoff[nPIDTypes][nWallTypes];
double mme_cutoff[nMMETypes][nWallTypes];
double pi0_cutoff[nPi0Types][nWallTypes];

double gaxis[37] = { -5.0, -4.33, -3.67, -3.0, -2.6, -2.2, -1.8, -1.4, -1.0,
                     -0.9, -0.8,  -0.7,  -0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0,
                      0.1,  0.2,   0.3,   0.4,  0.5,  0.6,  0.7,  0.8,  0.9, 
                      1.0,  1.4,   1.8,   2.2,  2.6,  3.0,  3.67, 4.33, 5.0  };

void LoadCardFile( const char * card );

template <int N>
void LoadCategory(SimpleCardReader* MasterCard,
                  const std::string sample_name,
                  std::string &sample_filename,
                  std::string &sample_prefix,
                  std::string &sample_output,
                  const std::string (&subsample_names)[N],
                  double (&lo)[N][nWallTypes][4],
                  double (&hi)[N][nWallTypes][4],
                  bool (&fit)[N][nWallTypes],
                  bool &gfit,
                  double (&cutoff)[N][nWallTypes] )
{
  std::string key; key.clear();
  key = sample_name + "_input";
  if ( !MasterCard->GetKey( key, sample_filename ) ) {
    std::cout << "\nCannot find input filename for " << sample_name
              << " sample. Skipping this sample." << std::endl;
    return;
  }
  std::cout << "\n" << key << "  " << sample_filename << std::endl;

  key.clear();
  key = sample_name + "_prefix";
  if ( !MasterCard->GetKey( key, sample_prefix ) ) {
    std::cout << "\nCannot find input prefix for " << sample_name
              << " sample. Skipping this sample." << std::endl;
    return;
  }
  std::cout << "\n" << key << "  " << sample_prefix << std::endl;

  std::string tmp_string;
  double * tmp_farray;

  key.clear();
  key = sample_name + "_cut";
  if ( MasterCard->GetKey( key, tmp_string ) ) {
    MasterCard->ParseArray(tmp_string, ",", &tmp_farray);
    for ( int i=0; i<N; ++i ) {
      for ( int j=0; j<nWallTypes; ++j ) {
        cutoff[i][j] = tmp_farray[i*nWallTypes+j+1];
        std::cout << key << "  " << subsample_names[i] << "_" << wall_names[j] << "  " << cutoff[i][j] << std::endl;
      }
    }
  } else {
    std::cout << "\nCannot find cutoff values for " << sample_name
              << " sample. Skipping this sample." << std::endl;
    return;
  }

  gfit = true;

  key.clear();
  key = sample_name + "_output";
  if ( !MasterCard->GetKey( key, sample_output ) ) {
    std::cout << "\nOutput filename not specified for " << sample_name
              << " sample. Using default " << sample_prefix << ".root" << std::endl;
    sample_output = sample_prefix + ".root";
  }
  std::cout << "\n" << key << "  " << sample_output << std::endl;

  static const char* suffixes[4] = { "_a1", "_b1", "_a2", "_b2" };

  for (int i=0; i<N; ++i)
  for (int j=0; j<nWallTypes; ++j)
  for ( int k=0; k<4; ++k ) {

    key.clear();
    key = subsample_names[i] + "_" + wall_names[j] + suffixes[k];
    if ( MasterCard->GetKey( key, tmp_string ) ) {
      MasterCard->ParseArray(tmp_string, ",", &tmp_farray);
      lo[i][j][k] = tmp_farray[1];
      hi[i][j][k] = tmp_farray[2];
      std::cout << key << "  " << tmp_string << std::endl;
    } else {
      fit[i][j] = false;
      std::cout << "Cannot find key " << key << ". Skipping this sample." << std::endl;
      break;
    }
  }

  return;
}

template <int N>
void GridScan( const std::string sample_filename,
               const std::string sample_prefix,
               const std::string sample_output,
               const std::string (&subsample_names)[N],
               double (&lo)[N][nWallTypes][4],
               double (&hi)[N][nWallTypes][4],
               bool (&fit)[N][nWallTypes],
               bool &gfit, double (&cutoff)[N][nWallTypes] );

void Scanner( TH1F * hdata, TH1F * hsig, TH1F * hbkg,
              double (&loc_axis)[4][37],
              double (&loc_chi2)[4][37], 
              double cutoff,
              TH1F * hposi, TH1F * hnega );
