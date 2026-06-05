#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "SimpleCardReader.h"
#include "SystErr.h"

#include "TH1.h"
#include "TGraph.h"
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"

#include "TMinFit.h"



int main( int argc, char * argv[] )
{
  if ( argc!=2 ) {
    std::cout << argv[0] << " [cardfile]" << std::endl;
  }

  std::string card = argv[1] ;

  LoadCardFile( card.c_str() );

  GridScan( rc_filename, 
            rc_prefix, 
            rc_output, 
            rc_names,  
            rc_lo, 
            rc_hi, 
            rc_fit, 
            rc_fit_global,
            rc_cutoff );
}



//============================================================
void LoadCardFile( const char * card ) 
{
  // Initialize fitting switches to false.
  for ( int j=0; j<nWallTypes; ++j ) {
    for ( int i=0; i<nRCTypes; ++i ) {
      rc_fit[i][j] = true;
    }
    for ( int i=0; i<nPIDTypes; ++i ) {
      pid1r_fit[i][j] = true;
      pidmr_fit[i][j] = true;
    }
    for ( int i=0; i<nMMETypes; ++i ) {
      mme_fit[i][j] = true;
    }
    for ( int i=0; i<nPi0Types; ++i ) {
      pi0_fit[i][j] = true;
    }
  }

  SimpleCardReader * MasterCard = new SimpleCardReader( card );

  LoadCategory( MasterCard, 
                "rc",
                rc_filename, 
                rc_prefix,
                rc_output, 
                rc_names, 
                rc_lo, 
                rc_hi, 
                rc_fit, 
                rc_fit_global,
                rc_cutoff );

  LoadCategory(MasterCard, 
               "pid1r", 
               pid1r_filename, 
               pid1r_prefix,
               pid1r_output,
               pid1r_names, 
               pid1r_lo, 
               pid1r_hi, 
               pid1r_fit, 
               pid1r_fit_global,
               pid1r_cutoff );

  LoadCategory(MasterCard, 
               "mrpid", 
               pidmr_filename, 
               pidmr_prefix,
               pidmr_output,
               pidmr_names, 
               pidmr_lo, 
               pidmr_hi, 
               pidmr_fit, 
               pidmr_fit_global,
               pidmr_cutoff );

  LoadCategory(MasterCard, 
               "mme", 
               mme_filename, 
               mme_prefix,
               mme_output,
               mme_names, 
               mme_lo, 
               mme_hi, 
               mme_fit, 
               mme_fit_global,
               mme_cutoff );

  LoadCategory(MasterCard, 
               "pi0", 
               pi0_filename, 
               pi0_prefix,
               pi0_output,
               pi0_names, 
               pi0_lo, 
               pi0_hi, 
               pi0_fit, 
               pi0_fit_global,
               pi0_cutoff );

  delete MasterCard;
  return;
}

//============================================================
template <int N>
void GridScan( const std::string sample_filename,
               const std::string sample_prefix,
               const std::string sample_output,
               const std::string (&subsample_names)[N],
               double (&lo)[N][nWallTypes][4],
               double (&hi)[N][nWallTypes][4],
               bool (&fit)[N][nWallTypes],
               bool &gfit, double (&cutoff)[N][nWallTypes] )
{
  if ( !gfit ) return;

  TFile * f_in = new TFile( sample_filename.c_str(), "read" );
  if ( !f_in->IsOpen() ) {
    std::cout << " Error in GridScan: failed to open " << sample_filename << std::endl;
    return;
  }

  TFile * f_out = new TFile( sample_output.c_str(), "recreate" );
  if ( !f_out->IsOpen() ) {
    std::cout << " Error in GridScan: failed to open " << sample_output << std::endl;
    return;
  }

  std::string histname;
  TH1F * hdata, * hsig, * hbkg;
  double loc_axis[4][37], loc_chi2[4][37];
  TH1F *** hpost, *** hnega;
  hpost = new TH1F ** [N];
  hnega = new TH1F ** [N];
  TGraph *** gchi2_p1, *** gchi2_p2, *** gchi2_p3, *** gchi2_p4;
  gchi2_p1 = new TGraph ** [N];
  gchi2_p2 = new TGraph ** [N];
  gchi2_p3 = new TGraph ** [N];
  gchi2_p4 = new TGraph ** [N];


  for ( int i=0; i<N; ++i ) {
    hpost[i] = new TH1F * [nWallTypes];
    hnega[i] = new TH1F * [nWallTypes];
    gchi2_p1[i] = new TGraph * [nWallTypes];
    gchi2_p2[i] = new TGraph * [nWallTypes];
    gchi2_p3[i] = new TGraph * [nWallTypes];
    gchi2_p4[i] = new TGraph * [nWallTypes];
  for ( int j=0; j<nWallTypes; ++j ) {
    if ( !fit[i][j] ) continue;

    histname = sample_prefix + "_data_" + subsample_names[i] + "_" + wall_names[j];
    hdata = (TH1F*) f_in->Get( histname.c_str() );

    histname = sample_prefix + "_sig_" + subsample_names[i] + "_" + wall_names[j];
    hsig = (TH1F*) f_in->Get( histname.c_str() );

    histname = sample_prefix + "_bkg_" + subsample_names[i] + "_" + wall_names[j];
    hbkg = (TH1F*) f_in->Get( histname.c_str() );

    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_positive";
    hpost[i][j] = new TH1F( histname.c_str(), histname.c_str(), 500, 0.0, 1E5 );
    hpost[i][j]->SetDirectory( f_out );

    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_negative";
    hnega[i][j] = new TH1F( histname.c_str(), histname.c_str(), 500, 0.0, 1E5 );
    hnega[i][j]->SetDirectory( f_out );

    for ( int k=0; k<37; ++k ) {
      for ( int l=0; l<4; ++l ) {
        loc_axis[l][k] = lo[i][j][l] + gaxis[k]*(hi[i][j][l]-lo[i][j][l])/5.0;
        loc_chi2[l][k] = 1E10;
      }
    }
    Scanner( hdata, hsig, hbkg, loc_axis, loc_chi2, cutoff[i][j], hpost[i][j], hnega[i][j] );

    gchi2_p1[i][j] = new TGraph( 37 );
    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_chi2_a1";
    gchi2_p1[i][j]->SetName( histname.c_str() );
    gchi2_p1[i][j]->SetTitle( ";a1;#Delta#chi^{2}" );

    gchi2_p2[i][j] = new TGraph( 37 );
    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_chi2_b1";
    gchi2_p2[i][j]->SetName( histname.c_str() );
    gchi2_p2[i][j]->SetTitle( ";b1;#Delta#chi^{2}" );

    gchi2_p3[i][j] = new TGraph( 37 );
    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_chi2_a2";
    gchi2_p3[i][j]->SetName( histname.c_str() );
    gchi2_p3[i][j]->SetTitle( ";a2;#Delta#chi^{2}" );

    gchi2_p4[i][j] = new TGraph( 37 );
    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_chi2_b2";
    gchi2_p4[i][j]->SetName( histname.c_str() );
    gchi2_p4[i][j]->SetTitle( ";b2;#Delta#chi^{2}" );

    for ( int k=0; k<37; ++k ) {
      gchi2_p1[i][j]->SetPoint( k, loc_axis[0][k], loc_chi2[0][k] );
      gchi2_p2[i][j]->SetPoint( k, loc_axis[1][k], loc_chi2[1][k] );
      gchi2_p3[i][j]->SetPoint( k, loc_axis[2][k], loc_chi2[2][k] );
      gchi2_p4[i][j]->SetPoint( k, loc_axis[3][k], loc_chi2[3][k] );
    }

    f_out->WriteTObject( gchi2_p1[i][j] );
    f_out->WriteTObject( gchi2_p2[i][j] );
    f_out->WriteTObject( gchi2_p3[i][j] );
    f_out->WriteTObject( gchi2_p4[i][j] );
  }
  }

  f_out->Write();
  f_out->Close();
  f_in->Close();

  for ( int i=0; i<N; ++i ) 
  for ( int j=0; j<nWallTypes; ++j ) {
    if ( fit[i][j] ) {
      delete hpost[i][j];
      delete hnega[i][j];
      delete gchi2_p1[i][j];
      delete gchi2_p2[i][j];
      delete gchi2_p3[i][j];
      delete gchi2_p4[i][j];
    }
  }

  return;
}

void Scanner( TH1F * hdata, TH1F * hsig, TH1F * hbkg,
              double (&loc_axis)[4][37],
              double (&loc_chi2)[4][37], 
              double cutoff,
              TH1F * hposi, TH1F * hnega
              )
{
  TMinFit * fitter = new TMinFit( hdata, hsig, hbkg );

  double p1, p2, p3, p4, chi2_tmp, nevt[2];
  for ( int ip1=0; ip1<37; ++ip1 )
  for ( int ip2=0; ip2<37; ++ip2 )
  for ( int ip3=0; ip3<37; ++ip3 )
  for ( int ip4=0; ip4<37; ++ip4 ) {
    p1 = loc_axis[0][ip1];
    p2 = loc_axis[1][ip2];
    p3 = loc_axis[2][ip3];
    p4 = loc_axis[3][ip4];

    chi2_tmp = fitter->Evaluate( p1, p2, p3, p4 );
    if ( chi2_tmp < loc_chi2[0][ip1] ) loc_chi2[0][ip1] = chi2_tmp;
    if ( chi2_tmp < loc_chi2[1][ip2] ) loc_chi2[1][ip2] = chi2_tmp;
    if ( chi2_tmp < loc_chi2[2][ip3] ) loc_chi2[2][ip3] = chi2_tmp;
    if ( chi2_tmp < loc_chi2[3][ip4] ) loc_chi2[3][ip4] = chi2_tmp;

    fitter->GetNevt( p1, p2, p3, p4, cutoff, nevt );
    hposi->Fill( nevt[0] );
    hnega->Fill( nevt[1] );
  }

  double chimin = *std::min_element( &loc_chi2[0][0], &loc_chi2[0][0]+4*37 );
  for ( int i=0; i<4; ++i )
  for ( int j=0; j<37; ++j ) {
        loc_chi2[i][j] -= chimin;
  }

  delete fitter;
  return;
}
