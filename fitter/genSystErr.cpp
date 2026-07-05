#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include "SimpleCardReader.h"
#include "SystErr.h"

#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TFile.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TTree.h"
#include "TParameter.h"
#include "TObject.h"
#include "TROOT.h"

#include "TMinFit.h"



int main( int argc, char * argv[] )
{
  if ( argc!=3 ) {
    std::cout << argv[0] << " [cardfile] [stage]\n";
    std::cout << "  [stage]   1: generate chi2 and # of event tree; 2: compute uncertainties" << std::endl;
    exit(-1);
  }

  std::string card = argv[1] ;

  LoadCardFile( card.c_str() );

  int stage = atoi(argv[2]);

  if (stage==1) {

    GridScan( rc_filename, 
              rc_prefix, 
              rc_output, 
              rc_names,  
              rc_lo, 
              rc_hi, 
              rc_fit, 
              rc_fit_global,
              rc_cutoff );

  } else if (stage==2) {

    GenSystErr( rc_filename,
                rc_prefix,
                rc_output,
                rc_names,
                rc_fit,
                rc_fit_global );
  }

  return 0;
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
  double loc_axis[4][Npoint], loc_chi2[4][Npoint];

  TTree *** tree;
  tree = new TTree ** [N];
  TGraph *** gchi2_p1, *** gchi2_p2, *** gchi2_p3, *** gchi2_p4;
  gchi2_p1 = new TGraph ** [N];
  gchi2_p2 = new TGraph ** [N];
  gchi2_p3 = new TGraph ** [N];
  gchi2_p4 = new TGraph ** [N];

  std::cout << " Starting GridScan for " << sample_prefix << " to generate # of events distribution" << std::endl;

  for ( int i=0; i<N; ++i ) {
    tree[i] = new TTree * [nWallTypes];
    gchi2_p1[i] = new TGraph * [nWallTypes];
    gchi2_p2[i] = new TGraph * [nWallTypes];
    gchi2_p3[i] = new TGraph * [nWallTypes];
    gchi2_p4[i] = new TGraph * [nWallTypes];
  for ( int j=0; j<nWallTypes; ++j ) {
    if ( !fit[i][j] ) continue;

    std::cout << " Processing " << sample_prefix << "_" << subsample_names[i] << "_" << wall_names[j] << std::endl;

    histname = sample_prefix + "_data_" + subsample_names[i] + "_" + wall_names[j];
    hdata = (TH1F*) f_in->Get( histname.c_str() );

    histname = sample_prefix + "_sig_" + subsample_names[i] + "_" + wall_names[j];
    hsig = (TH1F*) f_in->Get( histname.c_str() );
    hsig->Rebin(Nrebin);

    histname = sample_prefix + "_bkg_" + subsample_names[i] + "_" + wall_names[j];
    hbkg = (TH1F*) f_in->Get( histname.c_str() );
    hbkg->Rebin(Nrebin);
  
    histname = "t_" + sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j];
    tree[i][j] = new TTree( histname.c_str(), histname.c_str() );
    tree[i][j]->SetDirectory( f_out );

    for ( int k=0; k<Npoint; ++k ) {
      for ( int l=0; l<4; ++l ) {
        loc_axis[l][k] = (lo[i][j][l]+hi[i][j][l])/2.0 + gaxis[k]*(hi[i][j][l]-lo[i][j][l])/10.0;
        loc_chi2[l][k] = 1E10;
      }
    }
    Scanner( hdata, hsig, hbkg, loc_axis, loc_chi2, cutoff[i][j], tree[i][j] );

    gchi2_p1[i][j] = new TGraph( Npoint );
    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_chi2_a1";
    gchi2_p1[i][j]->SetName( histname.c_str() );
    gchi2_p1[i][j]->SetTitle( ";a1;#Delta#chi^{2}" );

    gchi2_p2[i][j] = new TGraph( Npoint );
    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_chi2_b1";
    gchi2_p2[i][j]->SetName( histname.c_str() );
    gchi2_p2[i][j]->SetTitle( ";b1;#Delta#chi^{2}" );

    gchi2_p3[i][j] = new TGraph( Npoint );
    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_chi2_a2";
    gchi2_p3[i][j]->SetName( histname.c_str() );
    gchi2_p3[i][j]->SetTitle( ";a2;#Delta#chi^{2}" );

    gchi2_p4[i][j] = new TGraph( Npoint );
    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_chi2_b2";
    gchi2_p4[i][j]->SetName( histname.c_str() );
    gchi2_p4[i][j]->SetTitle( ";b2;#Delta#chi^{2}" );

    for ( int k=0; k<Npoint; ++k ) {
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

  std::cout << " Exiting GridScan for " << sample_prefix << ". Trees written to " << f_out->GetName() << std::endl;

  return;
}

// scan grid defined by loc_axis, feed deltaChi2 to loc_chi2, fill # of event histograms (w/ & w/o weight)
void Scanner( TH1F * hdata, TH1F * hsig, TH1F * hbkg,
              double (&loc_axis)[4][Npoint],
              double (&loc_chi2)[4][Npoint], 
              double cutoff, TTree * tree
              )
{
  TMinFit * fitter = new TMinFit( hdata, hsig, hbkg );

  unsigned int iter = 0;
  double p1, p2, p3, p4, chi2_tmp, nevt[2], nevt_pos, nevt_neg, nevt_pos_nmn, nevt_neg_nmn;
  tree->Branch( "chi2",         &chi2_tmp,     "chi2/D" );
  tree->Branch( "nevt_pos",     &nevt_pos,     "nevt_pos/D" );
  tree->Branch( "nevt_neg",     &nevt_neg,     "nevt_neg/D" );

  fitter->GetNevt( 0, 0, 0, 0, 0, cutoff, nevt );
  nevt_pos_nmn = nevt[0];
  nevt_neg_nmn = nevt[1];
  
  TParameter<double>* nevt_pos_nmn_param = new TParameter<double>("NevtPosNominal", nevt_pos_nmn);
  TParameter<double>* nevt_neg_nmn_param = new TParameter<double>("NevtNegNominal", nevt_neg_nmn);
  tree->GetUserInfo()->Add(nevt_pos_nmn_param);
  tree->GetUserInfo()->Add(nevt_neg_nmn_param);

  for ( int ip1=0; ip1<Npoint; ++ip1 )
  for ( int ip2=0; ip2<Npoint; ++ip2 )
  for ( int ip3=0; ip3<Npoint; ++ip3 )
  for ( int ip4=0; ip4<Npoint; ++ip4 ) {
    p1 = loc_axis[0][ip1];
    p2 = loc_axis[1][ip2];
    p3 = loc_axis[2][ip3];
    p4 = loc_axis[3][ip4];

    chi2_tmp = fitter->Evaluate( p1, p2, p3, p4 );
    if ( chi2_tmp < loc_chi2[0][ip1] ) loc_chi2[0][ip1] = chi2_tmp;
    if ( chi2_tmp < loc_chi2[1][ip2] ) loc_chi2[1][ip2] = chi2_tmp;
    if ( chi2_tmp < loc_chi2[2][ip3] ) loc_chi2[2][ip3] = chi2_tmp;
    if ( chi2_tmp < loc_chi2[3][ip4] ) loc_chi2[3][ip4] = chi2_tmp;

    fitter->GetNevt( 0, p1, p2, p3, p4, cutoff, nevt );

    nevt_pos = nevt[0];
    nevt_neg = nevt[1];

    tree->Fill();

    ++iter;
    if ( iter%10000==0 ) {
      std::cout << " Iteration " << iter << " / " << Npoint*Npoint*Npoint*Npoint << " (" << 100.0*iter/(Npoint*Npoint*Npoint*Npoint) << "%)" << std::endl;
    }
  }

  double chimin = *std::min_element( &loc_chi2[0][0], &loc_chi2[0][0]+4*Npoint );
  for ( int i=0; i<4; ++i )
  for ( int j=0; j<Npoint; ++j ) {
    loc_chi2[i][j] -= chimin;
  }

  return;
}

// Read in the tree, generate histograms of # of events (w/ & w/o weight) and chi2 distribution,
// fit weighted histograms to get the systematic uncertainty, and write out the results to a new file.
template <int N>
void GenSystErr( const std::string sample_filename,
                 const std::string sample_prefix,
                 const std::string sample_output,
                 const std::string (&subsample_names)[N],
                 bool (&fit)[N][nWallTypes],
                 bool &gfit
               )
{
  if ( !gfit ) return;

  TFile * f_out = new TFile( sample_output.c_str(), "update" );
  if ( !f_out->IsOpen() ) {
    std::cout << " Error in GenSystErr: failed to open " << sample_output << std::endl;
    return;
  }

  std::string histname;
  TH1F *** hpos, *** hneg, *** hpos_wgt, *** hneg_wgt;
  hpos = new TH1F ** [N];
  hneg = new TH1F ** [N];
  hpos_wgt = new TH1F ** [N];
  hneg_wgt = new TH1F ** [N];
  TTree *** tree;
  tree = new TTree ** [N];
  TH1F *** histChi2;
  histChi2 = new TH1F ** [N];

  int foundBinId;
  double nevt_norm, prob_tmp;
  double chi2_tmp, chi2_min, nevt_pos, nevt_neg, nevt_pos_nmn, nevt_neg_nmn;

  std::cout << " Starting GenSystErr for " << sample_prefix << std::endl;
  std::ofstream outtxtfile;
  std::string outtxtname = sample_prefix + "_syst_err.txt";
  outtxtfile.open( outtxtname.c_str(), std::ios::out|std::ios::app );

  for ( int i=0; i<N; ++i ) {
    hpos[i] = new TH1F * [nWallTypes];
    hneg[i] = new TH1F * [nWallTypes];
    hpos_wgt[i] = new TH1F * [nWallTypes];
    hneg_wgt[i] = new TH1F * [nWallTypes];
    tree[i] = new TTree * [nWallTypes];
    histChi2[i] = new TH1F * [nWallTypes];
  for ( int j=0; j<nWallTypes; ++j ) {
    if ( !fit[i][j] ) continue;

    std::cout << " Processing " << sample_prefix << "_" << subsample_names[i] << "_" << wall_names[j] << std::endl;
  
    histname = "t_" + sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j];
    tree[i][j] = (TTree*)f_out->Get( histname.c_str() );

    if (tree[i][j]==nullptr) {
      std::cout << " Error in GenSystErr: failed to get tree " << histname << std::endl;
      continue;
    }
    double pos_max = tree[i][j]->GetMaximum("nevt_pos");
    double pos_min = tree[i][j]->GetMinimum("nevt_pos");
    double neg_max = tree[i][j]->GetMaximum("nevt_neg");
    double neg_min = tree[i][j]->GetMinimum("nevt_neg");
    double chi2_max = tree[i][j]->GetMaximum("chi2");
    double chi2_min = tree[i][j]->GetMinimum("chi2");

    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_positive";
    hpos[i][j] = new TH1F( histname.c_str(), histname.c_str(), 500, pos_min, pos_max );
    hpos[i][j]->SetDirectory( f_out );

    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_negative";
    hneg[i][j] = new TH1F( histname.c_str(), histname.c_str(), 500, neg_min, neg_max );
    hneg[i][j]->SetDirectory( f_out );

    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_positive_weighted";
    hpos_wgt[i][j] = new TH1F( histname.c_str(), histname.c_str(), 500, pos_min, pos_max );
    hpos_wgt[i][j]->SetDirectory( f_out );

    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_negative_weighted";
    hneg_wgt[i][j] = new TH1F( histname.c_str(), histname.c_str(), 500, neg_min, neg_max );
    hneg_wgt[i][j]->SetDirectory( f_out );

    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j] + "_chi2";
    histChi2[i][j] = new TH1F( histname.c_str(), histname.c_str(), NbinChi2, 0.0, chi2_max-chi2_min );
    histChi2[i][j]->SetDirectory( f_out );

    double ndf = 4;
    double probChi2[NbinChi2];
    for ( int k=0; k<NbinChi2; ++k ) {
      probChi2[k] = TMath::Prob(histChi2[i][j]->GetBinLowEdge(k+1), ndf) - TMath::Prob(histChi2[i][j]->GetBinLowEdge(k+2), ndf);
      if (k==NbinChi2-1) probChi2[k] = TMath::Prob(histChi2[i][j]->GetBinLowEdge(k+1), ndf);
    }

    double sumprob = 0.0;
    for ( int k=0; k<NbinChi2; ++k ) {
      sumprob += probChi2[k];
    }

    std::cout << "  Sum prob: " << sumprob << std::endl;

    tree[i][j]->SetBranchAddress( "chi2",         &chi2_tmp     );
    tree[i][j]->SetBranchAddress( "nevt_pos",     &nevt_pos     );
    tree[i][j]->SetBranchAddress( "nevt_neg",     &nevt_neg     );

    TObject* obj = tree[i][j]->GetUserInfo()->FindObject("NevtPosNominal");
    if (obj) {
      TParameter<double>* p_nevt_pos_nmn = (TParameter<double>*)obj;
      nevt_pos_nmn = p_nevt_pos_nmn->GetVal();
      std::cout << "NevtPosNominal: " << nevt_pos_nmn << std::endl;
    } else {
      std::cout << "NevtPosNominal object not found!" << std::endl;
    }
    obj = tree[i][j]->GetUserInfo()->FindObject("NevtNegNominal");
    if (obj) {
      TParameter<double>* p_nevt_neg_nmn = (TParameter<double>*)obj;
      nevt_neg_nmn = p_nevt_neg_nmn->GetVal();
      std::cout << "NevtNegNominal: " << nevt_neg_nmn << std::endl;
    } else {
      std::cout << "NevtNegNominal object not found!" << std::endl;
    }

    std::cout << "  Filling histograms..." << std::endl;
    for ( int entry=0; entry<tree[i][j]->GetEntries(); ++entry ) {
      tree[i][j]->GetEntry( entry );
      histChi2[i][j]->Fill( chi2_tmp-chi2_min );
      hpos[i][j]->Fill( nevt_pos );
      hneg[i][j]->Fill( nevt_neg );
    }

    // Must loop again to fill weighted histograms after histChi2 is filled
    for ( int entry=0; entry<tree[i][j]->GetEntries(); ++entry ) {
      tree[i][j]->GetEntry( entry );
      foundBinId = histChi2[i][j]->FindBin( chi2_tmp-chi2_min );
      nevt_norm = histChi2[i][j]->GetBinContent( foundBinId );
      if (nevt_norm>0) prob_tmp = probChi2[foundBinId-1] / nevt_norm;
      hpos_wgt[i][j]->Fill( nevt_pos, prob_tmp);
      hneg_wgt[i][j]->Fill( nevt_neg, prob_tmp);
    }

    // Fit the weighted histograms to get the systematic uncertainty
    histname = sample_prefix + "_" + subsample_names[i] + "_" + wall_names[j];
    DoGaussianFit( hpos_wgt[i][j], hneg_wgt[i][j], nevt_pos_nmn, nevt_neg_nmn, histname, outtxtfile );
  }
  }

  f_out->Write();
  f_out->Close();
  
  std::cout << " Exiting GenSystErr for " << sample_prefix << ". Histograms written to " << f_out->GetName() << std::endl;

  return;
}

void DoGaussianFit( TH1F * hpos_wgt, TH1F * hneg_wgt, double nevt_pos_nmn, double nevt_neg_nmn,
                    const std::string histname, std::ofstream & outtxtfile )
{
  double mean_range, half_range, peak, mean, sigma;
  double mean_pos, sigma_pos, mean_neg, sigma_neg;
  TCanvas * c1 = new TCanvas( "c1", "c1", 1200, 600 );
  c1->Divide(2,1);
  std::string padtitle;

  for ( int binId=1; binId<=500; ++binId ) {
    hpos_wgt->SetBinError( binId, 0.0 );
    hneg_wgt->SetBinError( binId, 0.0 );
  }

  // Gaussian fit for prob weighted positive distribution
  mean_range = hpos_wgt->GetMean();
  half_range = hpos_wgt->GetRMS() * 25.;
  hpos_wgt->GetXaxis()->SetRangeUser(mean_range-half_range,mean_range+half_range);
  padtitle = histname + " Positive Dist. (weighted)";
  hpos_wgt->SetTitle( padtitle.c_str() );
  c1->cd(1);
  hpos_wgt->GetYaxis()->SetTitleOffset(1.7);
  hpos_wgt->Draw();
  
  TF1 *fgaus_pos=new TF1("fgaus_pos","gaus",mean_range-half_range,mean_range+half_range);
  peak = hpos_wgt->GetMaximum()*sqrt(2.);
  mean = mean_range;
  sigma = hpos_wgt->GetRMS();
  fgaus_pos->SetParameters(peak, mean, sigma);
  hpos_wgt->Fit("fgaus_pos","RQW");
  peak = fgaus_pos->GetParameter(0);
  mean = fgaus_pos->GetParameter(1);
  sigma = fgaus_pos->GetParameter(2);
  fgaus_pos->SetRange(mean-3*sigma,mean+3*sigma);
  fgaus_pos->SetParameters(peak,mean,sigma);
  hpos_wgt->Fit("fgaus_pos","RQW");
  peak = fgaus_pos->GetParameter(0);
  mean = fgaus_pos->GetParameter(1);
  sigma = fgaus_pos->GetParameter(2);
  fgaus_pos->SetRange(mean-1.5*sigma,mean+1.5*sigma);
  fgaus_pos->SetParameters(peak,mean,sigma);
  hpos_wgt->Fit("fgaus_pos","RQW");
  mean_pos = fgaus_pos->GetParameter(1);
  sigma_pos = fgaus_pos->GetParameter(2);

  // Gaussian fit for prob weighted negative distribution
  mean_range = hneg_wgt->GetMean();
  half_range = hneg_wgt->GetRMS() * 25.;
  hneg_wgt->GetXaxis()->SetRangeUser(mean_range-half_range,mean_range+half_range);
  padtitle = histname + " Negative Dist. (weighted)";
  hneg_wgt->SetTitle( padtitle.c_str() );
  c1->cd(2);
  hneg_wgt->GetYaxis()->SetTitleOffset(1.7);
  hneg_wgt->Draw();
  
  TF1 *fgaus_neg=new TF1("fgaus_neg","gaus",mean_range-half_range,mean_range+half_range);
  peak = hneg_wgt->GetMaximum()*sqrt(2.);
  mean = mean_range;
  sigma = hneg_wgt->GetRMS();
  fgaus_neg->SetParameters(peak, mean, sigma);
  hneg_wgt->Fit("fgaus_neg","RQW");
  peak = fgaus_neg->GetParameter(0);
  mean = fgaus_neg->GetParameter(1);
  sigma = fgaus_neg->GetParameter(2);
  fgaus_neg->SetRange(mean-3*sigma,mean+3*sigma);
  fgaus_neg->SetParameters(peak,mean,sigma);
  hneg_wgt->Fit("fgaus_neg","RQW");
  peak = fgaus_neg->GetParameter(0);
  mean = fgaus_neg->GetParameter(1);
  sigma = fgaus_neg->GetParameter(2);
  fgaus_neg->SetRange(mean-1.5*sigma,mean+1.5*sigma);
  fgaus_neg->SetParameters(peak,mean,sigma);
  hneg_wgt->Fit("fgaus_neg","RQW");
  mean_neg = fgaus_neg->GetParameter(1);
  sigma_neg = fgaus_neg->GetParameter(2);

  c1->SaveAs( (histname+"_nhist_weighted.png").c_str() );

  outtxtfile << "=============================================================" << std::endl;
  outtxtfile << histname << std::endl;
  outtxtfile << "\tPositive: nominal = " << nevt_pos_nmn << ", mean = " << mean_pos << ", sigma = " << sigma_pos;
  outtxtfile << ", syst(abs sum) = " << (TMath::Abs(nevt_pos_nmn-mean_pos)+sigma_pos)/nevt_pos_nmn;
  outtxtfile << ", syst(qud sum) = " << TMath::Sqrt(TMath::Power(nevt_pos_nmn-mean_pos,2)+TMath::Power(sigma_pos,2))/nevt_pos_nmn << std::endl;
  outtxtfile << "\tNegative: nominal = " << nevt_neg_nmn << ", mean = " << mean_neg << ", sigma = " << sigma_neg;
  outtxtfile << ", syst(abs sum) = " << (TMath::Abs(nevt_neg_nmn-mean_neg)+sigma_neg)/nevt_neg_nmn;
  outtxtfile << ", syst(qud sum) = " << TMath::Sqrt(TMath::Power(nevt_neg_nmn-mean_neg,2)+TMath::Power(sigma_neg,2))/nevt_neg_nmn << std::endl;
  outtxtfile << "=============================================================\n" << std::endl;

  return;
}

//============================================================
// helper function to load values from card
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
