// TMinFit.cpp

#include "TMinFit.h"

#include <iostream>
#include "TMinuit.h"
#include "TMath.h"
#include "TH1.h"

#include <stdexcept>

#define DEBUG 0

namespace
{
  TMinFit* gCurrentFitter = nullptr;

  void MinuitFCN(
    Int_t& npar,
    Double_t* grad,
    Double_t& fval,
    Double_t* par,
    Int_t flag)
  {
    fval = gCurrentFitter->Evaluate(
      par[0],
      par[1],
      par[2],
      par[3]);
  }
}


//======================================================================
// Constructor
//======================================================================

TMinFit::TMinFit(
  TH1* data,
  TH1* signal,
  TH1* background)
  :
  data_(data),
  signal_(signal),
  background_(background),
  shiftsig_(nullptr),
  shiftbkg_(nullptr)
{
  std::cout << "Constructing TMinFit ..." << std::endl;
  if (!data_)
    throw std::runtime_error("data histogram is null");

  if (!signal_)
    throw std::runtime_error("signal histogram is null");

  if (!background_)
    throw std::runtime_error("background histogram is null");

  shiftsig_ = (TH1*)data_->Clone("shiftsig");
  shiftsig_->SetDirectory(nullptr);
  shiftsig_->Reset();
  shiftbkg_ = (TH1*)data_->Clone("shiftbkg");
  shiftbkg_->SetDirectory(nullptr);
  shiftbkg_->Reset();

  const int nbins = signal_->GetNbinsX();

  if ( background_->GetNbinsX()!=nbins )
    throw std::runtime_error("signal/background bin mismatch");

  cache_.resize(nbins+1);
  sigaxis_.resize(nbins+1);
  bkgaxis_.resize(nbins+1);

  for (int bin = 0; bin<=nbins+1; ++bin) {
    cache_[bin] = signal_->GetBinLowEdge(bin+1);
  }

  norm = data_->Integral()/(signal_->Integral()+background_->Integral());
  //norm = 1.0;
  xmin = data_->GetBinLowEdge(1);
  xmax = data_->GetBinLowEdge(data_->GetNbinsX()+1);
  binwidth = signal_->GetBinWidth(1);
}

//======================================================================
// Destructor
//======================================================================

TMinFit::~TMinFit()
{

  std::cout << "Destroying TMinFit ..." << std::endl;
  if (shiftsig_) {
    delete shiftsig_;
    shiftsig_ = nullptr;
  }
  if (shiftbkg_) {
    delete shiftbkg_;
    shiftbkg_ = nullptr;
  }
}

//======================================================================
// Fit
//======================================================================

TMinFit::Result TMinFit::RunFit(int mode)
{
  Result result;

  // Clear buffers
  best_chi2 = 0.0;
  best_pars.clear();
  best_errs.clear();
  best_pars.reserve(4);
  best_errs.reserve(4);

  gCurrentFitter = this;

  TMinuit minuit(4);

  minuit.SetFCN(MinuitFCN);

  Double_t arglist[10];
  Int_t ierflg = 0;
  arglist[0] = 1;

  minuit.mnexcm("SET ERR", arglist, 1, ierflg );

  minuit.mnparm( 0, "a1", 0.0, 1e-4, -1.0, 1.0, ierflg );
  minuit.mnparm( 1, "b1", 0.0, 1e-4, -10., 10., ierflg );
  minuit.mnparm( 2, "a2", 0.0, 1e-4, -1.0, 1.0, ierflg );
  minuit.mnparm( 3, "b2", 0.0, 1e-4, -10., 10., ierflg );

  arglist[0] = 2e5;
  arglist[1] = 1e-3;

  minuit.mnexcm("MIGRAD", arglist, 2, ierflg );
  if (mode>1)
  minuit.mnexcm("HESse",  arglist, 2, ierflg );
  if (mode>2)
  minuit.mnexcm("MINOS",  arglist, 2, ierflg );

  Double_t amin,edm,errdef;
  Int_t nvpar,nparx,icstat;
  minuit.mnstat(amin,edm,errdef,nvpar,nparx,icstat);

  best_chi2 = amin;
  result.chi2 = amin;

  for(int ipar=0; ipar<4; ++ipar) {
    TString chnam;
    Double_t val;
    Double_t err;
    Double_t xlolim;
    Double_t xuplim;
    Int_t iuint;
    // minos errors
    Double_t eplus;
    Double_t eminus;
    Double_t eparab;
    Double_t gcc;
    minuit.mnpout(ipar, chnam, val, err, xlolim, xuplim, iuint);
    minuit.mnerrs(ipar, eplus, eminus, eparab, gcc);

    best_pars.push_back(val);
    best_errs.push_back(err);
    result.pars.push_back(val);
    result.errs.push_back(err);
  }

  return result;
}

//======================================================================
// Evaluate likelihood
//======================================================================

double TMinFit::Evaluate(
    double p1,
    double p2,
    double p3,
    double p4)
{
    //------------------------------------------------------------------
    // Update transformed histogram
    //------------------------------------------------------------------

  TransformHistogram( p1, p2, p3, p4);

  double chi2 = 0.0;
  double obs, exp;

  for (int bin = 1; bin<=data_->GetNbinsX(); ++bin) {
    obs = data_->GetBinContent( bin );
    exp = shiftsig_->GetBinContent( bin ) + shiftbkg_->GetBinContent( bin );

    //std::cout << "obs = " << obs << "  exp = " << exp << std::endl;

    if (exp<0) {
#if DEBUG>0
      std::cerr << "Negative expectation in bin " << bin << ": " << exp << "\n";
      std::cerr << "Setting bin " << bin << " to zero." << std::endl;
#endif
      exp = 0;
    }

    if (obs<=0 || exp==0) {
      chi2 += 2.0 * TMath::Abs(exp - obs);
    } else {
      chi2 += 2.0 * (exp - obs + obs*log(obs/exp));
    }
  }

  return chi2;
}

//======================================================================
// Histogram transformation
//======================================================================

void TMinFit::TransformHistogram(
    double p1,
    double p2,
    double p3,
    double p4)
{
  //------------------------------------------------------------------
  // Apply parameter-dependent transformation
  //------------------------------------------------------------------

  p1 = p1 + 1.0;
  p3 = p3 + 1.0;

  shiftsig_->Reset();
  shiftbkg_->Reset();  

  for (int bin = 0; bin<=signal_->GetNbinsX(); ++bin) {
    sigaxis_[bin] = cache_[bin] * p1 + p2;
    bkgaxis_[bin] = cache_[bin] * p3 + p4;
  }

  double bc = 0.0;

  for (int bin=1; bin<=signal_->GetNbinsX(); ++bin) {
    bc = signal_->GetBinContent(bin);
    if (bc<=0.) continue;
    bc *= norm;

    lowedge = sigaxis_[bin-1];
    highedge = sigaxis_[bin];

#if DEBUG > 1
    std::cout << "sigaxis_[bin] = " << lowedge << "\n";
#endif

    // handle left overflow bin
    if (highedge<=xmin) {
      shiftsig_->AddBinContent( 1, bc );
      continue;
    }

    // handle right overflow bin
    if (lowedge>=xmax)  {
      shiftsig_->AddBinContent( shiftsig_->GetNbinsX(), bc );
      continue;
    }

    lowbin = shiftsig_->FindBin(lowedge);
    highbin = shiftsig_->FindBin(highedge);

#if DEUBG > 1
    if (highbin - lowbin>1) std::cout << "highbin - lowbin = " << highbin - lowbin << "\n";
#endif    

    if (lowbin==highbin) {
      shiftsig_->AddBinContent( lowbin, bc );

    } else { // This is based on an assumption that signal_ has finner binnings than shiftsig_.

      overlap = (shiftsig_->GetBinLowEdge( highbin ) - lowedge) / binwidth;
      if (lowedge>xmin)  shiftsig_->AddBinContent( lowbin, bc*overlap );
      //if (overlap>1.) continue;
      if (highedge<xmax) shiftsig_->AddBinContent( highbin, bc*(1-overlap) );
    }
  }

  bc = 0.0;

  for (int bin=1; bin<=background_->GetNbinsX(); ++bin) {
    bc = background_->GetBinContent(bin);
    if (bc<=0.) continue;
    bc *= norm;

    lowedge = bkgaxis_[bin-1];
    highedge = bkgaxis_[bin];

#if DEBUG > 1
    std::cout << "bkgaxis_[bin] = " << lowedge << "\n";
#endif

    // handle left overflow bin
    if (highedge<=xmin) {
      shiftbkg_->AddBinContent( 1, bc );
      continue;
    }

    // handle right overflow bin
    if (lowedge>=xmax)  {
      shiftbkg_->AddBinContent( shiftsig_->GetNbinsX(), bc );
      continue;
    }

    lowbin = shiftbkg_->FindBin(lowedge);
    highbin = shiftbkg_->FindBin(highedge);

#if DEBUG > 1
    if (highbin - lowbin>1) std::cout << "highbin - lowbin = " << highbin - lowbin << "\n";
#endif

    if (lowbin==highbin) {
      shiftbkg_->AddBinContent( lowbin, bc );

    } else { // This is based on an assumption that background_ has finner binnings than shiftbkg_.

      overlap = (shiftbkg_->GetBinLowEdge( highbin ) - lowedge) / binwidth;
      if (overlap>1.) overlap = 1.0;
      if (lowedge>xmin)  shiftbkg_->AddBinContent( lowbin, bc*overlap );
      if (highedge<xmax) shiftbkg_->AddBinContent( highbin, bc*(1-overlap) );
    }
  }

#if DEBUG > 1
  for (int bin=1; bin<=shiftsig_->GetNbinsX(); ++bin) {
    std::cout << "sig = " << shiftsig_->GetBinContent(bin) << " bkg = " << shiftbkg_->GetBinContent(bin) << "\n";
  }
#endif

#if DEBUG >0
  std::cout << "signal_ sum = " << signal_->Integral() << "\n";
  std::cout << "background_ sum = " << background_->Integral() << "\n";
  std::cout << "shiftsig_ sum = " << shiftsig_->Integral() << "\n";
  std::cout << "shiftbkg_ sum = " << shiftbkg_->Integral() << "\n";
#endif

  return;
}

//======================================================================
// Retrieve positive/negative event counts for given parameters
//======================================================================

void TMinFit::GetNevt(
    bool isNorm,
    double p1,
    double p2,
    double p3,
    double p4,
    double cutoff,
    double (&nevt)[2])
{
  double loc_norm = 1.;
  if (isNorm) loc_norm = norm;

  // inverted transformation
  p1 = p1 + 1.0;
  p3 = p3 + 1.0;
  double cutoff_sig = (cutoff - p2) / p1;
  double cutoff_bkg = (cutoff - p4) / p3;

  int nbins[2] = { signal_->GetNbinsX(), background_->GetNbinsX() };
  int cutbin[2] = { signal_->FindBin(cutoff_sig), background_->FindBin(cutoff_bkg) };
  if ( cutoff_sig < signal_->GetBinLowEdge(1) || cutoff_sig > signal_->GetBinLowEdge(nbins[0]+1) ) {
    std::cerr << "Cutoff " << cutoff_sig << " is out of signal histogram range.\n";
    nevt[0] = 0.0;
  } else
  if ( cutoff_bkg < background_->GetBinLowEdge(1) || cutoff_bkg > background_->GetBinLowEdge(nbins[1]+1) ) {
    std::cerr << "Cutoff " << cutoff_bkg << " is out of background histogram range.\n";
    nevt[0] = 0.0;
  } else {
    nevt[0] = signal_->Integral( 1, cutbin[0] );
    nevt[0] -= signal_->GetBinContent( cutbin[0] )
            * (signal_->GetBinLowEdge( cutbin[0]+1 ) - cutoff_sig) / binwidth;
    nevt[0] += background_->Integral( 1, cutbin[1] );
    nevt[0] -= background_->GetBinContent( cutbin[1] )
            * (background_->GetBinLowEdge( cutbin[1]+1 ) - cutoff_bkg) / binwidth;
    nevt[0] *= loc_norm;
  }

  nevt[1] = signal_->Integral() + background_->Integral();
  nevt[1] *= loc_norm;
  nevt[1] -= nevt[0];

  return;
}

