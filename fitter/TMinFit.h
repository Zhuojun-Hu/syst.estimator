#ifndef TMINFIT_H
#define TMINFIT_H

#include <vector>

class TH1;

class TMinFit
{
public:

  TMinFit(
    TH1* data,
    TH1* signal,
    TH1* background);

  ~TMinFit();

  // Prevent accidental copying
  TMinFit(const TMinFit&) = delete;
  TMinFit& operator=(const TMinFit&) = delete;

  struct Result {
    double chi2;
    std::vector<double> pars;
    std::vector<double> errs;
  };

  Result RunFit(int mode);

  double Evaluate(
    double p1,
    double p2,
    double p3,
    double p4);


  void GetNevt(
    bool isNorm,
    double p1,
    double p2,
    double p3,
    double p4,
    double cutoff,
    double (&nevt)[2]);

  double GetChi2() const { return best_chi2; }
  std::vector<double> GetPars() const { return best_pars; }
  std::vector<double> GetErrs() const { return best_errs; }

private:

  TH1* data_;
  TH1* signal_;
  TH1* background_;

  TH1* shiftsig_;
  TH1* shiftbkg_;

  std::vector<double> cache_;

  std::vector<double> sigaxis_;
  std::vector<double> bkgaxis_;

  int lowbin;
  int highbin;
  double lowedge;
  double highedge;
  double norm;
  double xmin;
  double xmax;
  double overlap;
  double binwidth;

  double best_chi2;
  std::vector<double> best_pars;
  std::vector<double> best_errs;

  void TransformHistogram(
    double p1,
    double p2,
    double p3,
    double p4);
};

#endif
