#include <chrono>
#include <cmath>
#include <iostream>

#include "TH1F.h"
#include "TMinFit.h"

int main()
{
  TH1::AddDirectory(kFALSE);

  const int nbins = 1000;
  const double xmin = 0.0;
  const double xmax = 100.0;

  TH1F data("hdata", "data", nbins, xmin, xmax);
  TH1F signal("hsig", "signal", nbins, xmin, xmax);
  TH1F background("hbkg", "background", nbins, xmin, xmax);

  for (int bin = 1; bin <= nbins; ++bin) {
    const double x = data.GetBinCenter(bin);
    const double s = 1000.0 * std::exp(-0.5 * std::pow((x - 40.0) / 8.0, 2));
    const double b = 200.0 + 40.0 * std::sin(x * 0.1);
    data.SetBinContent(bin, s + b + 10.0 * std::cos(x * 0.2));
    signal.SetBinContent(bin, s);
    background.SetBinContent(bin, b);
  }

  TMinFit fitter(&data, &signal, &background);

  const double p1 = 0.1;
  const double p2 = 1.0;
  const double p3 = -0.05;
  const double p4 = 0.5;
  const double cutoff = 50.0;
  double nevt[2] = {0.0, 0.0};
  double chi2 = 0.0;

  // Warm up once to avoid startup overhead.
  chi2 = fitter.Evaluate(p1, p2, p3, p4);
  fitter.GetNevt(p1, p2, p3, p4, cutoff, nevt);

  const int nIter = 1000;

  auto t0 = std::chrono::high_resolution_clock::now();
  for (int i = 0; i < nIter; ++i) {
    chi2 += fitter.Evaluate(p1, p2, p3, p4);
  }
  auto t1 = std::chrono::high_resolution_clock::now();

  for (int i = 0; i < nIter; ++i) {
    fitter.GetNevt(p1, p2, p3, p4, cutoff, nevt);
  }
  auto t2 = std::chrono::high_resolution_clock::now();

  const double durEvaluate = std::chrono::duration<double, std::micro>(t1 - t0).count();
  const double durGetNevt = std::chrono::duration<double, std::micro>(t2 - t1).count();

  std::cout << "Warm-up chi2=" << chi2 << " nevt0=" << nevt[0] << " nevt1=" << nevt[1] << "\n";
  std::cout << "Evaluate: total=" << durEvaluate << " us, avg=" << durEvaluate / nIter << " us/call\n";
  std::cout << "GetNevt: total=" << durGetNevt << " us, avg=" << durGetNevt / nIter << " us/call\n";

  return 0;
}
