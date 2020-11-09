#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

const double A = 35;
const double B = 6;
const double C = 5;

const int Nmax = A+B+C + 5*sqrt(A+B+C); // don't consider fluctuations beyond 5sigma
//const int Nmin = A+B+C - 5*sqrt(A+B+C); // don't consider fluctuations beyond 5sigma

const int invstep = 50; // step by .02 event. Do it this way around so we can work with arrays with integer indices

std::vector<double> DeltaScanValues()
{
  std::vector<double> ret(201);
  for(int i = 0; i <= 200; ++i) ret[i] = 2*M_PI*i/200.;
  return ret;
}

const std::vector<double> kDeltaScanValues = DeltaScanValues();

double Nexp(double delta)
{
  return A-B*sin(delta)+C; // hardcoded normal hierarchy here
}

double chisq(double obs, double exp)
{
  return (obs-exp)*(obs-exp)/exp;
}

double gaus(double x, double mu)
{
  const double sigmasq = mu;

  return exp(-(x-mu)*(x-mu)/(2*sigmasq))/sqrt(2*M_PI*sigmasq);
}

// Precompute what the best chisq would be for each number of observed events
std::vector<double> precompute_chisq_best()
{
  std::vector<double> chisq_best(Nmax*invstep);

  for(int i = 0; i < Nmax*invstep; ++i){
    const double Nobs = i/double(invstep);

    // Normal hierarchy assumption is currently hardcoded here in + sign on C
    if(Nobs < A-B+C) chisq_best[i] = chisq(Nobs, A-B-C);
    else if(Nobs < A+B+C) chisq_best[i] = 0; // find a perfect fit somewhere
    else chisq_best[i] = chisq(Nobs, A+B+C);
  }

  return chisq_best;
}

std::vector<double> wilks_critical_values()
{
  return std::vector<double>(kDeltaScanValues.size(), 1);
}

std::vector<double> fc_critical_values(const std::vector<double>& chisq_best)
{
  std::vector<double> dchisq_crit(kDeltaScanValues.size());

  for(int j = 0; j < kDeltaScanValues.size(); ++j){
    const double delta_true = kDeltaScanValues[j];

    struct Expt{
      Expt(double d, double p) : dchisq(d), prob(p) {}
      bool operator<(const Expt& e) const {return dchisq < e.dchisq;}
      double dchisq;
      double prob;
    };

    std::vector<Expt> expts;

    for(int i = 0; i < Nmax*invstep; ++i){
      const double Nobs = i/double(invstep);

      // Probablility of observing Nobs given delta_true
      const double prob = gaus(Nobs, Nexp(delta_true))/double(invstep);

      const double chisq_true = chisq(Nobs, Nexp(delta_true));

      const double dchisq = chisq_true - chisq_best[i];

      expts.emplace_back(dchisq, prob);
    } // end for Nobs (i)

    std::sort(expts.begin(), expts.end()); // from low to high dchisq

    double cov = 0;
    double crit = 0;
    for(const Expt& e: expts){
      cov += e.prob;
      if(cov > .6827){
        crit = e.dchisq; // TODO interpolate to previous value?
        break;
      }
    } // end for e
    dchisq_crit[j] = crit;
  } // end for delta_true (j)

  return dchisq_crit;
}

std::vector<double> evaluate_coverage(const std::vector<double>& dchisq_crit,
                                      const std::vector<double>& chisq_best)
{
  std::vector<double> cov(kDeltaScanValues.size());

  for(int j = 0; j < kDeltaScanValues.size(); ++j){
    const double delta_true = kDeltaScanValues[j];

    double coverage = 0;

    for(int i = 0; i < Nmax*invstep; ++i){
      const double Nobs = i/double(invstep);

      // Probablility of observing Nobs given delta_true
      const double prob = gaus(Nobs, Nexp(delta_true))/double(invstep);

      const double chisq_true = chisq(Nobs, Nexp(delta_true));

      const double dchisq = chisq_true - chisq_best[i];

      if(dchisq <= dchisq_crit[j]) coverage += prob;
    } // end for Nobs

    cov[j] = coverage;
  } // end for delta_true

  return cov;
}

enum EMethod{
  kInvalid,
  kWilks,
  kFC
};

int main(int argc, char** argv)
{
  EMethod method = kInvalid;
  if(argc > 1){
    if(std::string_view(argv[1]) == "wilks") method = kWilks;
    if(std::string_view(argv[1]) == "fc") method = kFC;
  }

  if(method == kInvalid){
    std::cerr << "Usage: fc METHOD" << std::endl
              << "  METHOD: 'wilks' or 'fc'" << std::endl;
    return 1;
  }

  // Precompute what the best chisq would be for each number of observed events
  const std::vector<double> chisq_best = precompute_chisq_best();

  std::cerr << "Computing critical values..." << std::endl;

  std::vector<double> dchisq_crit;
  if(method == kWilks) dchisq_crit = wilks_critical_values();
  if(method == kFC) dchisq_crit = fc_critical_values(chisq_best);

  std::cerr << "Evaluating coverage..." << std::endl;

  const std::vector<double> coverage = evaluate_coverage(dchisq_crit, chisq_best);

  std::cout << "delta_true\tcoverage\tdchisq_crit" << std::endl;

  for(int j = 0; j < kDeltaScanValues.size(); ++j){
    const double delta_true = kDeltaScanValues[j];
    std::cout << delta_true << "\t" << 100*coverage[j] << "\t" << dchisq_crit[j] << std::endl;
  } // end for delta_true
}
