#include <algorithm>
#undef NDEBUG // always want assertions to fire
#include <cassert>
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

enum Hierarchy
{
  kNH,
  kIH,
  kEither
};

double Nexp(double delta, Hierarchy hie)
{
  assert(hie != kEither); // meaningless here
  if(hie == kNH)
    return A-B*sin(delta)+C;
  else
    return A-B*sin(delta)-C;
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
std::vector<double> precompute_chisq_best(Hierarchy hie)
{
  assert(hie != kEither);

  std::vector<double> chisq_best(Nmax*invstep);

  // What's the maximum and minimum number of events that can be predicted?
  const double minExp = (hie == kNH) ? A-B+C : A-B-C;
  const double maxExp = (hie == kNH) ? A+B+C : A+B-C;
  // NB there's no number of events "between" the two curves so we don't have
  // to allow for that case specially.

  for(int i = 0; i < Nmax*invstep; ++i){
    const double Nobs = i/double(invstep);

    // Normal hierarchy assumption is currently hardcoded here in + sign on C
    if(Nobs < minExp) chisq_best[i] = chisq(Nobs, minExp);
    else if(Nobs < maxExp) chisq_best[i] = 0; // find a perfect fit somewhere
    else chisq_best[i] = chisq(Nobs, maxExp);
  }

  return chisq_best;
}

std::vector<double> wilks_critical_values()
{
  return std::vector<double>(kDeltaScanValues.size(), 1);
}

struct Expt
{
  Expt(double N, double p, double d = -1) : Nobs(N), prob(p), dchisq(d) {}
  bool operator<(const Expt& e) const {return dchisq < e.dchisq;}
  double Nobs;
  double prob;
  double dchisq;
};

std::vector<Expt> mock_expts(double delta_true,
                             Hierarchy hie_true,
                             const std::vector<double>& chisq_best)
{
  std::vector<Expt> expts;

  for(int i = 0; i < Nmax*invstep; ++i){
    const double Nobs = i/double(invstep);

    // Probablility of observing Nobs given truth
    const double prob = gaus(Nobs, Nexp(delta_true, hie_true))/double(invstep);

    const double chisq_true = chisq(Nobs, Nexp(delta_true, hie_true));

    const double dchisq = chisq_true - chisq_best[i];

    expts.emplace_back(Nobs, prob, dchisq);
  } // end for Nobs

  return expts;
}

double fc_critical_value_single(double delta_true,
                                Hierarchy hie_assumed,
                                const std::vector<double>& chisq_best)
{
  std::vector<Expt> expts = mock_expts(delta_true, hie_assumed, chisq_best);

  std::sort(expts.begin(), expts.end()); // from low to high dchisq

  double cov = 0;
  double crit = 0;
  for(const Expt& e: expts){
    cov += e.prob;
    if(cov > .6827) return e.dchisq; // TODO interpolate to previous value?
  }

  abort();
}

std::vector<double> fc_critical_values(Hierarchy assumed_hie,
                                       const std::vector<double>& chisq_best)
{
  std::vector<double> dchisq_crit(kDeltaScanValues.size());

  for(int j = 0; j < kDeltaScanValues.size(); ++j){
    const double delta_true = kDeltaScanValues[j];

    dchisq_crit[j] = fc_critical_value_single(delta_true, assumed_hie, chisq_best);
  } // end for delta_true (j)

  return dchisq_crit;
}

double evaluate_coverage_single(double delta_true,
                                Hierarchy hie_true,
                                double dchisq_crit,
                                const std::vector<double>& chisq_best)
{
  double coverage = 0;

  const std::vector<Expt> expts = mock_expts(delta_true, hie_true, chisq_best);

  for(const Expt& e: expts){
    if(e.dchisq <= dchisq_crit){
      coverage += e.prob;
    }
  }

  return coverage;
}

std::vector<double> evaluate_coverage(Hierarchy hie_true,
                                      const std::vector<double>& dchisq_crit,
                                      const std::vector<double>& chisq_best)
{
  std::vector<double> cov(kDeltaScanValues.size());

  for(int j = 0; j < kDeltaScanValues.size(); ++j){
    const double delta_true = kDeltaScanValues[j];

    cov[j] = evaluate_coverage_single(delta_true, hie_true, dchisq_crit[j], chisq_best);
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

  Hierarchy trueHie = kEither;
  if(argc > 2){
    if(std::string_view(argv[2]) == "nh") trueHie = kNH;
    if(std::string_view(argv[2]) == "ih") trueHie = kIH;
  }

  Hierarchy assumedHie = kEither;
  if(argc > 3){
    if(std::string_view(argv[3]) == "nh") assumedHie = kNH;
    if(std::string_view(argv[3]) == "ih") assumedHie = kIH;
  }

  if(method == kInvalid || trueHie == kEither || assumedHie == kEither){
    std::cerr << "Usage: fc METHOD TRUEHIE ASSUMEDHIE" << std::endl
              << "  METHOD:     'wilks' or 'fc'" << std::endl
              << "  TRUEHIE:    'nh' or 'ih'. True hierarchy (to evaluate coverage w.r.t)" << std::endl
              << "  ASSUMEDHIE: 'nh' or 'ih'. Assumed hierarchy in throwing and fitting expts" << std::endl;
    return 1;
  }

  // Precompute what the best chisq would be for each number of observed events
  const std::vector<double> chisq_best = precompute_chisq_best(assumedHie);

  std::cerr << "Computing critical values..." << std::endl;

  std::vector<double> dchisq_crit;
  if(method == kWilks) dchisq_crit = wilks_critical_values();
  if(method == kFC) dchisq_crit = fc_critical_values(assumedHie, chisq_best);

  std::cerr << "Evaluating coverage..." << std::endl;

  const std::vector<double> coverage = evaluate_coverage(trueHie, dchisq_crit, chisq_best);

  std::cout << "delta_true\tcoverage\tdchisq_crit" << std::endl;

  for(int j = 0; j < kDeltaScanValues.size(); ++j){
    const double delta_true = kDeltaScanValues[j];
    std::cout << delta_true << "\t" << 100*coverage[j] << "\t" << dchisq_crit[j] << std::endl;
  } // end for delta_true
}
