#include <algorithm>
#undef NDEBUG // always want assertions to fire
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

const double kTargetCoverage = .6827;

const double A = 35;
const double B = 6;
const double C = 5;

const int Nmax = A+B+C + 5*sqrt(A+B+C); // don't consider fluctuations beyond 5sigma
//const int Nmin = A+B+C - 5*sqrt(A+B+C); // don't consider fluctuations beyond 5sigma

const int invstep = 400; // step by .025 event. Do it this way around so we can work with arrays with integer indices. Increase for smoother plots

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
  std::vector<double> chisq_best(Nmax*invstep);

  // What's the maximum and minimum number of events that can be predicted?
  double minExp = A-B-C;
  double maxExp = A+B+C;
  if(hie != kEither){
    minExp = (hie == kNH) ? A-B+C : A-B-C;
    maxExp = (hie == kNH) ? A+B+C : A+B-C;
  }
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

    // TODO this should be configurable (is currently 'either')
    const double chisq_true = std::min(chisq(Nobs, Nexp(delta_true, kNH)),
                                       chisq(Nobs, Nexp(delta_true, kIH)));

    const double dchisq = chisq_true - chisq_best[i];

    expts.emplace_back(Nobs, prob, dchisq);
  } // end for Nobs

  return expts;
}

double interpolate(double x, double x0, double x1, double y0, double y1)
{
  return y0 + (x-x0)/(x1-x0) * (y1-y0);
}

double find_quantile(const std::vector<Expt>& expts)
{
  double prev_dchisq = 0;
  double cov = 0;

  for(const Expt& e: expts){
    if(cov+e.prob > kTargetCoverage){
      return interpolate(kTargetCoverage,
                         cov, cov+e.prob,
                         prev_dchisq, e.dchisq);
    }

    cov += e.prob;
    prev_dchisq = e.dchisq;
  }

  abort();
}

double fc_critical_value_single(double delta_mock,
                                Hierarchy hie_mock,
                                const std::vector<double>& chisq_best)
{
  std::vector<Expt> expts = mock_expts(delta_mock, hie_mock, chisq_best);

  std::sort(expts.begin(), expts.end()); // from low to high dchisq

  return find_quantile(expts);
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


double hc_critical_value_single(double delta_mock,
                                const std::vector<double>& chisq_best)
{
  // 50-50 mixture of NH and IH experiments (based on our prior)
  std::vector<Expt> exptsNH = mock_expts(delta_mock, kNH, chisq_best);
  std::vector<Expt> exptsIH = mock_expts(delta_mock, kIH, chisq_best);
  std::vector<Expt> expts;
  expts.insert(expts.end(), exptsNH.begin(), exptsNH.end());
  expts.insert(expts.end(), exptsIH.begin(), exptsIH.end());
  for(Expt& e: expts) e.prob *= .5;

  std::sort(expts.begin(), expts.end()); // from low to high dchisq

  return find_quantile(expts);
}

std::vector<double> hc_critical_values(const std::vector<double>& chisq_best)
{
  std::vector<double> dchisq_crit(kDeltaScanValues.size());

  for(int j = 0; j < kDeltaScanValues.size(); ++j){
    const double delta_true = kDeltaScanValues[j];

    dchisq_crit[j] = hc_critical_value_single(delta_true, chisq_best);
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

double evaluate_coverage_prof_single(double delta_true,
                                     Hierarchy hie_true,
                                     double dchisq_crit_nh,
                                     double dchisq_crit_ih,
                                     const std::vector<double>& chisq_best)
{
  double coverage = 0;

  const std::vector<Expt> expts = mock_expts(delta_true, hie_true, chisq_best);

  for(const Expt& e: expts){
    // Assuming this value of delta, which hierarchy has the best likelihood?
    const double chisq_NH = chisq(e.Nobs, Nexp(delta_true, kNH));
    const double chisq_IH = chisq(e.Nobs, Nexp(delta_true, kIH));

    // Use the appropriate critical value for the most likely hierarchy
    const double dchisq_crit = (chisq_NH < chisq_IH) ? dchisq_crit_nh : dchisq_crit_ih;

    if(e.dchisq <= dchisq_crit){
      coverage += e.prob;
    }
  }

  return coverage;
}

std::vector<double> evaluate_coverage_prof(Hierarchy hie_true,
                                           const std::vector<double>& dchisq_crit_nh,
                                           const std::vector<double>& dchisq_crit_ih,
                                           const std::vector<double>& chisq_best)
{
  std::vector<double> cov(kDeltaScanValues.size());

  for(int j = 0; j < kDeltaScanValues.size(); ++j){
    const double delta_true = kDeltaScanValues[j];

    cov[j] = evaluate_coverage_prof_single(delta_true, hie_true, dchisq_crit_nh[j], dchisq_crit_ih[j], chisq_best);
  } // end for delta_true

  return cov;
}

enum EMethod{
  kInvalid,
  kWilks,
  kFC,
  kHC, // Highland-Cousins
  kProf // "profiled"
};

int main(int argc, char** argv)
{
  EMethod method = kInvalid;
  if(argc > 1){
    if(std::string_view(argv[1]) == "wilks") method = kWilks;
    if(std::string_view(argv[1]) == "fc") method = kFC;
    if(std::string_view(argv[1]) == "hc") method = kHC;
    if(std::string_view(argv[1]) == "prof") method = kProf;
  }

  Hierarchy trueHie = kEither;
  if(argc > 2){
    if(std::string_view(argv[2]) == "nh") trueHie = kNH;
    if(std::string_view(argv[2]) == "ih") trueHie = kIH;
  }

  Hierarchy fitHie = kEither;
  if(argc > 3){
    if(std::string_view(argv[3]) == "nh") fitHie = kNH;
    if(std::string_view(argv[3]) == "ih") fitHie = kIH;
  }

  Hierarchy mockHie = kEither;
  if(argc > 4){
    if(std::string_view(argv[4]) == "nh") mockHie = kNH;
    if(std::string_view(argv[4]) == "ih") mockHie = kIH;
  }

  if(method == kInvalid || trueHie == kEither ||
     (mockHie == kEither && method == kFC)){
    std::cerr << "Usage: fc METHOD TRUEHIE FITHIE MOCKHIE" << std::endl
              << "  METHOD:  'wilks', 'fc', 'hc' or 'prof'" << std::endl
              << "  TRUEHIE: 'nh' or 'ih'. True hierarchy (to evaluate coverage w.r.t)" << std::endl
              << "  FITHIE:  'nh', 'ih' or 'either'. Hierarchy assumed when fitting (unused for 'prof')" << std::endl
              << "  MOCKHIE: 'nh' or 'ih'. Hierarchy used for mock experiments ('fc' only)" << std::endl;
    return 1;
  }

  // Precompute what the best chisq would be for each number of observed events
  const std::vector<double> chisq_best = precompute_chisq_best(fitHie);

  std::cerr << "Computing critical values..." << std::endl;

  std::vector<double> dchisq_crit;
  if(method == kWilks) dchisq_crit = wilks_critical_values();
  if(method == kFC) dchisq_crit = fc_critical_values(mockHie, chisq_best);
  if(method == kHC) dchisq_crit = hc_critical_values(chisq_best);

  std::vector<double> dchisq_crit_nh, dchisq_crit_ih;
  if(method == kProf){
    dchisq_crit = std::vector<double>(DeltaScanValues().size(), -1);
    dchisq_crit_nh = fc_critical_values(kNH, chisq_best);
    dchisq_crit_ih = fc_critical_values(kIH, chisq_best);
  }

  std::cerr << "Evaluating coverage..." << std::endl;

  std::vector<double> coverage;
  if(method != kProf){
    coverage = evaluate_coverage(trueHie, dchisq_crit, chisq_best);
  }
  else{
    coverage = evaluate_coverage_prof(trueHie, dchisq_crit_nh, dchisq_crit_ih, chisq_best);
  }

  std::cout << "delta_true\tcoverage\tdchisq_crit" << std::endl;

  for(int j = 0; j < kDeltaScanValues.size(); ++j){
    const double delta_true = kDeltaScanValues[j];
    std::cout << delta_true << "\t" << 100*coverage[j] << "\t" << dchisq_crit[j] << std::endl;
  } // end for delta_true
}
