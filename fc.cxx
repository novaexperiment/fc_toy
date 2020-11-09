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

int main()
{
  // Precompute what the best chisq would be for each number of observed events
  std::vector<double> chisq_best(Nmax*invstep);

  for(int i = 0; i < Nmax*invstep; ++i){
    const double Nobs = i/double(invstep);

    // Normal hierarchy assumption is currently hardcoded here in + sign on C
    if(Nobs < A-B+C) chisq_best[i] = chisq(Nobs, A-B-C);
    else if(Nobs < A+B+C) chisq_best[i] = 0; // find a perfect fit somewhere
    else chisq_best[i] = chisq(Nobs, A+B+C);
  }

  const double dchisq_crit = 1;

  std::cout << "delta_true\tcoverage" << std::endl;

  for(double delta_true: DeltaScanValues()){
    double coverage = 0;

    for(int i = 0; i < Nmax*invstep; ++i){
      const double Nobs = i/double(invstep);

      // Probablility of observing Nobs given delta_true
      const double prob = gaus(Nobs, Nexp(delta_true))/double(invstep);

      const double chisq_true = chisq(Nobs, Nexp(delta_true));

      const double dchisq = chisq_true - chisq_best[i];

      if(dchisq <= dchisq_crit) coverage += prob;
    } // end for Nobs

    std::cout << delta_true << "\t" << 100*coverage << std::endl;
  } // end for delta_true
}
