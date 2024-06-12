#include "TGraph.h"
#include "TH2.h"
#include "TPad.h"

#include <fstream>

std::string junks;

void plot(std::string fname, std::string covname = "", std::string critname = "")
{
  std::ifstream fin(fname);

  fin >> junks >> junks >> junks; // header line

  TGraph* gcov = new TGraph;
  TGraph* gcrit = new TGraph;

  while(fin){
    double delta, coverage, crit;
    fin >> delta >> coverage >> crit;

    gcov->SetPoint(gcov->GetN(), delta, coverage);
    gcrit->SetPoint(gcrit->GetN(), delta, crit);
  }

  (new TH2F("", ";#delta_{true};Coverage (%)", 100, 0, 2*M_PI, 100, 0, 100))->Draw();

  gcov->SetLineWidth(2);
  gcov->Draw("l same");

  TLine* lin = new TLine(0, 68.27, 2*M_PI, 68.27);
  lin->SetLineWidth(2);
  lin->SetLineStyle(7);
  lin->Draw();

  if(!covname.empty()) gPad->Print(covname.c_str());


  (new TH2F("", ";#delta_{true};#Delta#chi^{2}_{crit}", 100, 0, 2*M_PI, 100, 0, 2))->Draw();

  gcrit->SetLineWidth(2);
  gcrit->Draw("l same");

  lin = new TLine(0, 1, 2*M_PI, 1);
  lin->SetLineWidth(2);
  lin->SetLineStyle(7);
  lin->Draw();

 if(!critname.empty()) gPad->Print(critname.c_str());
}
