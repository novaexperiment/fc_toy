#include "TGraph.h"
#include "TH2.h"
#include "TPad.h"

#include <fstream>

std::string junks;

void plot(std::string fname, std::string outname = "")
{
  std::ifstream fin(fname);

  fin >> junks >> junks; // header line

  TGraph* g = new TGraph;

  while(fin){
    double delta, coverage;
    fin >> delta >> coverage;

    g->SetPoint(g->GetN(), delta, coverage);
  }

  (new TH2F("", ";#delta_{true};Coverage (%)", 100, 0, 2*M_PI, 100, 50, 100))->Draw();

  g->SetLineWidth(2);
  g->Draw("l same");

  TLine* lin = new TLine(0, 68.27, 2*M_PI, 68.27);
  lin->SetLineWidth(2);
  lin->SetLineStyle(7);
  lin->Draw();

  if(!outname.empty()) gPad->Print(outname.c_str());
}
