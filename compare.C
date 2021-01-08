#include "TGraph.h"
#include "TH2.h"
#include "TPad.h"

#include <fstream>

std::string junks;

TGraph* load_coverage(std::string fname)
{
  std::ifstream fin(fname);

  fin >> junks >> junks >> junks; // header line

  TGraph* gcov = new TGraph;

  while(fin){
    double delta, coverage, crit;
    fin >> delta >> coverage >> crit;

    gcov->SetPoint(gcov->GetN(), delta, coverage);
  }

  return gcov;
}

TGraph* load_crit(std::string fname)
{
  std::ifstream fin(fname);

  fin >> junks >> junks >> junks; // header line

  TGraph* gcrit = new TGraph;

  while(fin){
    double delta, coverage, crit;
    fin >> delta >> coverage >> crit;

    gcrit->SetPoint(gcrit->GetN(), delta, crit);
  }

  return gcrit;
}

void compare()
{
  TCanvas* c = new TCanvas;
  c->Divide(2, 2);

  for(int i = 0; i < 4; ++i){
    c->cd(i+1);

    std::string meth;
    int col;
    if(i == 0){meth = "wilks"; col = kBlack;}
    if(i == 1){meth = "fc"; col = kRed;}
    if(i == 2){meth = "hc"; col = kGreen+2;}
    if(i == 3){meth = "prof"; col = kBlue;}

    TGraph* gcov_nh = load_coverage(meth+"_nh_either_nh.txt");
    TGraph* gcov_ih = load_coverage(meth+"_ih_either_nh.txt");

    (new TH2F("", (meth+";#delta_{true};Coverage (%)").c_str(), 100, 0, 2*M_PI, 100, 0, 100))->Draw();

    gcov_nh->SetLineColor(col);
    gcov_nh->SetLineWidth(2);
    gcov_nh->Draw("l same");
    gcov_ih->SetLineColor(col);
    gcov_ih->SetLineStyle(7);
    gcov_ih->SetLineWidth(2);
    gcov_ih->Draw("l same");
  }

  c->cd(0);
  gPad->Print("comp_cov.png");

  /*
  TLine* lin = new TLine(0, 68.27, 2*M_PI, 68.27);
  lin->SetLineWidth(2);
  lin->SetLineStyle(7);
  lin->Draw();
  */

  new TCanvas;
  TGraph* gwilks = load_crit("wilks_nh_either_nh.txt");
  TGraph* gfcnh = load_crit("fc_nh_either_nh.txt");
  TGraph* gfcih = load_crit("fc_nh_either_ih.txt");
  TGraph* ghc = load_crit("hc_nh_either_nh.txt");
  //  TGraph* gprof = load_crit("prof_nh_either_nh.txt");

  gwilks->SetLineWidth(2);
  gfcnh->SetLineWidth(2);
  gfcih->SetLineWidth(2);
  ghc->SetLineWidth(2);

  gfcnh->SetLineColor(kRed);
  gfcih->SetLineColor(kRed);
  gfcih->SetLineStyle(7);
  ghc->SetLineColor(kGreen+2);

  (new TH2F("", ";#delta_{true};Critical value", 100, 0, 2*M_PI, 100, 0, 2))->Draw();
  gwilks->Draw("l same");
  gfcnh->Draw("l same");
  gfcih->Draw("l same");
  ghc->Draw("l same");
  //  gprof->Draw("l same");

  gPad->Print("comp_crit.png");
}
