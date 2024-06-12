#include "TCanvas.h"
#include "TGraph.h"
#include "TH2.h"
#include "TPad.h"

#include <fstream>

std::string junks;

TGraph* load_coverage(std::string fname, int col, int style = kSolid)
{
  std::ifstream fin(fname);

  fin >> junks >> junks >> junks; // header line

  TGraph* gcov = new TGraph;

  while(fin){
    double delta, coverage, crit;
    fin >> delta >> coverage >> crit;

    gcov->SetPoint(gcov->GetN(), delta, coverage);
  }

  gcov->SetLineWidth(2);
  gcov->SetLineColor(col);
  gcov->SetLineStyle(style);

  return gcov;
}

TGraph* load_crit(std::string fname, int col, int style = kSolid)
{
  std::ifstream fin(fname);

  fin >> junks >> junks >> junks; // header line

  TGraph* gcrit = new TGraph;

  while(fin){
    double delta, coverage, crit;
    fin >> delta >> coverage >> crit;

    gcrit->SetPoint(gcrit->GetN(), delta, crit);
  }

  gcrit->SetLineWidth(2);
  gcrit->SetLineColor(col);
  gcrit->SetLineStyle(style);

  return gcrit;
}

void compare()
{
  TCanvas* c = new TCanvas;
  c->Divide(2, 2);

  for(int i = 0; i < 4; ++i){
    c->cd(i+1);

    std::string meth;
    std::string suffix;
    int col;
    if(i == 0){meth = "wilks"; col = kBlack;   suffix = "_either";}
    if(i == 1){meth = "fc";    col = kRed;     suffix = "_either_nh";}
    if(i == 2){meth = "hc";    col = kGreen+2; suffix = "_either";}
    if(i == 3){meth = "prof";  col = kBlue;}

    TGraph* gcov_nh = load_coverage(meth+"_nh"+suffix+".txt", col);
    TGraph* gcov_ih = load_coverage(meth+"_ih"+suffix+".txt", col, 7);

    (new TH2F("", (meth+";#delta_{true};Coverage (%)").c_str(), 100, 0, 2*M_PI, 100, 0, 100))->Draw();

    gcov_nh->Draw("l same");
    gcov_ih->Draw("l same");
  }

  c->cd(0);
  gPad->Print("comp_cov.pdf");

  /*
  TLine* lin = new TLine(0, 68.27, 2*M_PI, 68.27);
  lin->SetLineWidth(2);
  lin->SetLineStyle(7);
  lin->Draw();
  */

  new TCanvas;
  (new TH2F("", ";#delta_{true};Coverage (%)", 100, 0, 2*M_PI, 100, 30, 100))->Draw();
  TGraph* gcov_wilks_nh = load_coverage("wilks_nh_either.txt", kBlack);
  TGraph* gcov_wilks_ih = load_coverage("wilks_ih_either.txt", kBlack, 7);
  TGraph* gcov_fc_nh = load_coverage("fc_nh_either_nh.txt", kRed);
  TGraph* gcov_fc_ih = load_coverage("fc_nh_either_ih.txt", kRed, 7);

  TGraph* target = new TGraph;
  target->SetPoint(0, 0, 68.27);
  target->SetPoint(1, 2*M_PI, 68.27);
  target->SetLineStyle(2);
  target->Draw("l same");

  gcov_wilks_nh->Draw("l same");
  gcov_wilks_ih->Draw("l same");
  gcov_fc_nh->Draw("l same");
  gcov_fc_ih->Draw("l same");

  gPad->Print("comp_cov_overlay_basic.pdf");

  new TCanvas;
  (new TH2F("", ";#delta_{true};Coverage (%)", 100, 0, 2*M_PI, 100, 50, 90))->Draw();
  TGraph* gcov_hc_nh = load_coverage("hc_nh_either.txt", kGreen+2);
  TGraph* gcov_hc_ih = load_coverage("hc_ih_either.txt", kGreen+2, 7);
  TGraph* gcov_prof_nh = load_coverage("prof_nh.txt", kBlue);
  TGraph* gcov_prof_ih = load_coverage("prof_ih.txt", kBlue, 7);

  target->Draw("l same");

  gcov_hc_nh->Draw("l same");
  gcov_hc_ih->Draw("l same");
  gcov_prof_nh->Draw("l same");
  gcov_prof_ih->Draw("l same");

  gPad->Print("comp_cov_overlay.pdf");

  new TCanvas;
  TGraph* gwilks = load_crit("wilks_nh_either.txt", kBlack);
  TGraph* gfcnh = load_crit("fc_nh_either_nh.txt", kRed);
  TGraph* gfcih = load_crit("fc_nh_either_ih.txt", kRed, 7);
  TGraph* ghc = load_crit("hc_nh_either.txt", kGreen+2);

  (new TH2F("", ";#delta_{true};Critical value", 100, 0, 2*M_PI, 100, 0, 2))->Draw();

  gwilks->Draw("l same");
  gfcnh->Draw("l same");
  gfcih->Draw("l same");
  ghc->Draw("l same");

  gPad->Print("comp_crit.pdf");
}
