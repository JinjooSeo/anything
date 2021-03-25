#include <iostream>
using std::cout;
using std::endl;

#include "TH1D.h"
#include "TF1.h"

void DrawHW2(){

  gStyle->SetPadTickX(0);
  gStyle->SetPadTickY(0);
  gStyle->Reset("Plain");
  gStyle->SetOptTitle(0);
  gStyle->SetOptStat(0);
  //gStyle->SetTextFont(42);

  TCanvas* canvas = new TCanvas("canvas","Dispersion",700,500);
  TPad *pad = new TPad("pad","",0,0,1,1);
  pad->Draw();
  pad->cd();

  TH1F *hr = pad->DrawFrame(0.74,1.74,0.9,1.95);
  pad->GetFrame()->SetFillColor(0);

  Double_t xmin = 5;
  Double_t xmax = 10000;

  TF1* gr1 = new TF1("f1","x/2",xmin,xmax);
  gr1->Draw("");
  gr1->GetXaxis()->SetTitle("#sqrt{s_{NN}} (GeV)");
  gr1->GetYaxis()->SetTitle("Lorentz Factor");

  gPad->SetLogy();
  gPad->SetLogx();
  gPad->SetGridy();
  gPad->SetFillStyle(4000);
  gPad->SetFrameFillStyle(4000);
  gStyle->SetTitleSize(0.04,"xyz");
   gStyle->SetTextFont(42);
  canvas->cd();

  TPad *overlay = new TPad("overlay","",0,0,1,1);
  overlay->SetFillStyle(4000);
  overlay->SetFillColor(0);
  overlay->SetFrameFillStyle(0);
  overlay->Draw();
  overlay->cd();
  Double_t xmin2 = gr1->Eval(pow(xmin,2)/2);
  Double_t ymin2 = pad->GetUymin();
  Double_t xmax2 = gr1->Eval(pow(xmax,2)/2);
  Double_t ymax2 = pad->GetUymax();
  TH1F *hframe = overlay->DrawFrame(xmin2,ymin2,xmax2,ymax2);
  hframe->GetXaxis()->SetLabelOffset(99);
  hframe->GetYaxis()->SetLabelOffset(99);
  hframe->GetYaxis()->SetTickLength(0);
  hframe->GetXaxis()->SetTickLength(0);
  gPad->SetLogx(1);
  gPad->SetGridx();


  TGaxis *axis = new TGaxis(xmin2,ymax2,xmax2,ymax2,xmin2,xmax2,510,"G-");
  axis->SetTitleOffset(0);
  axis->SetTitle("E_{lab}/A (GeV)");
  gStyle->SetTextFont(42);
  axis->Draw();
}
