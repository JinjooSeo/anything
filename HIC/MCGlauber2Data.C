#include "DrawTool.C"

//--Parameter Setting------------------------------------------------------
double k = 1.5; // at 5.02TeV, controls the width, fpar[0]
double mu = 46.4; // at 5.02TeV, mean multiplicity per ancestor, fpar[1]
double f = 0.801; // at 5.02TeV
bool DrawDist = true; //Draw spectrums on canvas

Double_t fNBD(double *x, double *fpar);

void MCGlauber2Data(){
  TFile *inputMCGlauber = TFile::Open("MCGlauber-PbPb-5020GeV-b0-18fm.root","read");
  TTree *GlauberTree = (TTree*) inputMCGlauber->Get("lemon");

  ofstream outFile;
  outFile.open("Geometry_PbPb_5020GeV.csv");
  outFile << "b\tNpart\tNcoll\n";

  Int_t Npart, Ncoll, Nanc;
  Float_t b;
  GlauberTree->SetBranchAddress("npart",&Npart);
  GlauberTree->SetBranchAddress("ncoll",&Ncoll);
  GlauberTree->SetBranchAddress("b",&b);

  TH1D *hbDist = new TH1D("hbDist","",100,0,18);
  TH1D *hNpart = new TH1D("hNpart","",100,0,420);
  TH1D *hnHitsDist = new TH1D("hnHitsDist","",100,0,13000);
  TH2D *hNpartvsb = new TH2D("hbvsNpart","",40,0,400,18,0,18);
  TH2D *hNcollvsb = new TH2D("hbvsNcoll","",20,0,2000,18,0,18);

  TF1 *pdf_NBD = new TF1("fNBD",fNBD,0,50,2);
  pdf_NBD->SetParameters(k,mu);

  for(int i=0; i<GlauberTree->GetEntriesFast(); i++){
    GlauberTree->GetEntry(i);
    hbDist->Fill(b);
    hNpart->Fill(Npart);
    hNpartvsb->Fill(Npart,b);
    hNcollvsb->Fill(Ncoll,b);

    Double_t nHits = 0;
    Nanc = round(f*Npart + (1-f)*Ncoll);
    for (int j=0; j<Nanc; j++){
      nHits += pdf_NBD->GetRandom();
    }
    hnHitsDist->Fill(nHits);

    outFile << b << "\t" << Npart << "\t" << Ncoll << "\n";
  }

  if(DrawDist){  //----------------------------------------------Draw histogram on canvas
    SetStyle();
    int nCan = 5;
    TCanvas **can = new TCanvas*[nCan];
    for(int i=0; i<nCan; i++) can[i] = new TCanvas(Form("Can%d",i),"",650,500);

    can[0]->cd();
    hbDist->Draw(); SetAxis(hbDist,"b","entries");  //Draw b dirtribution
    can[1]->cd();
    gPad->SetLogy();
    hnHitsDist->Draw(); SetAxis(hnHitsDist,"nHits","entries");  //Draw b dirtribution
    can[2]->cd();
    hNpartvsb->Draw("COLZ"); SetAxis(hNpartvsb,"Npart","b");  //Draw Npart vs b
    can[3]->cd();
    hNcollvsb->Draw("COLZ"); SetAxis(hNcollvsb,"Ncoll","b");  //Draw Ncoll vs b
    can[4]->cd();
    gPad->SetLogy();
    hNpart->Draw(); SetAxis(hNpart,"N_{part}","entries");

    can[0]->SaveAs("ImpactParameter.pdf");
    can[1]->SaveAs("propV0Amp.pdf");
    can[2]->SaveAs("Npartvsb.pdf");
    can[3]->SaveAs("Ncollvsb.pdf");
    can[4]->SaveAs("Npart.pdf");

  }

  outFile.close();
}

Double_t fNBD(double *x, double *fpar)  // n hits per ancestor
{
  Float_t xx = x[0];
  Double_t gamterm = tgamma(xx+fpar[0])/(tgamma(xx+1)*tgamma(fpar[0]));
  Double_t powterm = pow(fpar[1]/fpar[0],xx)/pow((fpar[1]/fpar[0])+1,xx+fpar[0]);
  Double_t nbd = gamterm*powterm;

  return nbd;
}
