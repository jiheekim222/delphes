//////////////////////////////////////////////////////////////
// 08/15/2023 Jihee Kim (jkim11@bnl.gov)
// Read HepMC file
//////////////////////////////////////////////////////////////
#include "HepMC3/GenEvent.h"
#include "HepMC3/ReaderAscii.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/Print.h"

#include "TH1F.h"
#include <iostream>
#include "TStyle.h"

using namespace HepMC3;

void read_hepmc(int pdgid = 11, const char* in_fname = "test.hepmc")
{
  // Setting for figures
  TStyle* kStyle = new TStyle("kStyle","Kim's Style");
  kStyle->SetOptStat("emr");
  kStyle->SetOptTitle(1);
  kStyle->SetOptFit(1);
  kStyle->SetStatColor(0);
  kStyle->SetStatW(0.15);
  kStyle->SetStatH(0.10);
  kStyle->SetStatX(0.85);
  kStyle->SetStatY(0.9);
  kStyle->SetStatBorderSize(1);
  kStyle->SetLabelSize(0.045,"xyz");
  kStyle->SetTitleSize(0.050,"xyz");
  kStyle->SetTitleOffset(1.2,"x");
  kStyle->SetTitleOffset(1.3,"y");
  kStyle->SetTitleOffset(1.2,"z");
  kStyle->SetLineWidth(2);
  kStyle->SetTitleFont(42,"xyz");
  kStyle->SetLabelFont(42,"xyz");
  kStyle->SetCanvasDefW(500);
  kStyle->SetCanvasDefH(500);
  kStyle->SetCanvasColor(0);
  kStyle->SetPadTickX(1);
  kStyle->SetPadTickY(1);
  kStyle->SetPadGridX(1);
  kStyle->SetPadGridY(1);
  kStyle->SetPadLeftMargin(0.15);
  kStyle->SetPadRightMargin(0.15);
  kStyle->SetPadTopMargin(0.1);
  kStyle->SetPadBottomMargin(0.15);
  TGaxis::SetMaxDigits(3);
  gStyle->SetPalette(1);
  gROOT->SetStyle("kStyle");

  ReaderAscii hepmc_input(in_fname);
  int         events_parsed = 0;
  GenEvent    evt(Units::GEV, Units::MM);

  // Book histograms
  TH1F* hEnergy = new TH1F("hEnergy",  "energy;E [GeV];Events",         100,-0.5,30.5);
  TH1F* hEta    = new TH1F("hEta",     "#eta;#eta;Events",              100,-10.0,10.0);
  TH1F* hTheta  = new TH1F("hTheta",   "#theta;#theta [degree];Events", 100,-0.5,180.5);
  TH1F* hPhi    = new TH1F("hPhi",     "#phi;#phi [degree];Events",     100,-180.5,180.5);
  TH2F* hPzPt   = new TH2F("hPzPt",    "pt vs pz;pt [GeV];pz [GeV]",    100,-0.5,30.5,100,-30.5,30.5);
  TH2F* hPxPy   = new TH2F("hPxPy",    "px vs py;px [GeV];py [GeV]",    100,-30.5,30.5,100,-30.5,30.5);
  TH3F* hP      = new TH3F("hP",       "p;px [GeV];py [GeV];pz [GeV]",  100,-30.5,30.5,100,-30.5,30.5,100,-30.5,30.5);
  TH1F* hMom    = new TH1F("hMom",     "momentum;p [GeV];Events",       100,-0.5,30.5);
	
  while(!hepmc_input.failed()) 
  {
    // Read event from input file
    hepmc_input.read_event(evt);
    // If reading failed - exit loop
    if( hepmc_input.failed() ) break;

    for(const auto& v : evt.vertices() ) 
    {
      // Loop over outgoing particles
      for(const auto& p : v->particles_out() ) 
      {
        if(p->pid() == pdgid) 
        {
          hEnergy->Fill(p->momentum().e());
	  hEta->Fill(p->momentum().eta());
	  hTheta->Fill(p->momentum().theta()*TMath::RadToDeg());
	  hPhi->Fill(p->momentum().phi()*TMath::RadToDeg());
	  hPzPt->Fill(TMath::Sqrt(p->momentum().px()*p->momentum().px() + p->momentum().py()*p->momentum().py()),p->momentum().pz());
	  hPxPy->Fill(p->momentum().px(),p->momentum().py());
	  hP->Fill(p->momentum().px(),p->momentum().py(),p->momentum().pz());
	  hMom->Fill(TMath::Sqrt(p->momentum().px()*p->momentum().px() + p->momentum().py()*p->momentum().py() + p->momentum().pz()*p->momentum().pz()));
	}
      }
    }
    evt.clear();
    events_parsed++;
  }
  std::cout << "Events parsed and written: " << events_parsed << std::endl;

  TCanvas *cnv1 = new TCanvas("cnv1", "cnv1", 3500, 1000);
  cnv1->Divide(4,2);
  cnv1->cd(1);
  hEnergy->SetLineWidth(2);
  hEnergy->Draw();
  cnv1->cd(2);
  hEta->SetLineWidth(2);
  hEta->Draw();
  cnv1->cd(3);
  hTheta->SetLineWidth(2);
  hTheta->Draw();
  cnv1->cd(4);
  hPhi->SetLineWidth(2);
  //hPhi->GetYaxis()->SetRangeUser(0.0,hPhi->GetMaximum()+300.0);
  hPhi->Draw();
  cnv1->cd(5);
  hPzPt->SetLineWidth(2);
  hPzPt->Draw("COLZ");
  cnv1->cd(6);
  hPxPy->SetLineWidth(2);
  hPxPy->Draw("COLZ");
  cnv1->cd(7);
  hP->SetLineWidth(2);
  hP->Draw();
  cnv1->cd(8);
  hMom->SetLineWidth(2);
  hMom->Draw();
}
