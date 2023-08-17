//////////////////////////////////////////
// 08/16/2023 Jihee Kim (jkim11@bnl.gov)
// Single 5GeV Pion+
//////////////////////////////////////////

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes.so)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

void singlePionPlus5GeV(const char *inputFile)
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

  gSystem->Load("libDelphes");
  // Create chain of root trees
  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);
  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  // Get pointers to branches used in this analysis
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray* branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");

  // Get total number of events
  Long64_t allEntries = treeReader->GetEntries();
  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Track *track, *eflowtrack;
  Tower *tower;
 
  // Book histograms
  TH1* hTrackDeltaEta = new TH1D("hTrackDeltaEta","-1.0<#eta<1.0;#eta^{rec}-#eta^{gen};Number of tracks for #pi+",100, -1.0, 1.0);
  TH1* hTrackDeltaPhi = new TH1D("hTrackDeltaPhi","-1.0<#eta<1.0;#phi^{rec}-#phi^{gen};Number of tracks for #pi+",100, -1.0, 1.0);
  TH1* hTrackDeltaPt = new TH1D("hTrackDeltaPt","-1.0<#eta<1.0;(Pt^{rec}-Pt^{gen})/Pt^{gen};Number of tracks for #pi+",100, -0.1, 0.1);
  TH1* hTrackDeltaP = new TH1D("hTrackDeltaP","-1.0<#eta<1.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks for #pi+",100, -0.1, 0.1);
  TH1* hTrackMissPIDRec = new TH1D("hTrackMissPIDRec","-1.0<#eta<1.0;PID;Number of tracks for #pi+",300,-300.,300.);
  TH1* hTrackGoodPIDRec = new TH1D("hTrackGoodPIDRec","-1.0<#eta<1.0;PID;Number of tracks for #pi+",300,-300.,300.);
    
  TH1* hEFlowTrackDeltaEta = new TH1D("hEFlowTrackDeltaEta","-1.0<#eta<1.0;#eta^{rec}-#eta^{gen};Number of eflowtracks for #pi+",100, -1.0, 1.0);
  TH1* hEFlowTrackDeltaPhi = new TH1D("hEFlowTrackDeltaPhi","-1.0<#eta<1.0;#phi^{rec}-#phi^{gen};Number of eflowtracks for #pi+",100, -1.0, 1.0);
  TH1* hEFlowTrackDeltaPt = new TH1D("hEFlowTrackDeltaPt","-1.0<#eta<1.0;(Pt^{rec}-Pt^{gen})/Pt^{gen};Number of eflowtracks for #pi+",100, -0.1, 0.1);
  TH1* hEFlowTrackDeltaP = new TH1D("hEFlowTrackDeltaP","-1.0<#eta<1.0;(P^{rec}-P^{gen})/P^{gen};Number of eflowtracks for #pi+",100, -0.1, 0.1);
    
  TH1* hTowerDeltaEta = new TH1D("hTowerDeltaEta","-1.0<#eta<1.0;#eta^{rec}-#eta^{gen};Number of towers for #pi+",100, -1.0, 1.0);
  TH1* hTowerDeltaPhi = new TH1D("hTowerDeltaPhi","-1.0<#eta<1.0;#phi^{rec}-#phi^{gen};Number of towers for #pi+",100, -1.0, 1.0);
  TH1* hTowerDeltaEem = new TH1D("hTowerDeltaEem","-1.0<#eta<1.0;Eem^{rec}/E^{gen};Number of towers for #pi+",100, 0., 3.0);
  TH1* hTowerDeltaEhad = new TH1D("hTowerDeltaEhad","-1.0<#eta<1.0;Ehad^{rec}/E^{gen};Number of towers for #pi+",100, 0., 3.0);
  TH1* hTowerDeltaE = new TH1D("hTowerDeltaE","-1.0<#eta<1.0;(E^{rec}-E^{gen})/E^{gen};Number of towers for #pi+",100, -0.5, 0.5);
  TH1* hTowerN = new TH1D("hTowerN", ";Number of Towers; Events",100,0.,100.);
  TH1* hTowerGenN = new TH1D("hTowerGenN", ";Number of Generated Particles associated with Towers; Events",100,0.,100.);
  TH1* hTowerGenPID = new TH1D("hTowerGenPID", ";PID of Generated Particles associated with Towers; Events",300,-300.,300.);
 
  Long64_t entry;
  Int_t i;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    //
    if (branchTrack->GetEntries() > 0)
    {
      // Loop on all tracks
      for (i = 0; i < branchTrack->GetEntries(); ++i)
      {
        track = (Track*) branchTrack->At(i);
        particle = (GenParticle*) track->Particle.GetObject();

	//
	if ((particle->PID == 211 && track->PID != 211) || (particle->PID != 211 && track->PID == 211))
	  hTrackMissPIDRec->Fill(track->PID);
        if (particle->PID == 211 && track->PID == 211)
          hTrackGoodPIDRec->Fill(track->PID);

	if (particle->PID == 211)
	{
	  hTrackDeltaEta->Fill((track->Eta - particle->Eta));
	  hTrackDeltaPhi->Fill((track->Phi - particle->Phi));
	  hTrackDeltaPt->Fill((track->PT - particle->PT)/particle->PT);
	  hTrackDeltaP->Fill((track->P - particle->P)/particle->P);
	}
      }
    }
    //
    if (branchEFlowTrack->GetEntries() > 0)
    {
      for (i = 0; i < branchEFlowTrack->GetEntries(); ++i)
      {
        eflowtrack = (Track*) branchEFlowTrack->At(i);
        particle = (GenParticle*) eflowtrack->Particle.GetObject();

        if (particle->PID == 211)
        {
          hEFlowTrackDeltaEta->Fill((eflowtrack->Eta - particle->Eta));
          hEFlowTrackDeltaPhi->Fill((eflowtrack->Phi - particle->Phi));
          hEFlowTrackDeltaPt->Fill((eflowtrack->PT - particle->PT)/particle->PT);
          hEFlowTrackDeltaP->Fill((eflowtrack->P - particle->P)/particle->P);
        }
      }
    }
    // 
    if (branchTower->GetEntries() > 0)
    {
      hTowerN->Fill(branchTower->GetEntriesFast());
      for (i = 0; i < branchTower->GetEntriesFast(); ++i)
      {
        tower = (Tower*) branchTower->At(i);
        hTowerGenN->Fill(tower->Particles.GetEntriesFast());
	for (int j=0; j< tower->Particles.GetEntriesFast(); ++j)
	{
          particle = (GenParticle*) tower->Particles.At(j);
          hTowerGenPID->Fill(particle->PID);

          if (particle->PID == 211)
          {
            hTowerDeltaEta->Fill((tower->Eta - particle->Eta));
            hTowerDeltaPhi->Fill((tower->Phi - particle->Phi));
            hTowerDeltaEem->Fill((tower->Eem/particle->E));
            hTowerDeltaEhad->Fill((tower->Ehad/particle->E));
            hTowerDeltaE->Fill((tower->E - particle->E)/particle->E);
	  }
        }
      }
    }
  }

  // Plot Track Figures
  TCanvas* c1 = new TCanvas("c1","c1",2000,1000);
  c1->Divide(4,2);
  c1->cd(1);
  hTrackDeltaEta->GetXaxis()->CenterTitle(true); 
  hTrackDeltaEta->GetYaxis()->CenterTitle(true); 
  hTrackDeltaEta->SetLineWidth(2);
  hTrackDeltaEta->Fit("gaus");
  TF1* gaus_hTrackDeltaEta = hTrackDeltaEta->GetFunction("gaus");
  gaus_hTrackDeltaEta->SetLineWidth(2);
  gaus_hTrackDeltaEta->SetLineColor(kRed);
  hTrackDeltaEta->Draw();
  //c1->SaveAs("./plots/hTrackDeltaEta.png");
  //TCanvas* c2 = new TCanvas("c2","c2");
  c1->cd(2);
  hTrackDeltaPhi->GetXaxis()->CenterTitle(true); 
  hTrackDeltaPhi->GetYaxis()->CenterTitle(true); 
  hTrackDeltaPhi->SetLineWidth(2);
  hTrackDeltaPhi->Fit("gaus");
  TF1* gaus_hTrackDeltaPhi = hTrackDeltaPhi->GetFunction("gaus");
  gaus_hTrackDeltaPhi->SetLineWidth(2);
  gaus_hTrackDeltaPhi->SetLineColor(kRed);  
  hTrackDeltaPhi->Draw();
  //c2->SaveAs("./plots/hTrackDeltaPhi.png");
  //TCanvas* c3 = new TCanvas("c3","c3");
  c1->cd(3);
  hTrackDeltaPt->GetXaxis()->CenterTitle(true); 
  hTrackDeltaPt->GetYaxis()->CenterTitle(true); 
  hTrackDeltaPt->SetLineWidth(2);
  hTrackDeltaPt->Fit("gaus");
  TF1* gaus_hTrackDeltaPt = hTrackDeltaPt->GetFunction("gaus");
  gaus_hTrackDeltaPt->SetLineWidth(2);
  gaus_hTrackDeltaPt->SetLineColor(kRed);
  hTrackDeltaPt->Draw();
  //c3->SaveAs("./plots/hTrackDeltaPt.png");
  //TCanvas* c4 = new TCanvas("c4","c4");
  c1->cd(4);
  hTrackDeltaP->GetXaxis()->CenterTitle(true); 
  hTrackDeltaP->GetYaxis()->CenterTitle(true); 
  hTrackDeltaP->SetLineWidth(2);
  hTrackDeltaP->Fit("gaus");
  TF1* gaus_hTrackDeltaP = hTrackDeltaP->GetFunction("gaus");
  gaus_hTrackDeltaP->SetLineWidth(2);
  gaus_hTrackDeltaP->SetLineColor(kRed);
  hTrackDeltaP->Draw();
  //c4->SaveAs("./plots/hTrackDeltaP.png");
  c1->cd(5);
  hTrackMissPIDRec->GetXaxis()->CenterTitle(true); 
  hTrackMissPIDRec->GetYaxis()->CenterTitle(true); 
  hTrackMissPIDRec->SetLineWidth(2);
  hTrackMissPIDRec->Draw();
  c1->cd(6);
  hTrackGoodPIDRec->GetXaxis()->CenterTitle(true); 
  hTrackGoodPIDRec->GetYaxis()->CenterTitle(true); 
  hTrackGoodPIDRec->SetLineWidth(2);
  hTrackGoodPIDRec->Draw();
  c1->SaveAs("./plots/hTrackSummary.png");

  TCanvas* c3 = new TCanvas("c3","c3",2000,500);
  c3->Divide(4,1);
  c3->cd(1);
  //TCanvas* c9 = new TCanvas("c9","c9");
  hEFlowTrackDeltaEta->GetXaxis()->CenterTitle(true); 
  hEFlowTrackDeltaEta->GetYaxis()->CenterTitle(true); 
  hEFlowTrackDeltaEta->SetLineWidth(2);
  hEFlowTrackDeltaEta->Fit("gaus");
  TF1* gaus_hEFlowDeltaEta = hEFlowTrackDeltaEta->GetFunction("gaus");
  gaus_hEFlowDeltaEta->SetLineWidth(2);
  gaus_hEFlowDeltaEta->SetLineColor(kRed);  
  hEFlowTrackDeltaEta->Draw();
  //c9->SaveAs("./plots/hEFlowTrackDeltaEta.png");
  c3->cd(2);
  //TCanvas* c10 = new TCanvas("c10","c10");
  hEFlowTrackDeltaPhi->GetXaxis()->CenterTitle(true); 
  hEFlowTrackDeltaPhi->GetYaxis()->CenterTitle(true); 
  hEFlowTrackDeltaPhi->SetLineWidth(2);
  hEFlowTrackDeltaPhi->Fit("gaus");
  TF1* gaus_hEFlowTrackDeltaPhi = hEFlowTrackDeltaPhi->GetFunction("gaus");
  gaus_hEFlowTrackDeltaPhi->SetLineWidth(2);
  gaus_hEFlowTrackDeltaPhi->SetLineColor(kRed);  
  hEFlowTrackDeltaPhi->Draw();
  //c10->SaveAs("./plots/hEFlowTrackDeltaPhi.png");
  c3->cd(3);
  //TCanvas* c11 = new TCanvas("c11","c11");
  hEFlowTrackDeltaPt->GetXaxis()->CenterTitle(true); 
  hEFlowTrackDeltaPt->GetYaxis()->CenterTitle(true); 
  hEFlowTrackDeltaPt->SetLineWidth(2);
  hEFlowTrackDeltaPt->Fit("gaus");
  TF1* gaus_hEFlowTrackDeltaPt = hEFlowTrackDeltaPt->GetFunction("gaus");
  gaus_hEFlowTrackDeltaPt->SetLineWidth(2);
  gaus_hEFlowTrackDeltaPt->SetLineColor(kRed);
  hEFlowTrackDeltaPt->Draw();
  //c11->SaveAs("./plots/hEFlowTrackDeltaPt.png");
  c3->cd(4);
  //TCanvas* c12 = new TCanvas("c12","c12");
  hEFlowTrackDeltaP->GetXaxis()->CenterTitle(true); 
  hEFlowTrackDeltaP->GetYaxis()->CenterTitle(true); 
  hEFlowTrackDeltaP->SetLineWidth(2);
  hEFlowTrackDeltaP->Fit("gaus");
  TF1* gaus_hEFlowTrackDeltaP = hEFlowTrackDeltaP->GetFunction("gaus");
  gaus_hEFlowTrackDeltaP->SetLineWidth(2);
  gaus_hEFlowTrackDeltaP->SetLineColor(kRed);
  hEFlowTrackDeltaP->Draw();
  //c12->SaveAs("./plots/hEFlowTrackDeltaP.png");
  c3->SaveAs("./plots/hEFlowTrackSummary.png");

  TCanvas* c4 = new TCanvas("c4","c4",2000,1000);
  c4->Divide(4,2);
  c4->cd(1);
  //TCanvas* c13 = new TCanvas("c13","c13");
  hTowerDeltaEta->GetXaxis()->CenterTitle(true); 
  hTowerDeltaEta->GetYaxis()->CenterTitle(true); 
  hTowerDeltaEta->SetLineWidth(2);
  hTowerDeltaEta->Fit("gaus");
  TF1* gaus_hTowerDeltaEta = hTowerDeltaEta->GetFunction("gaus");
  gaus_hTowerDeltaEta->SetLineWidth(2);
  gaus_hTowerDeltaEta->SetLineColor(kRed);
  hTowerDeltaEta->Draw();
  //c13->SaveAs("./plots/hTowerDeltaEta.png");
  c4->cd(2);
  //TCanvas* c14 = new TCanvas("c14","c14");
  hTowerDeltaPhi->GetXaxis()->CenterTitle(true); 
  hTowerDeltaPhi->GetYaxis()->CenterTitle(true); 
  hTowerDeltaPhi->SetLineWidth(2);
  hTowerDeltaPhi->Fit("gaus");
  TF1* gaus_hTowerDeltaPhi = hTowerDeltaPhi->GetFunction("gaus");
  gaus_hTowerDeltaPhi->SetLineWidth(2);
  gaus_hTowerDeltaPhi->SetLineColor(kRed);
  hTowerDeltaPhi->Draw();
  //c14->SaveAs("./plots/hTowerDeltaPhi.png");
  c4->cd(3);
  //TCanvas* c15 = new TCanvas("c15","c15");
  hTowerDeltaEem->GetXaxis()->CenterTitle(true); 
  hTowerDeltaEem->GetYaxis()->CenterTitle(true); 
  hTowerDeltaEem->SetLineWidth(2);
  hTowerDeltaEem->Draw();
  //c15->SaveAs("./plots/hTowerDeltaEem.png");
  c4->cd(4);
  //TCanvas* c16 = new TCanvas("c16","c16");
  hTowerDeltaEhad->GetXaxis()->CenterTitle(true); 
  hTowerDeltaEhad->GetYaxis()->CenterTitle(true); 
  hTowerDeltaEhad->SetLineWidth(2);
  hTowerDeltaEhad->Draw();
  //c16->SaveAs("./plots/hTowerDeltaEhad.png");
  c4->cd(5);
  //TCanvas* c17 = new TCanvas("c17","c17");
  hTowerDeltaE->GetXaxis()->CenterTitle(true); 
  hTowerDeltaE->GetYaxis()->CenterTitle(true); 
  hTowerDeltaE->SetLineWidth(2);
  hTowerDeltaE->Fit("gaus");
  TF1* gaus_hTowerDeltaE = hTowerDeltaE->GetFunction("gaus");
  gaus_hTowerDeltaE->SetLineWidth(2);
  gaus_hTowerDeltaE->SetLineColor(kRed);
  hTowerDeltaE->Draw();
  //c17->SaveAs("./plots/hTowerDeltaE.png");
  c4->cd(6);
  //TCanvas *cnv18 = new TCanvas("cnv18", "cnv18");
  hTowerN->GetXaxis()->CenterTitle(true);
  hTowerN->GetYaxis()->CenterTitle(true);
  hTowerN->SetLineWidth(2);
  hTowerN->Draw();
  //cnv18->SaveAs("./plots/hTowerN.png");
  c4->cd(7);
  //TCanvas *cnv19 = new TCanvas("cnv19", "cnv19");
  hTowerGenN->GetXaxis()->CenterTitle(true);
  hTowerGenN->GetYaxis()->CenterTitle(true);
  hTowerGenN->SetLineWidth(2);
  hTowerGenN->Draw();
  //cnv19->SaveAs("./plots/hTowerGenN.png");
  c4->cd(8);
  //TCanvas *cnv20 = new TCanvas("cnv20", "cnv20");
  hTowerGenPID->GetXaxis()->CenterTitle(true);
  hTowerGenPID->GetYaxis()->CenterTitle(true);
  hTowerGenPID->SetLineWidth(2);
  hTowerGenPID->Draw();
  //cnv20->SaveAs("./plots/hTowerGenPID.png");
  c4->SaveAs("./plots/hTowerSummary.png");

  cout << "** Done..." << endl;

  delete treeReader;
  delete chain;
}
