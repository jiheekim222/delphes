/*
  08-15-2023, Cheuk-Ping Wong
  ---------------------------
  A simple macro to study photon reconstruction in Delphes.
  This macro output a plots of energy and eta resolutions of 
  photon from the calorimeters. 
  It also output a root file called results.root
  
  root -l singlePhoton.C'("delphes_output.root")'
*/

#ifdef __CLING__
R__LOAD_LIBRARY(libDelphes)
#include "classes/DelphesClasses.h"
#include "external/ExRootAnalysis/ExRootTreeReader.h"
#include "external/ExRootAnalysis/ExRootResult.h"
#else
class ExRootTreeReader;
class ExRootResult;
#endif

#include <math.h>

//------------------------------------------------------------------------------  
// constant and binning
const float PI=M_PI;

const int nEtaBin=4;
float eta[nEtaBin+1]={-1,-0.5,0,0.5,1};

const int nEBin=17;
float E[nEBin+1]={0,0.2,0.4,0.6,0.8,1,2,3,4,5,6,7,8,9,10,13,16,20};

const int nPBin=17;
float p[nPBin+1]={0,0.2,0.4,0.6,0.8,1,2,3,4,5,6,7,8,9,10,13,16,20};

const int nPtBin=17;   
float pt[nPtBin+1]={0,0.2,0.4,0.6,0.8,1,2,3,4,5,6,7,8,9,10,13,16,20};

//------------------------------------------------------------------------------ 

int GetEtaBin(float myEta);
int GetEBin(float myE);
int GetPBin(float myP);
int GetPtBin(float myPt);
//------------------------------------------------------------------------------

struct TestPlots
{
  //track branch
  TH1 *fTrkEta;
  TH1 *fTrkP;
  TH1 *fTrkPt;
  TH1 *fTrkOuterPhi;
  TH1 *fTrkOuterEta;
  TH1 *fTrkDeltaEta[nEtaBin];
  TH1 *fTrkDeltaP[nEtaBin][nPBin];
  TH1 *fTrkDeltaPt[nEtaBin][nPtBin];
  TH1 *fTrkDeltaDCAt[nEtaBin][nPtBin];
  TH1 *fTrkDeltaDCAz[nEtaBin][nPtBin];

  //tower branch
  TH1 *fTowEta;
  TH1 *fTowPhi;
  TH1 *fTowE;
  TH1 *fTowEhad;
  TH1 *fTowEem;
  TH2 *fTowEfracHad;
  TH2 *fTowEfracEM;

  TH1 *fTowDeltaEta[nEtaBin];
  TH1 *fTowDeltaE[nEtaBin][nEBin];
  
  TH1 *fParticleEta;
  TH1 *fParticleE[nEtaBin];
  TH1 *fParticleP[nEtaBin];
  TH1 *fParticlePt[nEtaBin];

};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  // AddHisto1D(name,title,x-axis label,y-axis label,nBin,range1, range2);

  int i,j;
  
  //---------------
  // track branch
  //---------------
  plots->fTrkEta = result->AddHist1D("hEtaTrk","eta - track branch","#eta","",88,-1.1,1.1);
  plots->fTrkP = result->AddHist1D("hPTrk","p - track branch","p (GeV)","",55,0.,11.);
  plots->fTrkPt = result->AddHist1D("hPtTrk","pT - track branch","p_{T} (GeV)","",55,0.,11.); 
  plots->fTrkOuterPhi = result->AddHist1D("hOuterPhiTrk","outer phi - track branch","outer #phi","",360,-PI,PI);
  plots->fTrkOuterEta = result->AddHist1D("hOuterEtaTrk","outer eta - track branch","outer #eta","",20,-5.,5.);
  
  plots->fTrkEta->Sumw2();
  plots->fTrkP->Sumw2();
  plots->fTrkPt->Sumw2();
  plots->fTrkOuterPhi->Sumw2();
  plots->fTrkOuterEta->Sumw2();

  for (i=0;i<nEtaBin;i++) {
    plots->fTrkDeltaEta[i] = 
      result->AddHist1D(Form("Trk_delta_eta_eta%d",i),
                        Form("%.1f #leq #eta < %.1f",eta[i],eta[i+1]),
                        "#eta^{gen} - #eta^{rec}", "number of trk",
			20, -0.01, 0.01);
    plots->fTrkDeltaEta[i]->Sumw2();

    for (j=0;j<nPBin;j++) {
      plots->fTrkDeltaP[i][j] =
	result->AddHist1D(Form("Trk_delta_p_eta%d_p%d",i,j),
			  Form("%.1f #leq #eta < %.1f, %.0f - %.0f GeV",eta[i],eta[i+1],p[j],p[j+1]),
			  "(p^{gen} - p^{rec})/p^{gen}", "number of trk",
			  55, -0.055, 0.055);
      plots->fTrkDeltaP[i][j]->Sumw2();
    }
    
    for (j=0;j<nPtBin;j++) {
      plots->fTrkDeltaPt[i][j] = 
	result->AddHist1D(Form("Trk_delta_pt_eta%d_pt%d",i,j),
			  Form("%.1f #leq #eta < %.1f, %.0f - %.0f GeV",eta[i],eta[i+1],pt[j],pt[j+1]),
			  "(p_{T}^{gen} - p_{T}^{rec})/p_{T}^{gen}", "number of trk", 
			  55, -0.055, 0.055);  
      plots->fTrkDeltaPt[i][j]->Sumw2();

      plots->fTrkDeltaDCAt[i][j] =
        result->AddHist1D(Form("Trk_delta_DCAt_eta%d_pt%d",i,j),
                          Form("%.1f #leq #eta < %.1f, %.0f - %.0f GeV",eta[i],eta[i+1],pt[j],pt[j+1]),
                          "(DCA_{T}^{gen} - DCA_{T}^{rec})/DCA_{T}^{gen}", "number of trk",
                          80, -0.1, 0.1);
      plots->fTrkDeltaDCAt[i][j]->Sumw2();

      plots->fTrkDeltaDCAz[i][j] =
        result->AddHist1D(Form("Trk_delta_DCAz_eta%d_pt%d",i,j),
                          Form("%.1f #leq #eta < %.1f, %.0f - %.0f GeV",eta[i],eta[i+1],pt[j],pt[j+1]),
                          "(DCA_{z}^{gen} - DCA_{z}^{rec})/DCA_{z}^{gen}", "number of trk",
                          80, -0.1, 0.1);
      plots->fTrkDeltaDCAz[i][j]->Sumw2();
    }
  }
  
  //---------------
  // tower branch
  //---------------
  plots->fTowEta = result->AddHist1D("hEtaTow","eta - tower branch","#eta","",88,-1.1,1.1);
  plots->fTowPhi = result->AddHist1D("hPhiTow","phi - tower branch","#phi","",360,-PI,PI);
  plots->fTowE = result->AddHist1D("hETow","E - tower branch","E (GeV)","",55,0.,11.);
  plots->fTowEhad = result->AddHist1D("hEhadTow","E (HCal) - tower branch","E^{had} (GeV)","",55,0.,11.);
  plots->fTowEem = result->AddHist1D("hEemTow","E (EMCal) - tower branch","E^{em} (GeV)","",55,0.,11.);
  plots->fTowEfracHad = result->AddHist2D("hEfracHadTow","E fraction (HCal) - tower branch",
					  "E^{gen} (GeV)","fraction",
					  35,3.,10.,
					  20,0.,1.);
  plots->fTowEfracEM = result->AddHist2D("hEfracEMTow","E fraction (EMCal) - tower branch",
					 "E^{gen} (GeV)","fraction",
					 35,3.,10.,
					 20,0.,1.);

  plots->fTowEta->Sumw2();
  plots->fTowPhi->Sumw2();
  plots->fTowE->Sumw2();
  plots->fTowEhad->Sumw2();
  plots->fTowEem->Sumw2();
  plots->fTowEfracHad->Sumw2();
  plots->fTowEfracEM->Sumw2();

  for (i=0;i<nEtaBin;i++) {
    plots->fTowDeltaEta[i] =
      result->AddHist1D(Form("Tow_delta_eta_eta%d",i),
                        Form("%.1f #leq #eta < %.1f",eta[i],eta[i+1]),
                        "#eta^{gen} - #eta^{rec}", "number of tow",
                        100, -0.25, 0.25);
    plots->fTowDeltaEta[i]->Sumw2();

    for (j=0;j<nEBin;j++) {
      plots->fTowDeltaE[i][j] = 
	result->AddHist1D(Form("Tow_delta_E_eta%d_E%d",i,j),
                          Form("%.1f #leq #eta < %.1f, %.0f - %.0f GeV",eta[i],eta[i+1],E[j],E[j+1]),
                          "(E^{gen} - E^{rec})/E^{gen}", "number of tow",
                          101, -1.05, 1.05);
      plots->fTowDeltaE[i][j]->Sumw2();
    }
  }
  
  //---------------
  // gen particles
  //---------------
  plots->fParticleEta = result->AddHist1D("hParticleEta","particle","#eta_{gen}","",88,-1.1,1.1);
  plots->fParticleEta->Sumw2();
  
  for (i=0;i<nEtaBin;i++) {
    //full mom
    plots->fParticleP[i] = result->AddHist1D(Form("hParticleP_eta%d",i), 
					     Form("%.1f #leq #eta < %.1f",eta[i],eta[i+1]),
					     "p^{gen} (GeV)",
					     "",
					     55,0.,11.);
    plots->fParticleP[i]->Sumw2(); 

    //pt
    plots->fParticlePt[i] = result->AddHist1D(Form("hParticlePt_eta%d",i),
					      Form("%.1f #leq #eta < %.1f",eta[i],eta[i+1]),
					      "p_{T}^{gen} (GeV)",
					      "",
					      55,0.,11.);
    plots->fParticlePt[i]->Sumw2();

    //E
    plots->fParticleE[i] = result->AddHist1D(Form("hParticleE_eta%d",i),
					     Form("%.1f #leq #eta < %.1f",eta[i],eta[i+1]),
					     "E^{gen} (GeV)",
					     "",
					     55,0.,11.);
    plots->fParticleE[i]->Sumw2();
  }
}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *part1,*part2;
  Track *trk;
  Tower *tow;

  TObject *object;

  TLorentzVector momentum;

  Float_t Eem, Ehad;
  Bool_t skip;

  Long64_t entry;

  Int_t i, j, pdgCode;

  // Loop over all events
  //for(entry = 0; entry < 100; ++entry)
  for(entry = 0; entry < allEntries; ++entry)
    {
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);
  
      for(i = 0; i < branchTower->GetEntriesFast(); i++)
	{
	  //cout<<"evt "<<entry<<", tower "<<i<<endl;
	  //tower size and trk size are the same in this case
	  tow = (Tower*) branchTower->At(i); 
	  trk = (Track*) branchTrack->At(i);

	  if (!tow) {
	    cout<<"no tow, next evt."<<endl;
            continue;
	  }
	  if (!trk) {
	    //cout<<"no trk, next evt."<<endl;
	    //part1 = (GenParticle*) tow->Particles.At(0);	    
	    //cout<<"tower pid="<<part1->PID<<", mom="<<part1->P<<endl;
	    continue;
	  }
	  
	  // skip photons with references to multiple particles
	  //if(photon->Particles.GetEntriesFast() != 1) continue;
      
	  part1 = (GenParticle*) tow->Particles.At(0);
	  part2 = (GenParticle*) trk->Particle.GetObject();

      
	  if (!part1 || !part2) continue;
	  if (part1->PID!=211 || part2->PID!=211) continue;
	  //if (part1->fUniqueID!=part->fUniqueID) continue;

	  int iE = GetEBin(part1->E);
	  int iP = GetEBin(part1->P);
	  int iPt =GetEBin(part1->PT);
	  int iEta = GetEtaBin(part1->Eta);
      
	  if (iEta<0) continue;

	  //---------------
	  //fill tower histo
	  //---------------
	  plots->fTowEta->Fill(tow->Eta);
	  plots->fTowPhi->Fill(tow->Phi);
	  plots->fTowE->Fill(tow->E);
	  plots->fTowEhad->Fill(tow->Ehad);
	  plots->fTowEem->Fill(tow->Eem);
	  plots->fTowEfracHad->Fill(part1->E,tow->Ehad/part1->E);
	  plots->fTowEfracEM->Fill(part1->E,tow->Eem/part1->E);
	  plots->fTowDeltaEta[iEta]->Fill(part1->Eta - tow->Eta);
	  if (iE>=0)plots->fTowDeltaE[iEta][iE]->Fill((part1->E - tow->E)/part1->E);

	  //---------------
	  //fill trk histo
	  //---------------
	  plots->fTrkEta->Fill(trk->Eta);
	  plots->fTrkP->Fill(trk->P);
	  plots->fTrkPt->Fill(trk->PT);
	  plots->fTrkOuterPhi->Fill(trk->PhiOuter);
	  plots->fTrkOuterEta->Fill(trk->EtaOuter);
	  plots->fTrkDeltaEta[iEta]->Fill(part1->Eta - trk->Eta);
	  if (iP>=0) plots->fTrkDeltaP[iEta][iP]->Fill((part1->P - trk->P)/part1->P);
	  if (iPt>=0) {
	    plots->fTrkDeltaPt[iEta][iPt]->Fill((part1->PT - trk->PT)/part1->PT);
	    plots->fTrkDeltaDCAt[iEta][iPt]->Fill(trk->D0);
	    plots->fTrkDeltaDCAz[iEta][iPt]->Fill(trk->DZ);
	  }

	  //---------------
	  //fill part1 histo
	  //---------------
	  plots->fParticleEta->Fill(part1->Eta);
	  plots->fParticleE[iEta]->Fill(part1->E);
	  plots->fParticleP[iEta]->Fill(part1->P);
	  plots->fParticlePt[iEta]->Fill(part1->PT);
	}
    }
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("pdf");
}
//------------------------------------------------------------------------------

int GetEtaBin(float myEta)
{
  int i;

  for (i=0;i<nEtaBin;i++) {
    if (myEta>=eta[i] && myEta<eta[i+1]) return i;
  }

  return -1;
}

//------------------------------------------------------------------------------ 

int GetEBin(float myE)
{
  int i;
  
  if (myE<=0) return -1;

  for (i=0;i<nEBin;i++) {
    if (myE>=E[i] && myE<E[i+1]) return i;
  }
  
  return -1;
}
//------------------------------------------------------------------------------ 

int GetPBin(float myP)
{
  int i;

  if (myP<=0) return -1;

  for (i=0;i<nPBin;i++) {
    if (myP>=p[i] && myP<p[i+1]) return i;
  }

  return -1;
}

//------------------------------------------------------------------------------ 

int GetPtBin(float myPt)
{
  int i;

  if (myPt<=0) return -1;

  for (i=0;i<nPtBin;i++) {
    if (myPt>=pt[i] && myPt<pt[i+1]) return i;
  }

  return -1;
}

//------------------------------------------------------------------------------

void singlePi(const char *inputFile)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  TestPlots *plots = new TestPlots;

  //-------------------------
  // Book histo
  //-------------------------
  

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots);
  
  //PrintHistograms(result, plots);

  result->Write("results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}
