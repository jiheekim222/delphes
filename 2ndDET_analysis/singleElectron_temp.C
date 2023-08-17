/*
This macro shows how to access the particle-level reference for reconstructed objects.
It is also shown how to loop over the jet constituents.

root -l examples/Example3.C'("delphes_output.root")'
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


//------------------------------------------------------------------------------  
// constant and binning

const int nEtaBin=10;
float eta[nEtaBin+1]={-5,-4,-3,-2,-1,0,1,2,3,4,5};

const int nEBin=8;
float E[nEBin+1]={0,1,2,3,5,7,10,15,20};

const int nPtBin=8;
float Pt[nPtBin+1]={0,1,2,3,5,7,10,15,20};
//------------------------------------------------------------------------------ 

int GetEtaBin(float myEta);
int GetEBin(float myE);
int GetPtBin(float myPt);

//------------------------------------------------------------------------------

struct TestPlots
{
  TH1 *hTrackDeltaEta[nEtaBin];
  TH1 *hTrackDeltaPt[nEtaBin][nPtBin];
  TH1 *hElectronDeltaEta[nEtaBin];
  TH1 *hElectronDeltaPt[nEtaBin][nPtBin];
};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  int i,j,k;
  
  for (i=0;i<nEtaBin;i++) 
  {
    plots->hTrackDeltaEta[i] =
      result->AddHist1D(Form("track_delta_eta_eta%d",i), 
                        Form("%0.f #leq #eta < %.0f",eta[i],eta[i+1]),
                        "(#eta^{gen} - #eta^{rec})/#eta^{gen}", "number of track for electrons",
                        100, -0.1, 0.1);

    plots->hElectronDeltaEta[i] =
      result->AddHist1D(Form("electron_delta_eta_eta%d",i), 
			Form("%0.f #leq #eta < %.0f",eta[i],eta[i+1]),
			"(#eta^{rec} - #eta^{gen})/#eta^{gen}", "number of electrons",
			100, -0.1, 0.1);
  
    for (j=0;j<nPtBin;j++) {
      plots->hTrackDeltaPt[i][j] = 
        result->AddHist1D(Form("track_delta_transverse_momentum_eta%dE%d",i,j),
			  Form("%0.f #leq #eta < %.0f, %0.f - %0.f GeV",eta[i],eta[i+1],Pt[j],Pt[j+1]),
			  "(Pt^{rec} - Pt^{gen})/Pt^{gen}", "number of track for electrons",
			  100, -0.1, 0.1);
    }

    for (j=0;j<nPtBin;j++) {
      plots->hElectronDeltaPt[i][j] = 
        result->AddHist1D(Form("electron_delta_transverse_momentum_eta%dE%d",i,j),
                          Form("%0.f #leq #eta < %.0f, %0.f - %0.f GeV",eta[i],eta[i+1],Pt[j],Pt[j+1]),
                          "(Pt^{rec} - Pt^{gen})/Pt^{gen}", "number of electrons",
                          100, -0.1, 0.1);
    }
  }
}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchTrack = treeReader->UseBranch("Track");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Electron *electron;

  Track *track;
  Tower *tower;

  TObject *object;

  TLorentzVector momentum;

  Long64_t entry;

  Int_t i, j, pdgCode;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);

    // Loop over all tracks in event
    for (i = 0; i < branchTrack->GetEntriesFast(); ++i)
    {
      track = (Track*) branchTrack->At(i);
      if (track->PT <= 0) continue;
      particle = (GenParticle*) track->Particle.GetObject();
      
      int iPt = GetPtBin(particle->PT);
      int iEta = GetEtaBin(particle->Eta);
      if (iPt<0 || iEta<0) continue;

      if (particle->PID == 11)
      {
        plots->hTrackDeltaEta[iEta]->Fill((track->Eta - particle->Eta)/particle->Eta);
        plots->hTrackDeltaPt[iEta][iPt]->Fill((track->PT - particle->PT)/particle->PT);
      } 
    }

    // Loop over all electrons in event
    for(i = 0; i < branchElectron->GetEntriesFast(); ++i)
    {
      electron = (Electron*) branchElectron->At(i);
      if (electron->PT <= 0) continue;
      particle = (GenParticle*) electron->Particle.GetObject();
      
      int iPt = GetPtBin(particle->PT);
      int iEta = GetEtaBin(particle->Eta);
      if (iPt<0 || iEta<0) continue;

      plots->hElectronDeltaEta[iEta]->Fill((electron->Eta - particle->Eta)/particle->Eta);
      plots->hElectronDeltaPt[iEta][iPt]->Fill((electron->PT - particle->PT)/particle->PT);
    }
  }
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, TestPlots *plots)
{
  result->Print("png");
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

int GetPtBin(float myPt)
{
  int i;
  
  if (myPt<=0) return -1;

  for (i=0;i<nPtBin;i++) {
    if (myPt>=Pt[i] && myPt<Pt[i+1]) return i;
  }
  
  return -1;
}

//------------------------------------------------------------------------------



void singleElectron(const char *inputFile)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  TestPlots *plots = new TestPlots;

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots);

  PrintHistograms(result, plots);

  result->Write("results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}
