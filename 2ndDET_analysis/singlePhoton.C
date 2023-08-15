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


//------------------------------------------------------------------------------  
// constant and binning

const int nEtaBin=4;
float eta[nEtaBin+1]={-1,-0.5,0,0.5,1};

const int nEBin=5;
float E[nEBin+1]={5,6,7,8,9,10};
//------------------------------------------------------------------------------ 

int GetEtaBin(float myEta);
int GetEBin(float myE);

//------------------------------------------------------------------------------

struct TestPlots
{
  TH1 *fPhotonDeltaEta[nEtaBin];
  TH1 *fPhotonDeltaE[nEtaBin][nEBin];

  TH1 *fPhotonEta;
  TH1 *fParticleEta;
};

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, TestPlots *plots)
{
  // AddHisto1D(name,title,x-axis label,y-axis label,nBin,range1, range2);

  int i,j;
  
  for (i=0;i<nEtaBin;i++) {
    plots->fPhotonDeltaEta[i] =
      result->AddHist1D(Form("photon_delta_eta_eta%d",i), 
			Form("%.1f #leq #eta < %.1f",eta[i],eta[i+1]),
			"#eta^{gen} - #eta^{rec}", "number of photons",
			//"(#eta^{particle} - #eta^{photon})/#eta^{particle}", "number of photons",
			100, -0.25, 0.25);
    
    for (j=0;j<nEBin;j++) {
      plots->fPhotonDeltaE[i][j] = 
	result->AddHist1D(Form("photon_delta_energy_eta%dE%d",i,j),
			  Form("%.1f #leq #eta < %.1f, %.0f - %.0f GeV",eta[i],eta[i+1],E[j],E[j+1]),
			  "(E^{gen} - E^{rec})/E^{gen}", "number of photons",
			  100, -0.2, 0.2);
    }
  }

  plots->fPhotonEta = result->AddHist1D("hPhotonEta","photon",";#eta_{rec};","",88, -1.1, 1.1);
  plots->fParticleEta = result->AddHist1D("hParticleEta","particle",";#eta_{gen};","",88,-1.1,1.1);
}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, TestPlots *plots)
{
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
  TClonesArray *branchEFlowPhoton = treeReader->UseBranch("EFlowPhoton");

  Long64_t allEntries = treeReader->GetEntries();

  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle, *myPart1, *myPart2;
  Photon *photon;

  Track *track;
  Tower *tower;

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
  
    photon = (Photon*) branchPhoton->At(0);
    
    for(i = 0; i < branchPhoton->GetEntriesFast(); ++i)
    {
      photon = (Photon*) branchPhoton->At(i);
      // skip photons with references to multiple particles
      if(photon->Particles.GetEntriesFast() != 1) continue;
      
      particle = (GenParticle*) photon->Particles.At(0);
      if (!particle) continue;
      if (particle->PID!=22) continue;
      
      int iE = GetEBin(particle->E);
      int iEta = GetEtaBin(particle->Eta);
      if (iE<0 || iEta<0) continue;
      
      plots->fPhotonEta->Fill(photon->Eta);
      plots->fParticleEta->Fill(particle->Eta);
      
      /*
	if (fabs((particle->Eta - photon->Eta)/particle->Eta) >0.1) {
	cout<<"gen eta = "<<particle->Eta<<", rec eta="<<photon->Eta
	<<", resolution="<<(particle->Eta - photon->Eta)/particle->Eta<<endl;
	}
      */
      //cout<<"(Egen-Erec)/Egen="<<"("<<particle->E<<"-"<<photon->E<<")/"<<particle->E<<"="<<(particle->E - photon->E)/particle->E<<endl;
    
      plots->fPhotonDeltaEta[iEta]->Fill((particle->Eta - photon->Eta));///particle->Eta);
      plots->fPhotonDeltaE[iEta][iE]->Fill((particle->E - photon->E)/particle->E);

      // cout << "--  New event -- " << endl;
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

void singlePhoton(const char *inputFile)
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

  PrintHistograms(result, plots);

  result->Write("results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}
