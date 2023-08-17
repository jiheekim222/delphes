//////////////////////////////////////////
// 08/15/2023 Jihee Kim (jkim11@bnl.gov)
// Single Electron
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

void singleElectron(const char *inputFile)
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
  TClonesArray *branchElectron = treeReader->UseBranch("Electron");
  TClonesArray *branchTower = treeReader->UseBranch("Tower");
  TClonesArray *branchEFlowTrack = treeReader->UseBranch("EFlowTrack");

  // Get total number of events
  Long64_t allEntries = treeReader->GetEntries();
  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Track *track, *eflowtrack;
  Electron *electron;
  Tower *tower;
 
  // Book histograms
  const int nEtaBin=2;
  float eta[nEtaBin+1]={-1,0,1};
  
  const int nPBin=4;
  float p[nPBin+1]={3,5,7,10,15};
 
  const int nPtBin=4;
  float pt[nPtBin+1]={3,5,7,10,15};
  
  const int nEBin=4;
  float e[nEBin+1]={3,5,7,10,15};

  TH1D* hTrackDeltaEta[nEtaBin];
  TH1D* hTrackDeltaPt[nEtaBin][nPtBin];
  TH1D* hTrackDeltaP[nEtaBin][nPBin];
  TH1D* hElectronDeltaEta[nEtaBin];
  TH1D* hElectronDeltaPt[nEtaBin][nPtBin];
  TH1D* hElectronEhadOverEemPt[nEtaBin][nPtBin];
  TH1D* hEFlowTrackDeltaEta[nEtaBin];
  TH1D* hEFlowTrackDeltaPt[nEtaBin][nPtBin];
  TH1D* hEFlowTrackDeltaP[nEtaBin][nPBin];
  TH1D* hTowerDeltaEta[nEtaBin];
  TH1D* hTowerDeltaEem[nEtaBin];
  TH1D* hTowerDeltaEhad[nEtaBin];
  TH1D* hTowerDeltaE[nEtaBin][nEBin];
  for (int i = 0; i < nEtaBin; ++i)
  {
    hTrackDeltaEta[i] = new TH1D(Form("hTrackDeltaEta%d",i),Form("%.1f<#eta<%0.1f;#eta^{rec}-#eta^{gen};Number of tracks for electron",eta[i],eta[i+1]),100, -1.0, 1.0);
    hElectronDeltaEta[i] = new TH1D(Form("hElectronDeltaEta%d",i),Form("%.1f<#eta<%0.1f;#eta^{rec}-#eta^{gen};Number of electrons",eta[i],eta[i+1]),100, -1.0, 1.0);
    hEFlowTrackDeltaEta[i] = new TH1D(Form("hEFlowTrackDeltaEta%d",i),Form("%.1f<#eta<%0.1f;#eta^{rec}-#eta^{gen};Number of eflowtracks for electron",eta[i],eta[i+1]),100, -1.0, 1.0);
    hTowerDeltaEta[i] = new TH1D(Form("hTowerDeltaEta%d",i),Form("%.1f<#eta<%0.1f;#eta^{rec}-#eta^{gen};Number of towers for electron",eta[i],eta[i+1]),100, -1.0, 1.0);
    hTowerDeltaEem[i] = new TH1D(Form("hTowerDeltaEem%d",i),Form("%.1f<#eta<%0.1f;Eem^{rec}/E^{gen};Number of towers for electron",eta[i],eta[i+1]),100, -2., 2.0);
    hTowerDeltaEhad[i] = new TH1D(Form("hTowerDeltaEhad%d",i),Form("%.1f<#eta<%0.1f;Ehad^{rec}/E^{gen};Number of towers for electron",eta[i],eta[i+1]),100, -2., 2.0);
    for (int j = 0; j < nPtBin; ++j)
    {
      hTrackDeltaPt[i][j] = new TH1D(Form("hTrackDeltaEta%dPt%d",i,j),Form("%.1f<#eta<%0.1f, pt=%0.f-%0.fGeV;(Pt^{rec}-Pt^{gen})/Pt^{gen};Number of tracks for electron",eta[i],eta[i+1],pt[j],pt[j+1]),100, -0.1, 0.1);
      hElectronDeltaPt[i][j] = new TH1D(Form("hElectronDeltaEta%dPt%d",i,j),Form("%.1f<#eta<%0.1f, pt=%0.f-%0.fGeV;(Pt^{rec}-Pt^{gen})/Pt^{gen};Number of electrons",eta[i],eta[i+1],pt[j],pt[j+1]),100, -0.1, 0.1);
      hElectronEhadOverEemPt[i][j] = new TH1D(Form("hElectronEhadOverEemEta%dPt%d",i,j),Form("%.1f<#eta<%0.1f, pt=%0.f-%0.fGeV;Ratio of Ehad over Eem;Number of electrons",eta[i],eta[i+1],pt[j],pt[j+1]),100, -0.1, 1.0);
      hEFlowTrackDeltaPt[i][j] = new TH1D(Form("hEFlowTrackDeltaEta%dPt%d",i,j),Form("%.1f<#eta<%0.1f, pt=%0.f-%0.fGeV;(Pt^{rec}-Pt^{gen})/Pt^{gen};Number of eflowtracks for electron",eta[i],eta[i+1],pt[j],pt[j+1]),100, -0.1, 0.1);
    }
    for (int k = 0; k < nPBin; ++k)
    {
      hTrackDeltaP[i][k] = new TH1D(Form("hTrackDeltaEta%dP%d",i,k),Form("%.1f<#eta<%0.1f, p=%0.f-%0.fGeV;(P^{rec}-P^{gen})/P^{gen};Number of tracks for electron",eta[i],eta[i+1],p[k],p[k+1]),100, -0.1, 0.1);
      hEFlowTrackDeltaP[i][k] = new TH1D(Form("hEFlowTrackDeltaEta%dP%d",i,k),Form("%.1f<#eta<%0.1f, p=%0.f-%0.fGeV;(P^{rec}-P^{gen})/P^{gen};Number of eflowtracks for electron",eta[i],eta[i+1],p[k],p[k+1]),100, -0.1, 0.1);
    }
    for (int m = 0; m < nEBin; ++m)
    {
      hTowerDeltaE[i][m] = new TH1D(Form("hTowerDeltaEta%dE%d",i,m),Form("%.1f<#eta<%0.1f, E=%0.f-%0.fGeV;(E^{rec}-E^{gen})/E^{gen};Number of towers for electron",eta[i],eta[i+1],e[m],e[m+1]),100, -0.1, 0.1);
    }
  }
 
  TH1* hTowerN = new TH1D("hTowerN", ";Number of Towers; Events",100,0.,100.);
  TH1* hTowerGenN = new TH1D("hTowerGenN", ";Number of Generated Particles associated with Towers; Events",100,0.,100.);
  TH1* hTowerGenPID = new TH1D("hTowerGenPID", ";PID of Generated Particles associated with Towers; Events",100,-100.,100.);
 
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

	int iEta, iPt, iP;
        for (int i = 0; i < nEtaBin; ++i)
          if (particle->Eta >= eta[i] && particle->Eta < eta[i+1])
            iEta = i;
        for (int i = 0; i < nPtBin; ++i)
          if (particle->PT >= pt[i] && particle->PT < pt[i+1])
	    iPt = i;
        for (int i = 0; i < nPBin; ++i)
          if (particle->P >= p[i] && particle->P < p[i+1])
	    iP = i;

        if (iP< 0 || iPt < 0 || iEta < 0) continue;

	if (particle->PID == 11)
	{
	  hTrackDeltaEta[iEta]->Fill((track->Eta - particle->Eta));
	  hTrackDeltaPt[iEta][iPt]->Fill((track->PT - particle->PT)/particle->PT);
	  hTrackDeltaP[iEta][iP]->Fill((track->P - particle->P)/particle->P);
	  //cout << "Printing Track: eta " << track->Eta << endl;
          //cout << "Printing Track: eta gen " << particle->Eta << endl;
	  //cout << "Printing Track: eta outer" << track->EtaOuter << endl;
          //cout << "Printing Track: p " << track->P << endl;
          //cout << "Printing Track: p gen " << particle->P << endl;
        }
      }
    }
 
    if (branchElectron->GetEntries() > 0)
    {
      // Loop on all electrons
      for (i = 0; i < branchElectron->GetEntries(); ++i)
      {
	electron = (Electron*) branchElectron->At(i);
	particle = (GenParticle*) electron->Particle.GetObject();

        int iEta, iPt;
        for (int i = 0; i < nEtaBin; ++i)
          if (particle->Eta >= eta[i] && particle->Eta < eta[i+1])
            iEta = i;
        for (int i = 0; i < nPtBin; ++i)
          if (particle->PT >= pt[i] && particle->PT < pt[i+1])
            iPt = i;
        if (iPt < 0 || iEta < 0) continue;

        if (particle->PID == 11)
        {
          hElectronDeltaEta[iEta]->Fill((electron->Eta - particle->Eta));
          hElectronDeltaPt[iEta][iPt]->Fill((electron->PT - particle->PT)/particle->PT);
          hElectronEhadOverEemPt[iEta][iPt]->Fill(electron->EhadOverEem);
          //cout << "Printing Electron: eta " << electron->Eta << endl;
          //cout << "Printing Electron: eta gen " << particle->Eta << endl;
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

        int iEta, iPt, iP;
        for (int i = 0; i < nEtaBin; ++i)
          if (particle->Eta >= eta[i] && particle->Eta < eta[i+1])
            iEta = i;
        for (int i = 0; i < nPtBin; ++i)
          if (particle->PT >= pt[i] && particle->PT < pt[i+1])
            iPt = i;
        for (int i = 0; i < nPBin; ++i)
          if (particle->P >= p[i] && particle->P < p[i+1])
            iP = i;

        if (iP< 0 || iPt < 0 || iEta < 0) continue;

        if (particle->PID == 11)
        {
          hEFlowTrackDeltaEta[iEta]->Fill((eflowtrack->Eta - particle->Eta));
          hEFlowTrackDeltaPt[iEta][iPt]->Fill((eflowtrack->PT - particle->PT)/particle->PT);
          hEFlowTrackDeltaP[iEta][iP]->Fill((eflowtrack->P - particle->P)/particle->P);
          //cout << "Printing eflowtrack: eta " << eflowtrack->Eta << endl;
          //cout << "Printing eflowtrack: eta gen " << particle->Eta << endl;
          //cout << "Printing eflowtrack: eta outer" << eflowtrack->EtaOuter << endl;
          //cout << "Printing eflowtrack: p " << eflowtrack->P << endl;
          //cout << "Printing eflowtrack: p gen " << particle->P << endl;
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
          particle = (GenParticle*) tower->Particles.At(i);
          hTowerGenPID->Fill(particle->PID);
          int iEta, iE;
          for (int i = 0; i < nEtaBin; ++i)
            if (particle->Eta >= eta[i] && particle->Eta < eta[i+1])
              iEta = i;
          for (int i = 0; i < nEBin; ++i)
            if (particle->E >= e[i] && particle->E < e[i+1])
              iE = i;
          if (iE < 0 || iEta < 0) continue;

          if (particle->PID == 11)
          {
            hTowerDeltaEta[iEta]->Fill((tower->Eta - particle->Eta));
            hTowerDeltaEem[iEta]->Fill((tower->Eem/particle->E));
            hTowerDeltaEhad[iEta]->Fill((tower->Ehad/particle->E));
            hTowerDeltaE[iEta][iE]->Fill((tower->E - particle->E)/particle->E);
	  }
        }
        //cout << "Printing Tower: eta " << tower->Eta << endl;
        //cout << "Printing Tower: eta gen " << particle->Eta << endl;
        //cout << "Printing Tower: E " << tower->E << endl;
        //cout << "Printing Tower: Eem " << tower->Eem << endl;
        //cout << "Printing Tower: Ehad " << tower->Ehad << endl;
        //cout << "Printing Tower: E gen " << particle->E << endl;
        //for (int j =0; j<4; ++j)
          //cout << "Printing Tower: Edges: " << tower->Edges[j] << endl;
      }
    }
    //if (entry < 1) break;
  }

  const int n = 1e6;
  int icnv = 0;
  TCanvas* c[n];
  TF1* gaus_hTrackDeltaEta[nEtaBin];
  TF1* gaus_hTrackDeltaPt[nEtaBin][nPtBin];
  TF1* gaus_hTrackDeltaP[nEtaBin][nPBin];
  TF1* gaus_hElectronDeltaEta[nEtaBin];
  TF1* gaus_hElectronDeltaPt[nEtaBin][nPtBin];
  TF1* gaus_hEFlowTrackDeltaEta[nEtaBin];
  TF1* gaus_hEFlowTrackDeltaPt[nEtaBin][nPtBin];
  TF1* gaus_hEFlowTrackDeltaP[nEtaBin][nPBin];
  TF1* gaus_hTowerDeltaEta[nEtaBin];
  TF1* gaus_hTowerDeltaEem[nEtaBin];
  TF1* gaus_hTowerDeltaEhad[nEtaBin];
  TF1* gaus_hTowerDeltaE[nEtaBin][nEBin];
  // Plot Track Figures
  for (int i = 0; i < nEtaBin; ++i)
  {
    icnv++;
    c[icnv] = new TCanvas(Form("c%d",icnv),Form("c%d",icnv));
    hTrackDeltaEta[i]->GetXaxis()->CenterTitle(true); 
    hTrackDeltaEta[i]->GetYaxis()->CenterTitle(true); 
    hTrackDeltaEta[i]->SetLineWidth(2);
    hTrackDeltaEta[i]->Fit("gaus");
    gaus_hTrackDeltaEta[i] = hTrackDeltaEta[i]->GetFunction("gaus");
    gaus_hTrackDeltaEta[i]->SetLineWidth(2);
    gaus_hTrackDeltaEta[i]->SetLineColor(kRed);
    hTrackDeltaEta[i]->Draw();
    c[icnv]->SaveAs(Form("./plots/hTrackDeltaEta%d.png",i));

    for (int j =0; j < nPtBin; ++j)
    {
      icnv++;
      c[icnv] = new TCanvas(Form("c%d",icnv),Form("c%d",icnv));
      hTrackDeltaPt[i][j]->GetXaxis()->CenterTitle(true); 
      hTrackDeltaPt[i][j]->GetYaxis()->CenterTitle(true); 
      hTrackDeltaPt[i][j]->SetLineWidth(2);
      hTrackDeltaPt[i][j]->Fit("gaus");
      if (hTrackDeltaPt[i][j]->GetFunction("gaus"))
      {
        gaus_hTrackDeltaPt[i][j] = hTrackDeltaPt[i][j]->GetFunction("gaus");
        gaus_hTrackDeltaPt[i][j]->SetLineWidth(2);
        gaus_hTrackDeltaPt[i][j]->SetLineColor(kRed);
        hTrackDeltaPt[i][j]->Draw();
      }
      c[icnv]->SaveAs(Form("./plots/hTrackDeltaEta%dPt%d.png",i,j));
    }
  
    for (int k =0; k < nPBin; ++k)
    {
      icnv++;
      c[icnv] = new TCanvas(Form("c%d",icnv),Form("c%d",icnv));
      hTrackDeltaP[i][k]->GetXaxis()->CenterTitle(true); 
      hTrackDeltaP[i][k]->GetYaxis()->CenterTitle(true); 
      hTrackDeltaP[i][k]->SetLineWidth(2);
      hTrackDeltaP[i][k]->Fit("gaus");
      if (hTrackDeltaP[i][k]->GetFunction("gaus"))
      {
        gaus_hTrackDeltaP[i][k] = hTrackDeltaP[i][k]->GetFunction("gaus");
        gaus_hTrackDeltaP[i][k]->SetLineWidth(2);
        gaus_hTrackDeltaP[i][k]->SetLineColor(kRed);
        hTrackDeltaP[i][k]->Draw();
        c[icnv]->SaveAs(Form("./plots/hTrackDeltaEta%dP%d.png",i,k));
      }
    }
  }
  // Plot Electron Figures
  for (int i = 0; i < nEtaBin; ++i)
  {
    icnv++;
    c[icnv] = new TCanvas(Form("c%d",icnv),Form("c%d",icnv));
    hElectronDeltaEta[i]->GetXaxis()->CenterTitle(true); 
    hElectronDeltaEta[i]->GetYaxis()->CenterTitle(true); 
    hElectronDeltaEta[i]->SetLineWidth(2);
    hElectronDeltaEta[i]->Fit("gaus");
    gaus_hElectronDeltaEta[i] = hElectronDeltaEta[i]->GetFunction("gaus");
    gaus_hElectronDeltaEta[i]->SetLineWidth(2);
    gaus_hElectronDeltaEta[i]->SetLineColor(kRed);
    hElectronDeltaEta[i]->Draw();
    c[icnv]->SaveAs(Form("./plots/hTrackDeltaEta%d.png",i));

    for (int j =0; j < nPtBin; ++j)
    {
      icnv++;
      c[icnv] = new TCanvas(Form("c%d",icnv),Form("c%d",icnv));
      hElectronDeltaPt[i][j]->GetXaxis()->CenterTitle(true); 
      hElectronDeltaPt[i][j]->GetYaxis()->CenterTitle(true); 
      hElectronDeltaPt[i][j]->SetLineWidth(2);
      hElectronDeltaPt[i][j]->Fit("gaus");
      if (hElectronDeltaPt[i][j]->GetFunction("gaus"))
      {
        gaus_hElectronDeltaPt[i][j] = hElectronDeltaPt[i][j]->GetFunction("gaus");
        gaus_hElectronDeltaPt[i][j]->SetLineWidth(2);
        gaus_hElectronDeltaPt[i][j]->SetLineColor(kRed);
        hElectronDeltaPt[i][j]->Draw();
      }
      c[icnv]->SaveAs(Form("./plots/hElectronDeltaEta%dPt%d.png",i,j));
      
      icnv++;
      c[icnv] = new TCanvas(Form("c%d",icnv),Form("c%d",icnv));
      hElectronEhadOverEemPt[i][j]->GetXaxis()->CenterTitle(true); 
      hElectronEhadOverEemPt[i][j]->GetYaxis()->CenterTitle(true); 
      hElectronEhadOverEemPt[i][j]->SetLineWidth(2);
      hElectronEhadOverEemPt[i][j]->Draw();
      c[icnv]->SaveAs(Form("./plots/hElectronEhadOverEemEta%dPt%d.png",i,j));
    }
  }
  // Plot EFlowTrack figures
  for (int i = 0; i < nEtaBin; ++i)
  {
    icnv++;
    c[icnv] = new TCanvas(Form("c%d",icnv),Form("c%d",icnv));
    hEFlowTrackDeltaEta[i]->GetXaxis()->CenterTitle(true); 
    hEFlowTrackDeltaEta[i]->GetYaxis()->CenterTitle(true); 
    hEFlowTrackDeltaEta[i]->SetLineWidth(2);
    hEFlowTrackDeltaEta[i]->Fit("gaus");
    gaus_hEFlowTrackDeltaEta[i] = hEFlowTrackDeltaEta[i]->GetFunction("gaus");
    gaus_hEFlowTrackDeltaEta[i]->SetLineWidth(2);
    gaus_hEFlowTrackDeltaEta[i]->SetLineColor(kRed);
    hEFlowTrackDeltaEta[i]->Draw();
    c[icnv]->SaveAs(Form("./plots/hEFlowTrackDeltaEta%d.png",i));

    for (int j =0; j < nPtBin; ++j)
    {
      icnv++;
      c[icnv] = new TCanvas(Form("c%d",icnv),Form("c%d",icnv));
      hEFlowTrackDeltaPt[i][j]->GetXaxis()->CenterTitle(true); 
      hEFlowTrackDeltaPt[i][j]->GetYaxis()->CenterTitle(true); 
      hEFlowTrackDeltaPt[i][j]->SetLineWidth(2);
      hEFlowTrackDeltaPt[i][j]->Fit("gaus");
      if (hEFlowTrackDeltaPt[i][j]->GetFunction("gaus"))
      {
        gaus_hEFlowTrackDeltaPt[i][j] = hEFlowTrackDeltaPt[i][j]->GetFunction("gaus");
        gaus_hEFlowTrackDeltaPt[i][j]->SetLineWidth(2);
        gaus_hEFlowTrackDeltaPt[i][j]->SetLineColor(kRed);
        hEFlowTrackDeltaPt[i][j]->Draw();
      }
      c[icnv]->SaveAs(Form("./plots/hEFlowTrackDeltaEta%dPt%d.png",i,j));
    }
    for (int k =0; k < nPBin; ++k)
    {
      icnv++;
      c[icnv] = new TCanvas(Form("c%d",icnv),Form("c%d",icnv));
      hEFlowTrackDeltaP[i][k]->GetXaxis()->CenterTitle(true); 
      hEFlowTrackDeltaP[i][k]->GetYaxis()->CenterTitle(true); 
      hEFlowTrackDeltaP[i][k]->SetLineWidth(2);
      hEFlowTrackDeltaP[i][k]->Fit("gaus");
      if (hEFlowTrackDeltaP[i][k]->GetFunction("gaus"))
      {
        gaus_hEFlowTrackDeltaP[i][k] = hEFlowTrackDeltaP[i][k]->GetFunction("gaus");
        gaus_hEFlowTrackDeltaP[i][k]->SetLineWidth(2);
        gaus_hEFlowTrackDeltaP[i][k]->SetLineColor(kRed);
        hEFlowTrackDeltaP[i][k]->Draw();
        c[icnv]->SaveAs(Form("./plots/hEFlowTrackDeltaEta%dP%d.png",i,k));
      }
    }
  }
  // Plot Tower Figures
  for (int i = 0; i < nEtaBin; ++i)
  {
    icnv++;
    c[icnv] = new TCanvas(Form("c%d",icnv),Form("c%d",icnv));
    hTowerDeltaEta[i]->GetXaxis()->CenterTitle(true); 
    hTowerDeltaEta[i]->GetYaxis()->CenterTitle(true); 
    hTowerDeltaEta[i]->SetLineWidth(2);
    hTowerDeltaEta[i]->Fit("gaus");
    gaus_hTowerDeltaEta[i] = hTowerDeltaEta[i]->GetFunction("gaus");
    gaus_hTowerDeltaEta[i]->SetLineWidth(2);
    gaus_hTowerDeltaEta[i]->SetLineColor(kRed);
    hTowerDeltaEta[i]->Draw();
    c[icnv]->SaveAs(Form("./plots/hTowerDeltaEta%d.png",i));

    icnv++;
    c[icnv] = new TCanvas(Form("c%d",icnv),Form("c%d",icnv));
    hTowerDeltaEem[i]->GetXaxis()->CenterTitle(true); 
    hTowerDeltaEem[i]->GetYaxis()->CenterTitle(true); 
    hTowerDeltaEem[i]->SetLineWidth(2);
    hTowerDeltaEem[i]->Fit("gaus");
    gaus_hTowerDeltaEem[i] = hTowerDeltaEem[i]->GetFunction("gaus");
    gaus_hTowerDeltaEem[i]->SetLineWidth(2);
    gaus_hTowerDeltaEem[i]->SetLineColor(kRed);
    hTowerDeltaEem[i]->Draw();
    c[icnv]->SaveAs(Form("./plots/hTowerDeltaEem%d.png",i));

    icnv++;
    c[icnv] = new TCanvas(Form("c%d",icnv),Form("c%d",icnv));
    hTowerDeltaEhad[i]->GetXaxis()->CenterTitle(true); 
    hTowerDeltaEhad[i]->GetYaxis()->CenterTitle(true); 
    hTowerDeltaEhad[i]->SetLineWidth(2);
    hTowerDeltaEhad[i]->Fit("gaus");
    gaus_hTowerDeltaEhad[i] = hTowerDeltaEhad[i]->GetFunction("gaus");
    gaus_hTowerDeltaEhad[i]->SetLineWidth(2);
    gaus_hTowerDeltaEhad[i]->SetLineColor(kRed);
    hTowerDeltaEhad[i]->Draw();
    c[icnv]->SaveAs(Form("./plots/hTowerDeltaEhad%d.png",i));

    for (int k =0; k < nEBin; ++k)
    {
      icnv++;
      c[icnv] = new TCanvas(Form("c%d",icnv),Form("c%d",icnv));
      hTowerDeltaE[i][k]->GetXaxis()->CenterTitle(true); 
      hTowerDeltaE[i][k]->GetYaxis()->CenterTitle(true); 
      hTowerDeltaE[i][k]->SetLineWidth(2);
      hTowerDeltaE[i][k]->Fit("gaus");
      if (hTowerDeltaE[i][k]->GetFunction("gaus"))
      {
        gaus_hTowerDeltaE[i][k] = hTowerDeltaE[i][k]->GetFunction("gaus");
        gaus_hTowerDeltaE[i][k]->SetLineWidth(2);
        gaus_hTowerDeltaE[i][k]->SetLineColor(kRed);
        hTowerDeltaE[i][k]->Draw();
        c[icnv]->SaveAs(Form("./plots/hTowerDeltaEta%dE%d.png",i,k));
      }
    }
  }

  TCanvas *cnv1000 = new TCanvas("cnv1000", "cnv1000");
  hTowerN->GetXaxis()->CenterTitle(true);
  hTowerN->GetYaxis()->CenterTitle(true);
  hTowerN->SetLineWidth(2);
  hTowerN->Draw();
  cnv1000->SaveAs("./plots/hTowerN.png");

  TCanvas *cnv1001 = new TCanvas("cnv1001", "cnv1001");
  hTowerGenN->GetXaxis()->CenterTitle(true);
  hTowerGenN->GetYaxis()->CenterTitle(true);
  hTowerGenN->SetLineWidth(2);
  hTowerGenN->Draw();
  cnv1001->SaveAs("./plots/hTowerGenN.png");

  TCanvas *cnv1002 = new TCanvas("cnv1002", "cnv1002");
  hTowerGenPID->GetXaxis()->CenterTitle(true);
  hTowerGenPID->GetYaxis()->CenterTitle(true);
  hTowerGenPID->SetLineWidth(2);
  hTowerGenPID->Draw();
  cnv1002->SaveAs("./plots/hTowerGenPID.png");

  cout << "** Done..." << endl;

  delete treeReader;
  delete chain;
}
