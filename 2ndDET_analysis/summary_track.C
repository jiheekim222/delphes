//////////////////////////////////////////
// 08/11/2023 Jihee Kim (jkim11@bnl.gov)
// Tracking resolutions
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

void summary_track(const char *inputFile)
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
  // Get total number of events
  Long64_t allEntries = treeReader->GetEntries();
  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Track *track;  
  
  // Book histograms
  // Based on eta and momentum
  // -3.5 < eta <= -3.0
  TH1* hTrackPRes_Eta_N35_N30_P_00_05 = new TH1D("hTrackPRes_Eta_N35_N30_P_00_05","0<P[GeV]<=5   & -3.5<#eta<=-3.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N35_N30_P_05_10 = new TH1D("hTrackPRes_Eta_N35_N30_P_05_10","5<P[GeV]<=10  & -3.5<#eta<=-3.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N35_N30_P_10_15 = new TH1D("hTrackPRes_Eta_N35_N30_P_10_15","10<P[GeV]<=15 & -3.5<#eta<=-3.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N35_N30_P_15_20 = new TH1D("hTrackPRes_Eta_N35_N30_P_15_20","15<P[GeV]<=20 & -3.5<#eta<=-3.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  // -3.0 < eta <= -2.5 
  TH1* hTrackPRes_Eta_N30_N25_P_00_05 = new TH1D("hTrackPRes_Eta_N30_N25_P_00_05","0<P[GeV]<=5   & -3.0<#eta<=-2.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N30_N25_P_05_10 = new TH1D("hTrackPRes_Eta_N30_N25_P_05_10","5<P[GeV]<=10  & -3.0<#eta<=-2.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N30_N25_P_10_15 = new TH1D("hTrackPRes_Eta_N30_N25_P_10_15","10<P[GeV]<=15 & -3.0<#eta<=-2.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N30_N25_P_15_20 = new TH1D("hTrackPRes_Eta_N30_N25_P_15_20","15<P[GeV]<=20 & -3.0<#eta<=-2.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1);
  // -2.5 < eta <= -2.0
  TH1* hTrackPRes_Eta_N25_N20_P_00_05 = new TH1D("hTrackPRes_Eta_N25_N20_P_00_05","0<P[GeV]<=5   & -2.5<#eta<=-2.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N25_N20_P_05_10 = new TH1D("hTrackPRes_Eta_N25_N20_P_05_10","5<P[GeV]<=10  & -2.5<#eta<=-2.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N25_N20_P_10_15 = new TH1D("hTrackPRes_Eta_N25_N20_P_10_15","10<P[GeV]<=15 & -2.5<#eta<=-2.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N25_N20_P_15_20 = new TH1D("hTrackPRes_Eta_N25_N20_P_15_20","15<P[GeV]<=20 & -2.5<#eta<=-2.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  // -2.0 < eta <= -1.5
  TH1* hTrackPRes_Eta_N20_N15_P_00_05 = new TH1D("hTrackPRes_Eta_N20_N15_P_00_05","0<P[GeV]<=5   & -2.0<#eta<=-1.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N20_N15_P_05_10 = new TH1D("hTrackPRes_Eta_N20_N15_P_05_10","5<P[GeV]<=10  & -2.0<#eta<=-1.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N20_N15_P_10_15 = new TH1D("hTrackPRes_Eta_N20_N15_P_10_15","10<P[GeV]<=15 & -2.0<#eta<=-1.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N20_N15_P_15_20 = new TH1D("hTrackPRes_Eta_N20_N15_P_15_20","15<P[GeV]<=20 & -2.0<#eta<=-1.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1);
  // -1.5 < eta <= -1.0
  TH1* hTrackPRes_Eta_N15_N10_P_00_05 = new TH1D("hTrackPRes_Eta_N15_N10_P_00_05","0<P[GeV]<=5   & -1.5<#eta<=-1.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N15_N10_P_05_10 = new TH1D("hTrackPRes_Eta_N15_N10_P_05_10","5<P[GeV]<=10  & -1.5<#eta<=-1.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N15_N10_P_10_15 = new TH1D("hTrackPRes_Eta_N15_N10_P_10_15","10<P[GeV]<=15 & -1.5<#eta<=-1.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N15_N10_P_15_20 = new TH1D("hTrackPRes_Eta_N15_N10_P_15_20","15<P[GeV]<=20 & -1.5<#eta<=-1.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1);
  // -1.0 < eta <= -0.5
  TH1* hTrackPRes_Eta_N10_N05_P_00_05 = new TH1D("hTrackPRes_Eta_N10_N05_P_00_05","0<P[GeV]<=5   & -1.0<#eta<=-0.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N10_N05_P_05_10 = new TH1D("hTrackPRes_Eta_N10_N05_P_05_10","5<P[GeV]<=10  & -1.0<#eta<=-0.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N10_N05_P_10_15 = new TH1D("hTrackPRes_Eta_N10_N05_P_10_15","10<P[GeV]<=15 & -1.0<#eta<=-0.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N10_N05_P_15_20 = new TH1D("hTrackPRes_Eta_N10_N05_P_15_20","15<P[GeV]<=20 & -1.0<#eta<=-0.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  // -0.5 < eta <= 0.0
  TH1* hTrackPRes_Eta_N05_000_P_00_05 = new TH1D("hTrackPRes_Eta_N05_000_P_00_05","0<P[GeV]<=5   & -0.5<#eta<=-0.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N05_000_P_05_10 = new TH1D("hTrackPRes_Eta_N05_000_P_05_10","5<P[GeV]<=10  & -0.5<#eta<=-0.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N05_000_P_10_15 = new TH1D("hTrackPRes_Eta_N05_000_P_10_15","10<P[GeV]<=15 & -0.5<#eta<=-0.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_N05_000_P_15_20 = new TH1D("hTrackPRes_Eta_N05_000_P_15_20","15<P[GeV]<=20 & -0.5<#eta<=-0.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  // 0.0 < eta <= 0.5
  TH1* hTrackPRes_Eta_000_P05_P_00_05 = new TH1D("hTrackPRes_Eta_000_P05_P_00_05","0<P[GeV]<=5   & 0.0<#eta<=0.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_000_P05_P_05_10 = new TH1D("hTrackPRes_Eta_000_P05_P_05_10","5<P[GeV]<=10  & 0.0<#eta<=0.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_000_P05_P_10_15 = new TH1D("hTrackPRes_Eta_000_P05_P_10_15","10<P[GeV]<=15 & 0.0<#eta<=0.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_000_P05_P_15_20 = new TH1D("hTrackPRes_Eta_000_P05_P_15_20","15<P[GeV]<=20 & 0.0<#eta<=0.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1);
  // 0.5 < eta <= 1.0
  TH1* hTrackPRes_Eta_P05_P10_P_00_05 = new TH1D("hTrackPRes_Eta_P05_P10_P_00_05","0<P[GeV]<=5   & 0.5<#eta<=1.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P05_P10_P_05_10 = new TH1D("hTrackPRes_Eta_P05_P10_P_05_10","5<P[GeV]<=10  & 0.5<#eta<=1.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P05_P10_P_10_15 = new TH1D("hTrackPRes_Eta_P05_P10_P_10_15","10<P[GeV]<=15 & 0.5<#eta<=1.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P05_P10_P_15_20 = new TH1D("hTrackPRes_Eta_P05_P10_P_15_20","15<P[GeV]<=20 & 0.5<#eta<=1.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1);
  // 1.0 < eta <= 1.5
  TH1* hTrackPRes_Eta_P10_P15_P_00_05 = new TH1D("hTrackPRes_Eta_P10_P15_P_00_05","0<P[GeV]<=5   & 1.0<#eta<=1.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P10_P15_P_05_10 = new TH1D("hTrackPRes_Eta_P10_P15_P_05_10","5<P[GeV]<=10  & 1.0<#eta<=1.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P10_P15_P_10_15 = new TH1D("hTrackPRes_Eta_P10_P15_P_10_15","10<P[GeV]<=15 & 1.0<#eta<=1.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P10_P15_P_15_20 = new TH1D("hTrackPRes_Eta_P10_P15_P_15_20","15<P[GeV]<=20 & 1.0<#eta<=1.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1);
  // 1.5 < eta <= 2.0
  TH1* hTrackPRes_Eta_P15_P20_P_00_05 = new TH1D("hTrackPRes_Eta_P15_P20_P_00_05","0<P[GeV]<=5   & 1.5<#eta<=2.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P15_P20_P_05_10 = new TH1D("hTrackPRes_Eta_P15_P20_P_05_10","5<P[GeV]<=10  & 1.5<#eta<=2.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P15_P20_P_10_15 = new TH1D("hTrackPRes_Eta_P15_P20_P_10_15","10<P[GeV]<=15 & 1.5<#eta<=2.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P15_P20_P_15_20 = new TH1D("hTrackPRes_Eta_P15_P20_P_15_20","15<P[GeV]<=20 & 1.5<#eta<=2.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1);
  // 2.0 < eta <= 2.5
  TH1* hTrackPRes_Eta_P20_P25_P_00_05 = new TH1D("hTrackPRes_Eta_P20_P25_P_00_05","0<P[GeV]<=5   & 2.0<#eta<=2.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P20_P25_P_05_10 = new TH1D("hTrackPRes_Eta_P20_P25_P_05_10","5<P[GeV]<=10  & 2.0<#eta<=2.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P20_P25_P_10_15 = new TH1D("hTrackPRes_Eta_P20_P25_P_10_15","10<P[GeV]<=15 & 2.0<#eta<=2.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P20_P25_P_15_20 = new TH1D("hTrackPRes_Eta_P20_P25_P_15_20","15<P[GeV]<=20 & 2.0<#eta<=2.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1);
  // 2.5 < eta <= 3.0
  TH1* hTrackPRes_Eta_P25_P30_P_00_05 = new TH1D("hTrackPRes_Eta_P25_P30_P_00_05","0<P[GeV]<=5   & 2.5<#eta<=3.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P25_P30_P_05_10 = new TH1D("hTrackPRes_Eta_P25_P30_P_05_10","5<P[GeV]<=10  & 2.5<#eta<=3.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P25_P30_P_10_15 = new TH1D("hTrackPRes_Eta_P25_P30_P_10_15","10<P[GeV]<=15 & 2.5<#eta<=3.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P25_P30_P_15_20 = new TH1D("hTrackPRes_Eta_P25_P30_P_15_20","15<P[GeV]<=20 & 2.5<#eta<=3.0;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1);
  // 3.0 < eta <= 3.5
  TH1* hTrackPRes_Eta_P30_P35_P_00_05 = new TH1D("hTrackPRes_Eta_P30_P35_P_00_05","0<P[GeV]<=5   & 3.0<#eta<=3.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P30_P35_P_05_10 = new TH1D("hTrackPRes_Eta_P30_P35_P_05_10","5<P[GeV]<=10  & 3.0<#eta<=3.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P30_P35_P_10_15 = new TH1D("hTrackPRes_Eta_P30_P35_P_10_15","10<P[GeV]<=15 & 3.0<#eta<=3.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1); 
  TH1* hTrackPRes_Eta_P30_P35_P_15_20 = new TH1D("hTrackPRes_Eta_P30_P35_P_15_20","15<P[GeV]<=20 & 3.0<#eta<=3.5;(P^{rec}-P^{gen})/P^{gen};Number of tracks",100, -0.1, 0.1);

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
	
	// Track momentum resolution
        if (particle->P > 0.0 && particle->P <= 5.0)
	{
          if (particle->Eta > -3.5 && particle->Eta <= -3.0)
	    hTrackPRes_Eta_N35_N30_P_00_05->Fill((track->P - particle->P)/particle->P);
	  else if (particle->Eta > -3.0 && particle->Eta <= -2.5)
	    hTrackPRes_Eta_N30_N25_P_00_05->Fill((track->P - particle->P)/particle->P);
	  else if (particle->Eta > -2.5 && particle->Eta <= -2.0)
	    hTrackPRes_Eta_N25_N20_P_00_05->Fill((track->P - particle->P)/particle->P);
	  else if (particle->Eta > -2.0 && particle->Eta <= -1.5)
	    hTrackPRes_Eta_N20_N15_P_00_05->Fill((track->P - particle->P)/particle->P);
	  else if (particle->Eta > -1.5 && particle->Eta <= -1.0)
	    hTrackPRes_Eta_N15_N10_P_00_05->Fill((track->P - particle->P)/particle->P);
	  else if (particle->Eta > -1.0 && particle->Eta <= -0.5)
	    hTrackPRes_Eta_N10_N05_P_00_05->Fill((track->P - particle->P)/particle->P);
	  else if (particle->Eta > -0.5 && particle->Eta <=  0.0)
	    hTrackPRes_Eta_N05_000_P_00_05->Fill((track->P - particle->P)/particle->P);
	  else if (particle->Eta >  0.0 && particle->Eta <=  0.5)
	    hTrackPRes_Eta_000_P05_P_00_05->Fill((track->P - particle->P)/particle->P);
	  else if (particle->Eta >  0.5 && particle->Eta <=  1.0)
	    hTrackPRes_Eta_P05_P10_P_00_05->Fill((track->P - particle->P)/particle->P);
	  else if (particle->Eta >  1.0 && particle->Eta <=  1.5)
	    hTrackPRes_Eta_P10_P15_P_00_05->Fill((track->P - particle->P)/particle->P);
	  else if (particle->Eta >  1.5 && particle->Eta <=  2.0)
	    hTrackPRes_Eta_P15_P20_P_00_05->Fill((track->P - particle->P)/particle->P);
	  else if (particle->Eta >  2.0 && particle->Eta <=  2.5)
	    hTrackPRes_Eta_P20_P25_P_00_05->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  2.5 && particle->Eta <=  3.0)
	    hTrackPRes_Eta_P25_P30_P_00_05->Fill((track->P - particle->P)/particle->P);
	  else if (particle->Eta >  3.0 && particle->Eta <=  3.5)
	    hTrackPRes_Eta_P30_P35_P_00_05->Fill((track->P - particle->P)/particle->P);
	  else
	    continue;
	}
	else if (particle->P > 5.0 && particle->P <= 10.0)
	{
          if (particle->Eta > -3.5 && particle->Eta <= -3.0)
            hTrackPRes_Eta_N35_N30_P_05_10->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -3.0 && particle->Eta <= -2.5)
            hTrackPRes_Eta_N30_N25_P_05_10->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -2.5 && particle->Eta <= -2.0)
            hTrackPRes_Eta_N25_N20_P_05_10->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -2.0 && particle->Eta <= -1.5)
            hTrackPRes_Eta_N20_N15_P_05_10->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -1.5 && particle->Eta <= -1.0)
            hTrackPRes_Eta_N15_N10_P_05_10->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -1.0 && particle->Eta <= -0.5)
            hTrackPRes_Eta_N10_N05_P_05_10->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -0.5 && particle->Eta <=  0.0)
            hTrackPRes_Eta_N05_000_P_05_10->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  0.0 && particle->Eta <=  0.5)
            hTrackPRes_Eta_000_P05_P_05_10->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  0.5 && particle->Eta <=  1.0)
            hTrackPRes_Eta_P05_P10_P_05_10->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  1.0 && particle->Eta <=  1.5)
            hTrackPRes_Eta_P10_P15_P_05_10->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  1.5 && particle->Eta <=  2.0)
            hTrackPRes_Eta_P15_P20_P_05_10->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  2.0 && particle->Eta <=  2.5)
            hTrackPRes_Eta_P20_P25_P_05_10->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  2.5 && particle->Eta <=  3.0)
            hTrackPRes_Eta_P25_P30_P_05_10->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  3.0 && particle->Eta <=  3.5)
            hTrackPRes_Eta_P30_P35_P_05_10->Fill((track->P - particle->P)/particle->P);
          else
            continue;
	}
	else if (particle->P > 10.0 && particle->P <= 15.0)
	{
          if (particle->Eta > -3.5 && particle->Eta <= -3.0)
            hTrackPRes_Eta_N35_N30_P_10_15->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -3.0 && particle->Eta <= -2.5)
            hTrackPRes_Eta_N30_N25_P_10_15->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -2.5 && particle->Eta <= -2.0)
            hTrackPRes_Eta_N25_N20_P_10_15->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -2.0 && particle->Eta <= -1.5)
            hTrackPRes_Eta_N20_N15_P_10_15->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -1.5 && particle->Eta <= -1.0)
            hTrackPRes_Eta_N15_N10_P_10_15->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -1.0 && particle->Eta <= -0.5)
            hTrackPRes_Eta_N10_N05_P_10_15->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -0.5 && particle->Eta <=  0.0)
            hTrackPRes_Eta_N05_000_P_10_15->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  0.0 && particle->Eta <=  0.5)
            hTrackPRes_Eta_000_P05_P_10_15->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  0.5 && particle->Eta <=  1.0)
            hTrackPRes_Eta_P05_P10_P_10_15->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  1.0 && particle->Eta <=  1.5)
            hTrackPRes_Eta_P10_P15_P_10_15->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  1.5 && particle->Eta <=  2.0)
            hTrackPRes_Eta_P15_P20_P_10_15->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  2.0 && particle->Eta <=  2.5)
            hTrackPRes_Eta_P20_P25_P_10_15->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  2.5 && particle->Eta <=  3.0)
            hTrackPRes_Eta_P25_P30_P_10_15->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  3.0 && particle->Eta <=  3.5)
            hTrackPRes_Eta_P30_P35_P_10_15->Fill((track->P - particle->P)/particle->P);
          else
            continue;
	}
	else if (particle->P > 15.0 && particle->P <= 20.0)
	{
          if (particle->Eta > -3.5 && particle->Eta <= -3.0)
            hTrackPRes_Eta_N35_N30_P_15_20->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -3.0 && particle->Eta <= -2.5)
            hTrackPRes_Eta_N30_N25_P_15_20->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -2.5 && particle->Eta <= -2.0)
            hTrackPRes_Eta_N25_N20_P_15_20->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -2.0 && particle->Eta <= -1.5)
            hTrackPRes_Eta_N20_N15_P_15_20->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -1.5 && particle->Eta <= -1.0)
            hTrackPRes_Eta_N15_N10_P_15_20->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -1.0 && particle->Eta <= -0.5)
            hTrackPRes_Eta_N10_N05_P_15_20->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta > -0.5 && particle->Eta <=  0.0)
            hTrackPRes_Eta_N05_000_P_15_20->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  0.0 && particle->Eta <=  0.5)
            hTrackPRes_Eta_000_P05_P_15_20->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  0.5 && particle->Eta <=  1.0)
            hTrackPRes_Eta_P05_P10_P_15_20->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  1.0 && particle->Eta <=  1.5)
            hTrackPRes_Eta_P10_P15_P_15_20->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  1.5 && particle->Eta <=  2.0)
            hTrackPRes_Eta_P15_P20_P_15_20->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  2.0 && particle->Eta <=  2.5)
            hTrackPRes_Eta_P20_P25_P_15_20->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  2.5 && particle->Eta <=  3.0)
            hTrackPRes_Eta_P25_P30_P_15_20->Fill((track->P - particle->P)/particle->P);
          else if (particle->Eta >  3.0 && particle->Eta <=  3.5)
            hTrackPRes_Eta_P30_P35_P_15_20->Fill((track->P - particle->P)/particle->P);
          else
            continue;
	}
	else
	  continue;
      }
    }
  }

  // Plot figures
  TCanvas *cnv1 = new TCanvas("cnv1", "cnv1", 3500, 1000);
  cnv1->Divide(7,2);
  cnv1->cd(1);
  hTrackPRes_Eta_N35_N30_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N35_N30_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N35_N30_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_N35_N30_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N35_N30_P_00_05 = hTrackPRes_Eta_N35_N30_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N35_N30_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N35_N30_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_N35_N30_P_00_05->Draw();
  cnv1->cd(2);
  hTrackPRes_Eta_N30_N25_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N30_N25_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N30_N25_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_N30_N25_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N30_N25_P_00_05 = hTrackPRes_Eta_N30_N25_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N30_N25_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N30_N25_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_N30_N25_P_00_05->Draw();
  cnv1->cd(3);
  hTrackPRes_Eta_N25_N20_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N25_N20_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N25_N20_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_N25_N20_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N25_N20_P_00_05 =hTrackPRes_Eta_N25_N20_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N25_N20_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N25_N20_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_N25_N20_P_00_05->Draw();
  cnv1->cd(4);
  hTrackPRes_Eta_N20_N15_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N20_N15_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N20_N15_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_N20_N15_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N20_N15_P_00_05 = hTrackPRes_Eta_N20_N15_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N20_N15_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N20_N15_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_N20_N15_P_00_05->Draw();
  cnv1->cd(5);
  hTrackPRes_Eta_N15_N10_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N15_N10_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N15_N10_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_N15_N10_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N15_N10_P_00_05 = hTrackPRes_Eta_N15_N10_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N15_N10_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N15_N10_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_N15_N10_P_00_05->Draw();
  cnv1->cd(6);
  hTrackPRes_Eta_N10_N05_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N10_N05_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N10_N05_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_N10_N05_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N10_N05_P_00_05 = hTrackPRes_Eta_N10_N05_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N10_N05_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N10_N05_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_N10_N05_P_00_05->Draw();
  cnv1->cd(7);
  hTrackPRes_Eta_N05_000_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N05_000_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N05_000_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_N05_000_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N05_000_P_00_05 = hTrackPRes_Eta_N05_000_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N05_000_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N05_000_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_N05_000_P_00_05->Draw();
  cnv1->cd(8);
  hTrackPRes_Eta_000_P05_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_000_P05_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_000_P05_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_000_P05_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_000_P05_P_00_05 = hTrackPRes_Eta_000_P05_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_000_P05_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_000_P05_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_000_P05_P_00_05->Draw();
  cnv1->cd(9);
  hTrackPRes_Eta_P05_P10_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P05_P10_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P05_P10_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_P05_P10_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P05_P10_P_00_05 = hTrackPRes_Eta_P05_P10_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P05_P10_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P05_P10_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_P05_P10_P_00_05->Draw();
  cnv1->cd(10);
  hTrackPRes_Eta_P10_P15_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P10_P15_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P10_P15_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_P10_P15_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P10_P15_P_00_05 = hTrackPRes_Eta_P10_P15_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P10_P15_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P10_P15_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_P10_P15_P_00_05->Draw();
  cnv1->cd(11);
  hTrackPRes_Eta_P15_P20_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P15_P20_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P15_P20_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_P15_P20_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P15_P20_P_00_05 = hTrackPRes_Eta_P15_P20_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P15_P20_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P15_P20_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_P15_P20_P_00_05->Draw();
  cnv1->cd(12);
  hTrackPRes_Eta_P20_P25_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P20_P25_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P20_P25_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_P20_P25_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P20_P25_P_00_05 = hTrackPRes_Eta_P20_P25_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P20_P25_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P20_P25_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_P20_P25_P_00_05->Draw();
  cnv1->cd(13);
  hTrackPRes_Eta_P25_P30_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P25_P30_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P25_P30_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_P25_P30_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P25_P30_P_00_05 = hTrackPRes_Eta_P25_P30_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P25_P30_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P25_P30_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_P25_P30_P_00_05->Draw();
  cnv1->cd(14);
  hTrackPRes_Eta_P30_P35_P_00_05->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P30_P35_P_00_05->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P30_P35_P_00_05->SetLineWidth(2);
  hTrackPRes_Eta_P30_P35_P_00_05->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P30_P35_P_00_05 = hTrackPRes_Eta_P30_P35_P_00_05->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P30_P35_P_00_05->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P30_P35_P_00_05->SetLineColor(kRed);
  hTrackPRes_Eta_P30_P35_P_00_05->Draw();
  cnv1->SaveAs("./plots/hTrackPRes_P_00_05.png");
  
  TCanvas *cnv2 = new TCanvas("cnv2", "cnv2", 3500, 1000);
  cnv2->Divide(7,2);
  cnv2->cd(1);
  hTrackPRes_Eta_N35_N30_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N35_N30_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N35_N30_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_N35_N30_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N35_N30_P_05_10 = hTrackPRes_Eta_N35_N30_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N35_N30_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N35_N30_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_N35_N30_P_05_10->Draw();
  cnv2->cd(2);
  hTrackPRes_Eta_N30_N25_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N30_N25_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N30_N25_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_N30_N25_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N30_N25_P_05_10 = hTrackPRes_Eta_N30_N25_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N30_N25_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N30_N25_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_N30_N25_P_05_10->Draw();
  cnv2->cd(3);
  hTrackPRes_Eta_N25_N20_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N25_N20_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N25_N20_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_N25_N20_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N25_N20_P_05_10 =hTrackPRes_Eta_N25_N20_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N25_N20_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N25_N20_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_N25_N20_P_05_10->Draw();
  cnv2->cd(4);
  hTrackPRes_Eta_N20_N15_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N20_N15_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N20_N15_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_N20_N15_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N20_N15_P_05_10 = hTrackPRes_Eta_N20_N15_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N20_N15_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N20_N15_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_N20_N15_P_05_10->Draw();
  cnv2->cd(5);
  hTrackPRes_Eta_N15_N10_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N15_N10_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N15_N10_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_N15_N10_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N15_N10_P_05_10 = hTrackPRes_Eta_N15_N10_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N15_N10_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N15_N10_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_N15_N10_P_05_10->Draw();
  cnv2->cd(6);
  hTrackPRes_Eta_N10_N05_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N10_N05_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N10_N05_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_N10_N05_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N10_N05_P_05_10 = hTrackPRes_Eta_N10_N05_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N10_N05_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N10_N05_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_N10_N05_P_05_10->Draw();
  cnv2->cd(7);
  hTrackPRes_Eta_N05_000_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N05_000_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N05_000_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_N05_000_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N05_000_P_05_10 = hTrackPRes_Eta_N05_000_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N05_000_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N05_000_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_N05_000_P_05_10->Draw();
  cnv2->cd(8);
  hTrackPRes_Eta_000_P05_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_000_P05_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_000_P05_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_000_P05_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_000_P05_P_05_10 = hTrackPRes_Eta_000_P05_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_000_P05_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_000_P05_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_000_P05_P_05_10->Draw();
  cnv2->cd(9);
  hTrackPRes_Eta_P05_P10_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P05_P10_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P05_P10_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_P05_P10_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P05_P10_P_05_10 = hTrackPRes_Eta_P05_P10_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P05_P10_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P05_P10_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_P05_P10_P_05_10->Draw();
  cnv2->cd(10);
  hTrackPRes_Eta_P10_P15_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P10_P15_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P10_P15_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_P10_P15_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P10_P15_P_05_10 = hTrackPRes_Eta_P10_P15_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P10_P15_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P10_P15_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_P10_P15_P_05_10->Draw();
  cnv2->cd(11);
  hTrackPRes_Eta_P15_P20_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P15_P20_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P15_P20_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_P15_P20_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P15_P20_P_05_10 = hTrackPRes_Eta_P15_P20_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P15_P20_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P15_P20_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_P15_P20_P_05_10->Draw();
  cnv2->cd(12);
  hTrackPRes_Eta_P20_P25_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P20_P25_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P20_P25_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_P20_P25_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P20_P25_P_05_10 = hTrackPRes_Eta_P20_P25_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P20_P25_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P20_P25_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_P20_P25_P_05_10->Draw();
  cnv2->cd(13);
  hTrackPRes_Eta_P25_P30_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P25_P30_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P25_P30_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_P25_P30_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P25_P30_P_05_10 = hTrackPRes_Eta_P25_P30_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P25_P30_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P25_P30_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_P25_P30_P_05_10->Draw();
  cnv2->cd(14);
  hTrackPRes_Eta_P30_P35_P_05_10->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P30_P35_P_05_10->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P30_P35_P_05_10->SetLineWidth(2);
  hTrackPRes_Eta_P30_P35_P_05_10->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P30_P35_P_05_10 = hTrackPRes_Eta_P30_P35_P_05_10->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P30_P35_P_05_10->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P30_P35_P_05_10->SetLineColor(kRed);
  hTrackPRes_Eta_P30_P35_P_05_10->Draw();
  cnv2->SaveAs("./plots/hTrackPRes_P_05_10.png");

  TCanvas *cnv3 = new TCanvas("cnv3", "cnv3", 3500, 1000);
  cnv3->Divide(7,2);
  cnv3->cd(1);
  hTrackPRes_Eta_N35_N30_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N35_N30_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N35_N30_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_N35_N30_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N35_N30_P_10_15 = hTrackPRes_Eta_N35_N30_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N35_N30_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N35_N30_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_N35_N30_P_10_15->Draw();
  cnv3->cd(2);
  hTrackPRes_Eta_N30_N25_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N30_N25_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N30_N25_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_N30_N25_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N30_N25_P_10_15 = hTrackPRes_Eta_N30_N25_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N30_N25_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N30_N25_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_N30_N25_P_10_15->Draw();
  cnv3->cd(3);
  hTrackPRes_Eta_N25_N20_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N25_N20_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N25_N20_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_N25_N20_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N25_N20_P_10_15 =hTrackPRes_Eta_N25_N20_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N25_N20_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N25_N20_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_N25_N20_P_10_15->Draw();
  cnv3->cd(4);
  hTrackPRes_Eta_N20_N15_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N20_N15_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N20_N15_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_N20_N15_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N20_N15_P_10_15 = hTrackPRes_Eta_N20_N15_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N20_N15_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N20_N15_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_N20_N15_P_10_15->Draw();
  cnv3->cd(5);
  hTrackPRes_Eta_N15_N10_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N15_N10_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N15_N10_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_N15_N10_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N15_N10_P_10_15 = hTrackPRes_Eta_N15_N10_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N15_N10_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N15_N10_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_N15_N10_P_10_15->Draw();
  cnv3->cd(6);
  hTrackPRes_Eta_N10_N05_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N10_N05_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N10_N05_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_N10_N05_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N10_N05_P_10_15 = hTrackPRes_Eta_N10_N05_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N10_N05_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N10_N05_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_N10_N05_P_10_15->Draw();
  cnv3->cd(7);
  hTrackPRes_Eta_N05_000_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N05_000_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N05_000_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_N05_000_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N05_000_P_10_15 = hTrackPRes_Eta_N05_000_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N05_000_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N05_000_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_N05_000_P_10_15->Draw();
  cnv3->cd(8);
  hTrackPRes_Eta_000_P05_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_000_P05_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_000_P05_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_000_P05_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_000_P05_P_10_15 = hTrackPRes_Eta_000_P05_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_000_P05_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_000_P05_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_000_P05_P_10_15->Draw();
  cnv3->cd(9);
  hTrackPRes_Eta_P05_P10_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P05_P10_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P05_P10_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_P05_P10_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P05_P10_P_10_15 = hTrackPRes_Eta_P05_P10_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P05_P10_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P05_P10_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_P05_P10_P_10_15->Draw();
  cnv3->cd(10);
  hTrackPRes_Eta_P10_P15_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P10_P15_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P10_P15_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_P10_P15_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P10_P15_P_10_15 = hTrackPRes_Eta_P10_P15_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P10_P15_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P10_P15_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_P10_P15_P_10_15->Draw();
  cnv3->cd(11);
  hTrackPRes_Eta_P15_P20_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P15_P20_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P15_P20_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_P15_P20_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P15_P20_P_10_15 = hTrackPRes_Eta_P15_P20_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P15_P20_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P15_P20_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_P15_P20_P_10_15->Draw();
  cnv3->cd(12);
  hTrackPRes_Eta_P20_P25_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P20_P25_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P20_P25_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_P20_P25_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P20_P25_P_10_15 = hTrackPRes_Eta_P20_P25_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P20_P25_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P20_P25_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_P20_P25_P_10_15->Draw();
  cnv3->cd(13);
  hTrackPRes_Eta_P25_P30_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P25_P30_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P25_P30_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_P25_P30_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P25_P30_P_10_15 = hTrackPRes_Eta_P25_P30_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P25_P30_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P25_P30_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_P25_P30_P_10_15->Draw();
  cnv3->cd(14);
  hTrackPRes_Eta_P30_P35_P_10_15->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P30_P35_P_10_15->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P30_P35_P_10_15->SetLineWidth(2);
  hTrackPRes_Eta_P30_P35_P_10_15->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P30_P35_P_10_15 = hTrackPRes_Eta_P30_P35_P_10_15->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P30_P35_P_10_15->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P30_P35_P_10_15->SetLineColor(kRed);
  hTrackPRes_Eta_P30_P35_P_10_15->Draw();
  cnv3->SaveAs("./plots/hTrackPRes_P_10_15.png");
  
  TCanvas *cnv4 = new TCanvas("cnv4", "cnv4", 3500, 1000);
  cnv4->Divide(7,2);
  cnv4->cd(1);
  hTrackPRes_Eta_N35_N30_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N35_N30_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N35_N30_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_N35_N30_P_15_20->Fit("gaus");
  if (hTrackPRes_Eta_N35_N30_P_15_20->GetFunction("gaus"))
  {
    TF1 *gaus_hTrackPRes_Eta_N35_N30_P_15_20 = hTrackPRes_Eta_N35_N30_P_15_20->GetFunction("gaus");
    gaus_hTrackPRes_Eta_N35_N30_P_15_20->SetLineWidth(2);
    gaus_hTrackPRes_Eta_N35_N30_P_15_20->SetLineColor(kRed);
    hTrackPRes_Eta_N35_N30_P_15_20->Draw();
  }
  cnv4->cd(2);
  hTrackPRes_Eta_N30_N25_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N30_N25_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N30_N25_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_N30_N25_P_15_20->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N30_N25_P_15_20 = hTrackPRes_Eta_N30_N25_P_15_20->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N30_N25_P_15_20->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N30_N25_P_15_20->SetLineColor(kRed);
  hTrackPRes_Eta_N30_N25_P_15_20->Draw();
  cnv4->cd(3);
  hTrackPRes_Eta_N25_N20_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N25_N20_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N25_N20_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_N25_N20_P_15_20->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N25_N20_P_15_20 =hTrackPRes_Eta_N25_N20_P_15_20->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N25_N20_P_15_20->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N25_N20_P_15_20->SetLineColor(kRed);
  hTrackPRes_Eta_N25_N20_P_15_20->Draw();
  cnv4->cd(4);
  hTrackPRes_Eta_N20_N15_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N20_N15_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N20_N15_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_N20_N15_P_15_20->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N20_N15_P_15_20 = hTrackPRes_Eta_N20_N15_P_15_20->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N20_N15_P_15_20->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N20_N15_P_15_20->SetLineColor(kRed);
  hTrackPRes_Eta_N20_N15_P_15_20->Draw();
  cnv4->cd(5);
  hTrackPRes_Eta_N15_N10_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N15_N10_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N15_N10_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_N15_N10_P_15_20->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N15_N10_P_15_20 = hTrackPRes_Eta_N15_N10_P_15_20->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N15_N10_P_15_20->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N15_N10_P_15_20->SetLineColor(kRed);
  hTrackPRes_Eta_N15_N10_P_15_20->Draw();
  cnv4->cd(6);
  hTrackPRes_Eta_N10_N05_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N10_N05_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N10_N05_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_N10_N05_P_15_20->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N10_N05_P_15_20 = hTrackPRes_Eta_N10_N05_P_15_20->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N10_N05_P_15_20->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N10_N05_P_15_20->SetLineColor(kRed);
  hTrackPRes_Eta_N10_N05_P_15_20->Draw();
  cnv4->cd(7);
  hTrackPRes_Eta_N05_000_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_N05_000_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_N05_000_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_N05_000_P_15_20->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_N05_000_P_15_20 = hTrackPRes_Eta_N05_000_P_15_20->GetFunction("gaus");
  gaus_hTrackPRes_Eta_N05_000_P_15_20->SetLineWidth(2);
  gaus_hTrackPRes_Eta_N05_000_P_15_20->SetLineColor(kRed);
  hTrackPRes_Eta_N05_000_P_15_20->Draw();
  cnv4->cd(8);
  hTrackPRes_Eta_000_P05_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_000_P05_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_000_P05_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_000_P05_P_15_20->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_000_P05_P_15_20 = hTrackPRes_Eta_000_P05_P_15_20->GetFunction("gaus");
  gaus_hTrackPRes_Eta_000_P05_P_15_20->SetLineWidth(2);
  gaus_hTrackPRes_Eta_000_P05_P_15_20->SetLineColor(kRed);
  hTrackPRes_Eta_000_P05_P_15_20->Draw();
  cnv4->cd(9);
  hTrackPRes_Eta_P05_P10_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P05_P10_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P05_P10_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_P05_P10_P_15_20->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P05_P10_P_15_20 = hTrackPRes_Eta_P05_P10_P_15_20->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P05_P10_P_15_20->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P05_P10_P_15_20->SetLineColor(kRed);
  hTrackPRes_Eta_P05_P10_P_15_20->Draw();
  cnv4->cd(10);
  hTrackPRes_Eta_P10_P15_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P10_P15_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P10_P15_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_P10_P15_P_15_20->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P10_P15_P_15_20 = hTrackPRes_Eta_P10_P15_P_15_20->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P10_P15_P_15_20->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P10_P15_P_15_20->SetLineColor(kRed);
  hTrackPRes_Eta_P10_P15_P_15_20->Draw();
  cnv4->cd(11);
  hTrackPRes_Eta_P15_P20_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P15_P20_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P15_P20_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_P15_P20_P_15_20->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P15_P20_P_15_20 = hTrackPRes_Eta_P15_P20_P_15_20->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P15_P20_P_15_20->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P15_P20_P_15_20->SetLineColor(kRed);
  hTrackPRes_Eta_P15_P20_P_15_20->Draw();
  cnv4->cd(12);
  hTrackPRes_Eta_P20_P25_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P20_P25_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P20_P25_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_P20_P25_P_15_20->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P20_P25_P_15_20 = hTrackPRes_Eta_P20_P25_P_15_20->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P20_P25_P_15_20->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P20_P25_P_15_20->SetLineColor(kRed);
  hTrackPRes_Eta_P20_P25_P_15_20->Draw();
  cnv4->cd(13);
  hTrackPRes_Eta_P25_P30_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P25_P30_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P25_P30_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_P25_P30_P_15_20->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P25_P30_P_15_20 = hTrackPRes_Eta_P25_P30_P_15_20->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P25_P30_P_15_20->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P25_P30_P_15_20->SetLineColor(kRed);
  hTrackPRes_Eta_P25_P30_P_15_20->Draw();
  cnv4->cd(14);
  hTrackPRes_Eta_P30_P35_P_15_20->GetXaxis()->CenterTitle(true);
  hTrackPRes_Eta_P30_P35_P_15_20->GetYaxis()->CenterTitle(true);
  hTrackPRes_Eta_P30_P35_P_15_20->SetLineWidth(2);
  hTrackPRes_Eta_P30_P35_P_15_20->Fit("gaus");
  TF1 *gaus_hTrackPRes_Eta_P30_P35_P_15_20 = hTrackPRes_Eta_P30_P35_P_15_20->GetFunction("gaus");
  gaus_hTrackPRes_Eta_P30_P35_P_15_20->SetLineWidth(2);
  gaus_hTrackPRes_Eta_P30_P35_P_15_20->SetLineColor(kRed);
  hTrackPRes_Eta_P30_P35_P_15_20->Draw();
  cnv4->SaveAs("./plots/hTrackPRes_P_15_20.png");
  
  const int n = 4;
  Double_t x[n] = {2.5, 7.5, 12.5, 17.5};
  Double_t x_error[n] = {0., 0., 0., 0.};
  
  Double_t y_Eta_N35_N30[n] = {gaus_hTrackPRes_Eta_N35_N30_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_N35_N30_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_N35_N30_P_10_15->GetParameter(2)*100., 0.0};
  Double_t y_Eta_N35_N30_error[n] = {gaus_hTrackPRes_Eta_N35_N30_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_N35_N30_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_N35_N30_P_10_15->GetParError(2)*100., 0.0};
  Double_t y_Eta_N30_N25[n] = {gaus_hTrackPRes_Eta_N30_N25_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_N30_N25_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_N30_N25_P_10_15->GetParameter(2)*100., gaus_hTrackPRes_Eta_N30_N25_P_05_10->GetParameter(2)*100.};
  Double_t y_Eta_N30_N25_error[n] = {gaus_hTrackPRes_Eta_N30_N25_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_N30_N25_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_N30_N25_P_10_15->GetParError(2)*100., gaus_hTrackPRes_Eta_N30_N25_P_05_10->GetParError(2)*100.};
  Double_t y_Eta_N25_N20[n] = {gaus_hTrackPRes_Eta_N25_N20_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_N25_N20_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_N25_N20_P_10_15->GetParameter(2)*100., gaus_hTrackPRes_Eta_N25_N20_P_15_20->GetParameter(2)*100.};
  Double_t y_Eta_N25_N20_error[n] = {gaus_hTrackPRes_Eta_N25_N20_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_N25_N20_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_N25_N20_P_10_15->GetParError(2)*100., gaus_hTrackPRes_Eta_N25_N20_P_15_20->GetParError(2)*100.};
  Double_t y_Eta_N20_N15[n] = {gaus_hTrackPRes_Eta_N20_N15_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_N20_N15_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_N20_N15_P_10_15->GetParameter(2)*100., gaus_hTrackPRes_Eta_N20_N15_P_15_20->GetParameter(2)*100.};
  Double_t y_Eta_N20_N15_error[n] = {gaus_hTrackPRes_Eta_N20_N15_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_N20_N15_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_N20_N15_P_10_15->GetParError(2)*100., gaus_hTrackPRes_Eta_N20_N15_P_15_20->GetParError(2)*100.};
  Double_t y_Eta_N15_N10[n] = {gaus_hTrackPRes_Eta_N15_N10_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_N15_N10_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_N15_N10_P_10_15->GetParameter(2)*100., gaus_hTrackPRes_Eta_N15_N10_P_15_20->GetParameter(2)*100.};
  Double_t y_Eta_N15_N10_error[n] = {gaus_hTrackPRes_Eta_N15_N10_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_N15_N10_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_N15_N10_P_10_15->GetParError(2)*100., gaus_hTrackPRes_Eta_N15_N10_P_15_20->GetParError(2)*100.};
  Double_t y_Eta_N10_N05[n] = {gaus_hTrackPRes_Eta_N10_N05_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_N10_N05_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_N10_N05_P_10_15->GetParameter(2)*100., gaus_hTrackPRes_Eta_N10_N05_P_15_20->GetParameter(2)*100.};
  Double_t y_Eta_N10_N05_error[n] = {gaus_hTrackPRes_Eta_N10_N05_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_N10_N05_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_N10_N05_P_10_15->GetParError(2)*100., gaus_hTrackPRes_Eta_N10_N05_P_15_20->GetParError(2)*100.};
  Double_t y_Eta_N05_000[n] = {gaus_hTrackPRes_Eta_N05_000_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_N05_000_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_N05_000_P_10_15->GetParameter(2)*100., gaus_hTrackPRes_Eta_N05_000_P_15_20->GetParameter(2)*100.};
  Double_t y_Eta_N05_000_error[n] = {gaus_hTrackPRes_Eta_N05_000_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_N05_000_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_N05_000_P_10_15->GetParError(2)*100., gaus_hTrackPRes_Eta_N05_000_P_15_20->GetParError(2)*100.};
  Double_t y_Eta_000_P05[n] = {gaus_hTrackPRes_Eta_000_P05_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_000_P05_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_000_P05_P_10_15->GetParameter(2)*100., gaus_hTrackPRes_Eta_000_P05_P_15_20->GetParameter(2)*100.};
  Double_t y_Eta_000_P05_error[n] = {gaus_hTrackPRes_Eta_000_P05_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_000_P05_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_000_P05_P_10_15->GetParError(2)*100., gaus_hTrackPRes_Eta_000_P05_P_15_20->GetParError(2)*100.};
  Double_t y_Eta_P05_P10[n] = {gaus_hTrackPRes_Eta_P05_P10_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_P05_P10_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_P05_P10_P_10_15->GetParameter(2)*100., gaus_hTrackPRes_Eta_P05_P10_P_15_20->GetParameter(2)*100.};
  Double_t y_Eta_P05_P10_error[n] = {gaus_hTrackPRes_Eta_P05_P10_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_P05_P10_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_P05_P10_P_10_15->GetParError(2)*100., gaus_hTrackPRes_Eta_P05_P10_P_15_20->GetParError(2)*100.};
  Double_t y_Eta_P10_P15[n] = {gaus_hTrackPRes_Eta_P10_P15_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_P10_P15_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_P10_P15_P_10_15->GetParameter(2)*100., gaus_hTrackPRes_Eta_P10_P15_P_15_20->GetParameter(2)*100.};
  Double_t y_Eta_P10_P15_error[n] = {gaus_hTrackPRes_Eta_P10_P15_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_P10_P15_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_P10_P15_P_10_15->GetParError(2)*100., gaus_hTrackPRes_Eta_P10_P15_P_15_20->GetParError(2)*100.};
  Double_t y_Eta_P15_P20[n] = {gaus_hTrackPRes_Eta_P15_P20_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_P15_P20_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_P15_P20_P_10_15->GetParameter(2)*100., gaus_hTrackPRes_Eta_P15_P20_P_15_20->GetParameter(2)*100.};
  Double_t y_Eta_P15_P20_error[n] = {gaus_hTrackPRes_Eta_P15_P20_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_P15_P20_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_P15_P20_P_10_15->GetParError(2)*100., gaus_hTrackPRes_Eta_P15_P20_P_15_20->GetParError(2)*100.};
  Double_t y_Eta_P20_P25[n] = {gaus_hTrackPRes_Eta_P20_P25_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_P20_P25_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_P20_P25_P_10_15->GetParameter(2)*100., gaus_hTrackPRes_Eta_P20_P25_P_15_20->GetParameter(2)*100.};
  Double_t y_Eta_P20_P25_error[n] = {gaus_hTrackPRes_Eta_P20_P25_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_P20_P25_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_P20_P25_P_10_15->GetParError(2)*100., gaus_hTrackPRes_Eta_P20_P25_P_15_20->GetParError(2)*100.};
  Double_t y_Eta_P25_P30[n] = {gaus_hTrackPRes_Eta_P25_P30_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_P25_P30_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_P25_P30_P_10_15->GetParameter(2)*100., gaus_hTrackPRes_Eta_P25_P30_P_15_20->GetParameter(2)*100.};
  Double_t y_Eta_P25_P30_error[n] = {gaus_hTrackPRes_Eta_P25_P30_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_P25_P30_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_P25_P30_P_10_15->GetParError(2)*100., gaus_hTrackPRes_Eta_P25_P30_P_15_20->GetParError(2)*100.};
  Double_t y_Eta_P30_P35[n] = {gaus_hTrackPRes_Eta_P30_P35_P_00_05->GetParameter(2)*100., gaus_hTrackPRes_Eta_P30_P35_P_05_10->GetParameter(2)*100.,
			       gaus_hTrackPRes_Eta_P30_P35_P_10_15->GetParameter(2)*100., gaus_hTrackPRes_Eta_P30_P35_P_15_20->GetParameter(2)*100.};
  Double_t y_Eta_P30_P35_error[n] = {gaus_hTrackPRes_Eta_P30_P35_P_00_05->GetParError(2)*100., gaus_hTrackPRes_Eta_P30_P35_P_05_10->GetParError(2)*100.,
				     gaus_hTrackPRes_Eta_P30_P35_P_10_15->GetParError(2)*100., gaus_hTrackPRes_Eta_P30_P35_P_15_20->GetParError(2)*100.};

  TGraphErrors* gPvsPRes_Eta_N35_N30 = new TGraphErrors(n, x, y_Eta_N35_N30, x_error, y_Eta_N35_N30_error);
  TGraphErrors* gPvsPRes_Eta_N30_N25 = new TGraphErrors(n, x, y_Eta_N30_N25, x_error, y_Eta_N30_N25_error);
  TGraphErrors* gPvsPRes_Eta_N25_N20 = new TGraphErrors(n, x, y_Eta_N25_N20, x_error, y_Eta_N25_N20_error);
  TGraphErrors* gPvsPRes_Eta_N20_N15 = new TGraphErrors(n, x, y_Eta_N20_N15, x_error, y_Eta_N20_N15_error);
  TGraphErrors* gPvsPRes_Eta_N15_N10 = new TGraphErrors(n, x, y_Eta_N15_N10, x_error, y_Eta_N15_N10_error);
  TGraphErrors* gPvsPRes_Eta_N10_N05 = new TGraphErrors(n, x, y_Eta_N10_N05, x_error, y_Eta_N10_N05_error);
  TGraphErrors* gPvsPRes_Eta_N05_000 = new TGraphErrors(n, x, y_Eta_N05_000, x_error, y_Eta_N05_000_error);
  TGraphErrors* gPvsPRes_Eta_000_P05 = new TGraphErrors(n, x, y_Eta_000_P05, x_error, y_Eta_000_P05_error);
  TGraphErrors* gPvsPRes_Eta_P05_P10 = new TGraphErrors(n, x, y_Eta_P05_P10, x_error, y_Eta_P05_P10_error);
  TGraphErrors* gPvsPRes_Eta_P10_P15 = new TGraphErrors(n, x, y_Eta_P10_P15, x_error, y_Eta_P10_P15_error);
  TGraphErrors* gPvsPRes_Eta_P15_P20 = new TGraphErrors(n, x, y_Eta_P15_P20, x_error, y_Eta_P15_P20_error);
  TGraphErrors* gPvsPRes_Eta_P20_P25 = new TGraphErrors(n, x, y_Eta_P20_P25, x_error, y_Eta_P20_P25_error);
  TGraphErrors* gPvsPRes_Eta_P25_P30 = new TGraphErrors(n, x, y_Eta_P25_P30, x_error, y_Eta_P25_P30_error);
  TGraphErrors* gPvsPRes_Eta_P30_P35 = new TGraphErrors(n, x, y_Eta_P30_P35, x_error, y_Eta_P30_P35_error);
 
  // ePIC tracking resolution
  const int m = 7;
  Double_t x_epic[m] = {0.5, 1.0, 2.0, 5.0, 10.0, 15., 20.};
  Double_t x_epic_error[m] = {0., 0., 0., 0., 0., 0., 0.};
  
  Double_t y_epic_Eta_N35_N30[m] = {8.20643, 7.56995, 7.69945, 7.50791, 8.96818, 10.8031, 11.9965};
  Double_t y_epic_Eta_N35_N30_error[m] = {0.124932, 0.108154, 0.116841, 0.113985, 0.145546, 0.206098, 0.262865};
  Double_t y_epic_Eta_N30_N25[m] = {6.18447, 5.18337, 4.59104, 3.81117, 3.8943, 4.70526, 5.18833};
  Double_t y_epic_Eta_N30_N25_error[m] = {0.0954384, 0.0738454, 0.0652719, 0.0535285, 0.0570558, 0.0730889, 0.0781465};
  Double_t y_epic_Eta_N25_N20[m] = {5.23584, 3.79213, 2.85988, 2.10586, 1.90826, 2.07517, 2.27234};
  Double_t y_epic_Eta_N25_N20_error[m] = {0.084421, 0.0600886, 0.0414147, 0.0314495, 0.0265339, 0.0289302, 0.0311505};
  Double_t y_epic_Eta_N20_N15[m] = {7.22785, 6.15815, 2.44739, 1.48566, 1.20496, 1.17855, 1.2609};
  Double_t y_epic_Eta_N20_N15_error[m] = {0.12228, 0.129706, 0.043368, 0.0225332, 0.0172614, 0.0162043, 0.0183713};
  Double_t y_epic_Eta_N15_N10[m] = {1.97556, 1.35238, 1.20394, 1.09134, 1.05168, 1.08447, 1.21477};
  Double_t y_epic_Eta_N15_N10_error[m] = {0.0402112, 0.0223593, 0.0182101, 0.0150474, 0.0154879, 0.0160599, 0.0179248};
  Double_t y_epic_Eta_N10_N05[m] = {0.524789, 0.525903, 0.541316, 0.649561, 0.781043, 0.880509, 1.05389};
  Double_t y_epic_Eta_N10_N05_error[m] = {0.00755207, 0.00714204, 0.00760855, 0.0092659, 0.0114611, 0.0129635, 0.0152524};
  Double_t y_epic_Eta_N05_000[m] = {0.421372, 0.416316, 0.444236, 0.535543, 0.716816, 0.909378, 1.10203};
  Double_t y_epic_Eta_N05_000_error[m] = {0.00713129, 0.00657915, 0.0053638, 0.00688407, 0.010973, 0.013214, 0.0162811};
  Double_t y_epic_Eta_000_P05[m] = {0.433855, 0.428103, 0.446722, 0.536589, 0.685991, 0.921926, 1.11138};
  Double_t y_epic_Eta_000_P05_error[m] = {0.00573178, 0.00691924, 0.00559723, 0.00891355, 0.00980929, 0.0135328, 0.0163765};
  Double_t y_epic_Eta_P05_P10[m] = {0.532086, 0.521685, 0.547007, 0.637458, 0.77645, 0.929363, 1.12052};
  Double_t y_epic_Eta_P05_P10_error[m] = {0.00829484, 0.00750643, 0.00794825, 0.00895842, 0.0117291, 0.0128793, 0.0158172};
  Double_t y_epic_Eta_P10_P15[m] = {4.20079, 1.32766, 1.21846, 1.10273, 1.06113, 1.13194, 1.23432};
  Double_t y_epic_Eta_P10_P15_error[m] = {0.14938, 0.0222423, 0.017656, 0.015621, 0.0153349, 0.0156933, 0.0172399};
  Double_t y_epic_Eta_P15_P20[m] = {7.48416, 6.20959, 2.00866, 1.39446, 1.16597, 1.15743, 1.20808};
  Double_t y_epic_Eta_P15_P20_error[m] = {0.136198, 0.147771, 0.032618, 0.0206594, 0.0163857, 0.0167934, 0.0166522};
  Double_t y_epic_Eta_P20_P25[m] = {4.31769, 2.10402, 2.16885, 1.79735, 1.57311, 1.58057, 1.61354};
  Double_t y_epic_Eta_P20_P25_error[m] = {0.113272, 0.0331813, 0.031794, 0.0249298, 0.0214576, 0.0219766, 0.0234497};
  Double_t y_epic_Eta_P25_P30[m] = {4.16891, 3.43897, 3.24534, 2.99746, 2.73549, 2.96044, 3.19435};
  Double_t y_epic_Eta_P25_P30_error[m] = {0.0672504, 0.0538235, 0.0472295, 0.0429689, 0.0392788, 0.0419626, 0.0482197};
  Double_t y_epic_Eta_P30_P35[m] = {6.74298, 5.88147, 5.49895, 5.65444, 5.73812, 6.42789, 7.23576};
  Double_t y_epic_Eta_P30_P35_error[m] = {0.108379, 0.0855475, 0.077388, 0.0830848, 0.0841084, 0.0976485, 0.110392};
 
  TF1* func = new TF1("func", "pol1", 0.8, 21.0);
  TF1* func1 = new TF1("func1", "pol1", 0.8, 21.0);
  TF1* func2 = new TF1("func2", "pol1", 0.8, 21.0);
  TF1* func3 = new TF1("func3", "pol1", 0.8, 21.0);
  TF1* func4 = new TF1("func4", "pol1", 0.8, 21.0);
  TF1* func5 = new TF1("func5", "pol1", 0.8, 21.0);
  TF1* func6 = new TF1("func6", "pol1", 0.8, 21.0);
  TF1* func7 = new TF1("func7", "pol1", 0.8, 21.0);
  TF1* func8 = new TF1("func8", "pol1", 0.8, 21.0);
  TF1* func9 = new TF1("func9", "pol1", 0.8, 21.0);
  TF1* func10 = new TF1("func10", "pol1", 0.8, 21.0);
  TF1* func11 = new TF1("func11", "pol1", 0.8, 21.0);
  TF1* func12 = new TF1("func12", "pol1", 0.8, 21.0);
  TF1* func13 = new TF1("func13", "pol1", 0.8, 21.0);
  func->SetLineColor(kRed); func->SetLineWidth(2); func->SetLineStyle(2);
  func1->SetLineColor(kRed); func1->SetLineWidth(2); func1->SetLineStyle(2);
  func2->SetLineColor(kRed); func2->SetLineWidth(2); func2->SetLineStyle(2);
  func3->SetLineColor(kRed); func3->SetLineWidth(2); func3->SetLineStyle(2);
  func4->SetLineColor(kRed); func4->SetLineWidth(2); func4->SetLineStyle(2);
  func5->SetLineColor(kRed); func5->SetLineWidth(2); func5->SetLineStyle(2);
  func6->SetLineColor(kRed); func6->SetLineWidth(2); func6->SetLineStyle(2);
  func7->SetLineColor(kRed); func7->SetLineWidth(2); func7->SetLineStyle(2);
  func8->SetLineColor(kRed); func8->SetLineWidth(2); func8->SetLineStyle(2);
  func9->SetLineColor(kRed); func9->SetLineWidth(2); func9->SetLineStyle(2);
  func10->SetLineColor(kRed); func10->SetLineWidth(2); func10->SetLineStyle(2);
  func11->SetLineColor(kRed); func11->SetLineWidth(2); func11->SetLineStyle(2);
  func12->SetLineColor(kRed); func12->SetLineWidth(2); func12->SetLineStyle(2);
  func13->SetLineColor(kRed); func13->SetLineWidth(2); func13->SetLineStyle(2);

  TGraphErrors* gPvsPRes_Eta_N35_N30_epic = new TGraphErrors(m, x_epic, y_epic_Eta_N35_N30, x_epic_error, y_epic_Eta_N35_N30_error);
  TGraphErrors* gPvsPRes_Eta_N30_N25_epic = new TGraphErrors(m, x_epic, y_epic_Eta_N30_N25, x_epic_error, y_epic_Eta_N30_N25_error);
  TGraphErrors* gPvsPRes_Eta_N25_N20_epic = new TGraphErrors(m, x_epic, y_epic_Eta_N25_N20, x_epic_error, y_epic_Eta_N25_N20_error);
  TGraphErrors* gPvsPRes_Eta_N20_N15_epic = new TGraphErrors(m, x_epic, y_epic_Eta_N20_N15, x_epic_error, y_epic_Eta_N20_N15_error);
  TGraphErrors* gPvsPRes_Eta_N15_N10_epic = new TGraphErrors(m, x_epic, y_epic_Eta_N15_N10, x_epic_error, y_epic_Eta_N15_N10_error);
  TGraphErrors* gPvsPRes_Eta_N10_N05_epic = new TGraphErrors(m, x_epic, y_epic_Eta_N10_N05, x_epic_error, y_epic_Eta_N10_N05_error);
  TGraphErrors* gPvsPRes_Eta_N05_000_epic = new TGraphErrors(m, x_epic, y_epic_Eta_N05_000, x_epic_error, y_epic_Eta_N05_000_error);
  TGraphErrors* gPvsPRes_Eta_000_P05_epic = new TGraphErrors(m, x_epic, y_epic_Eta_000_P05, x_epic_error, y_epic_Eta_000_P05_error);
  TGraphErrors* gPvsPRes_Eta_P05_P10_epic = new TGraphErrors(m, x_epic, y_epic_Eta_P05_P10, x_epic_error, y_epic_Eta_P05_P10_error);
  TGraphErrors* gPvsPRes_Eta_P10_P15_epic = new TGraphErrors(m, x_epic, y_epic_Eta_P10_P15, x_epic_error, y_epic_Eta_P10_P15_error);
  TGraphErrors* gPvsPRes_Eta_P15_P20_epic = new TGraphErrors(m, x_epic, y_epic_Eta_P15_P20, x_epic_error, y_epic_Eta_P15_P20_error);
  TGraphErrors* gPvsPRes_Eta_P20_P25_epic = new TGraphErrors(m, x_epic, y_epic_Eta_P20_P25, x_epic_error, y_epic_Eta_P20_P25_error);
  TGraphErrors* gPvsPRes_Eta_P25_P30_epic = new TGraphErrors(m, x_epic, y_epic_Eta_P25_P30, x_epic_error, y_epic_Eta_P25_P30_error);
  TGraphErrors* gPvsPRes_Eta_P30_P35_epic = new TGraphErrors(m, x_epic, y_epic_Eta_P30_P35, x_epic_error, y_epic_Eta_P30_P35_error);

  TCanvas *cnv5 = new TCanvas("cnv5", "cnv5");
  gPvsPRes_Eta_N35_N30_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_N35_N30_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_N35_N30_epic->SetTitle("-3.5<#eta<-3.0");
  gPvsPRes_Eta_N35_N30_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_N35_N30_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_N35_N30_epic->GetYaxis()->SetRangeUser(0.0, 13.0);
  gPvsPRes_Eta_N35_N30_epic->SetLineColor(kRed);
  gPvsPRes_Eta_N35_N30_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_N35_N30_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_N35_N30_epic->Draw("AEP");
  func->SetParameter(0, 7.4); 
  func->SetParameter(1, 0.492);
  func->Draw("SAME"); 
  gPvsPRes_Eta_N35_N30->SetMarkerStyle(21);
  gPvsPRes_Eta_N35_N30->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_N35_N30 = gPvsPRes_Eta_N35_N30->GetFunction("pol1");
  fit_gPvsPRes_Eta_N35_N30->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_N35_N30->SetLineStyle(2);
  gPvsPRes_Eta_N35_N30->Draw("EP, SAME");
  cnv5->SaveAs("./plots/hTrackPRes_N35_N30.png");

  TCanvas *cnv6 = new TCanvas("cnv6", "cnv6");
  gPvsPRes_Eta_N30_N25_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_N30_N25_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_N30_N25_epic->SetTitle("-3.0<#eta<-2.5");
  gPvsPRes_Eta_N30_N25_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_N30_N25_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_N30_N25_epic->GetYaxis()->SetRangeUser(0.0,13.0);
  gPvsPRes_Eta_N30_N25_epic->SetLineColor(kRed);
  gPvsPRes_Eta_N30_N25_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_N30_N25_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_N30_N25_epic->Draw("AEP");
  func1->SetParameter(0, 4.0); 
  func1->SetParameter(1, 0.155);
  func1->Draw("SAME"); 
  gPvsPRes_Eta_N30_N25->SetMarkerStyle(21);
  gPvsPRes_Eta_N30_N25->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_N30_N25 = gPvsPRes_Eta_N30_N25->GetFunction("pol1");
  fit_gPvsPRes_Eta_N30_N25->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_N30_N25->SetLineStyle(2);
  gPvsPRes_Eta_N30_N25->Draw("EP, SAME");
  cnv6->SaveAs("./plots/hTrackPRes_N30_N25.png");
 
  TCanvas *cnv7 = new TCanvas("cnv7", "cnv7");
  gPvsPRes_Eta_N25_N20_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_N25_N20_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_N25_N20_epic->SetTitle("-2.5<#eta<-2.0");
  gPvsPRes_Eta_N25_N20_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_N25_N20_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_N25_N20_epic->GetYaxis()->SetRangeUser(0.0,13.0);
  gPvsPRes_Eta_N25_N20_epic->SetLineColor(kRed);
  gPvsPRes_Eta_N25_N20_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_N25_N20_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_N25_N20_epic->Draw("AEP");
  func2->SetParameter(0, 2.3); 
  func2->SetParameter(1, 0.000);
  func2->Draw("SAME"); 
  gPvsPRes_Eta_N25_N20->SetMarkerStyle(21);
  gPvsPRes_Eta_N25_N20->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_N25_N20 = gPvsPRes_Eta_N25_N20->GetFunction("pol1");
  fit_gPvsPRes_Eta_N25_N20->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_N25_N20->SetLineStyle(2);
  gPvsPRes_Eta_N25_N20->Draw("EP, SAME");
  cnv7->SaveAs("./plots/hTrackPRes_N25_N20.png");

  TCanvas *cnv8 = new TCanvas("cnv8", "cnv8");
  gPvsPRes_Eta_N20_N15_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_N20_N15_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_N20_N15_epic->SetTitle("-2.0<#eta<-1.5");
  gPvsPRes_Eta_N20_N15_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_N20_N15_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_N20_N15_epic->GetYaxis()->SetRangeUser(0.0, 13.0);
  gPvsPRes_Eta_N20_N15_epic->SetLineColor(kRed);
  gPvsPRes_Eta_N20_N15_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_N20_N15_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_N20_N15_epic->Draw("AEP");
  func3->SetParameter(0, 1.4); 
  func3->SetParameter(1, 0.000);
  func3->Draw("SAME");
  gPvsPRes_Eta_N20_N15->SetMarkerStyle(21);
  gPvsPRes_Eta_N20_N15->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_N20_N15 = gPvsPRes_Eta_N20_N15->GetFunction("pol1");
  fit_gPvsPRes_Eta_N20_N15->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_N20_N15->SetLineStyle(2);
  gPvsPRes_Eta_N20_N15->Draw("EP, SAME");
  cnv8->SaveAs("./plots/hTrackPRes_N20_N15.png");

  TCanvas *cnv9 = new TCanvas("cnv9", "cnv9");
  gPvsPRes_Eta_N15_N10_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_N15_N10_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_N15_N10_epic->SetTitle("-1.5<#eta<-1.0");
  gPvsPRes_Eta_N15_N10_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_N15_N10_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_N15_N10_epic->GetYaxis()->SetRangeUser(0.0,13.0);
  gPvsPRes_Eta_N15_N10_epic->SetLineColor(kRed);
  gPvsPRes_Eta_N15_N10_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_N15_N10_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_N15_N10_epic->Draw("AEP");
  func4->SetParameter(0, 1.2); 
  func4->SetParameter(1, 0.000);
  func4->Draw("SAME"); 
  gPvsPRes_Eta_N15_N10->SetMarkerStyle(21);
  gPvsPRes_Eta_N15_N10->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_N15_N10 = gPvsPRes_Eta_N15_N10->GetFunction("pol1");
  fit_gPvsPRes_Eta_N15_N10->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_N15_N10->SetLineStyle(2);
  gPvsPRes_Eta_N15_N10->Draw("EP, SAME");
  cnv9->SaveAs("./plots/hTrackPRes_N15_N10.png");

  TCanvas *cnv10 = new TCanvas("cnv10", "cnv10");
  gPvsPRes_Eta_N10_N05_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_N10_N05_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_N10_N05_epic->SetTitle("-1.0<#eta<-0.5");
  gPvsPRes_Eta_N10_N05_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_N10_N05_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_N10_N05_epic->GetYaxis()->SetRangeUser(0.0,13.0);
  gPvsPRes_Eta_N10_N05_epic->SetLineColor(kRed);
  gPvsPRes_Eta_N10_N05_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_N10_N05_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_N10_N05_epic->Draw("AEP");
  func5->SetParameter(0, 0.5); 
  func5->SetParameter(1, 0.048);
  func5->Draw("SAME"); 
  gPvsPRes_Eta_N10_N05->SetMarkerStyle(21);
  gPvsPRes_Eta_N10_N05->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_N10_N05 = gPvsPRes_Eta_N10_N05->GetFunction("pol1");
  fit_gPvsPRes_Eta_N10_N05->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_N10_N05->SetLineStyle(2);
  gPvsPRes_Eta_N10_N05->Draw("EP, SAME");
  cnv10->SaveAs("./plots/hTrackPRes_N10_N05.png");

  TCanvas *cnv11 = new TCanvas("cnv11", "cnv11");
  gPvsPRes_Eta_N05_000_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_N05_000_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_N05_000_epic->SetTitle("-0.5<#eta<0.0");
  gPvsPRes_Eta_N05_000_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_N05_000_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_N05_000_epic->GetYaxis()->SetRangeUser(0.0,13.);
  gPvsPRes_Eta_N05_000_epic->SetLineColor(kRed);
  gPvsPRes_Eta_N05_000_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_N05_000_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_N05_000_epic->Draw("AEP");
  func6->SetParameter(0, 0.4); 
  func6->SetParameter(1, 0.053);
  func6->Draw("SAME"); 
  gPvsPRes_Eta_N05_000->SetMarkerStyle(21);
  gPvsPRes_Eta_N05_000->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_N05_000 = gPvsPRes_Eta_N05_000->GetFunction("pol1");
  fit_gPvsPRes_Eta_N05_000->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_N05_000->SetLineStyle(2);
  gPvsPRes_Eta_N05_000->Draw("EP, SAME");
  cnv11->SaveAs("./plots/hTrackPRes_N05_000.png");

  TCanvas *cnv12 = new TCanvas("cnv12", "cnv12");
  gPvsPRes_Eta_000_P05_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_000_P05_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_000_P05_epic->SetTitle("0.0<#eta<0.5");
  gPvsPRes_Eta_000_P05_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_000_P05_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_000_P05_epic->GetYaxis()->SetRangeUser(0.,13.);
  gPvsPRes_Eta_000_P05_epic->SetLineColor(kRed);
  gPvsPRes_Eta_000_P05_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_000_P05_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_000_P05_epic->Draw("AEP");
  func7->SetParameter(0, 0.4); 
  func7->SetParameter(1, 0.053);
  func7->Draw("SAME"); 
  gPvsPRes_Eta_000_P05->SetMarkerStyle(21);
  gPvsPRes_Eta_000_P05->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_000_P05 = gPvsPRes_Eta_000_P05->GetFunction("pol1");
  fit_gPvsPRes_Eta_000_P05->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_000_P05->SetLineStyle(2);
  gPvsPRes_Eta_000_P05->Draw("EP, SAME");
  cnv12->SaveAs("./plots/hTrackPRes_000_P05.png");

  TCanvas *cnv13 = new TCanvas("cnv13", "cnv13");
  gPvsPRes_Eta_P05_P10_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_P05_P10_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_P05_P10_epic->SetTitle("0.5<#eta<1.0");
  gPvsPRes_Eta_P05_P10_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_P05_P10_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_P05_P10_epic->GetYaxis()->SetRangeUser(0.,13.);
  gPvsPRes_Eta_P05_P10_epic->SetLineColor(kRed);
  gPvsPRes_Eta_P05_P10_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_P05_P10_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_P05_P10_epic->Draw("AEP");
  func8->SetParameter(0, 0.5);
  func8->SetParameter(1, 0.051);
  func8->Draw("SAME"); 
  gPvsPRes_Eta_P05_P10->SetMarkerStyle(21);
  gPvsPRes_Eta_P05_P10->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_P05_P10 = gPvsPRes_Eta_P05_P10->GetFunction("pol1");
  fit_gPvsPRes_Eta_P05_P10->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_P05_P10->SetLineStyle(2);
  gPvsPRes_Eta_P05_P10->Draw("EP, SAME");
  cnv13->SaveAs("./plots/hTrackPRes_P05_P10.png");

  TCanvas *cnv14 = new TCanvas("cnv14", "cnv14");
  gPvsPRes_Eta_P10_P15_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_P10_P15_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_P10_P15_epic->SetTitle("1.0<#eta<1.5");
  gPvsPRes_Eta_P10_P15_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_P10_P15_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_P10_P15_epic->GetYaxis()->SetRangeUser(0.,13.);
  gPvsPRes_Eta_P10_P15_epic->SetLineColor(kRed);
  gPvsPRes_Eta_P10_P15_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_P10_P15_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_P10_P15_epic->Draw("AEP");
  func9->SetParameter(0, 1.2); //1.176); 
  func9->SetParameter(1, 0.007); //-0.002);
  func9->Draw("SAME"); 
  gPvsPRes_Eta_P10_P15->SetMarkerStyle(21);
  gPvsPRes_Eta_P10_P15->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_P10_P15 = gPvsPRes_Eta_P10_P15->GetFunction("pol1");
  fit_gPvsPRes_Eta_P10_P15->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_P10_P15->SetLineStyle(2);
  gPvsPRes_Eta_P10_P15->Draw("EP, SAME");
  cnv14->SaveAs("./plots/hTrackPRes_P10_P15.png");

  TCanvas *cnv15 = new TCanvas("cnv15", "cnv15");
  gPvsPRes_Eta_P15_P20_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_P15_P20_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_P15_P20_epic->SetTitle("1.5<#eta<2.0");
  gPvsPRes_Eta_P15_P20_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_P15_P20_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_P15_P20_epic->GetYaxis()->SetRangeUser(0.,13.);
  gPvsPRes_Eta_P15_P20_epic->SetLineColor(kRed);
  gPvsPRes_Eta_P15_P20_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_P15_P20_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_P15_P20_epic->Draw("AEP");
  func10->SetParameter(0, 1.3);
  func10->SetParameter(1, 0.000);
  func10->Draw("SAME"); 
  gPvsPRes_Eta_P15_P20->SetMarkerStyle(21);
  gPvsPRes_Eta_P15_P20->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_P15_P20 = gPvsPRes_Eta_P15_P20->GetFunction("pol1");
  fit_gPvsPRes_Eta_P15_P20->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_P15_P20->SetLineStyle(2);
  gPvsPRes_Eta_P15_P20->Draw("EP, SAME");
  cnv15->SaveAs("./plots/hTrackPRes_P15_P20.png");

  TCanvas *cnv16 = new TCanvas("cnv16", "cnv16");
  gPvsPRes_Eta_P20_P25_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_P20_P25_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_P20_P25_epic->SetTitle("2.0<#eta<2.5");
  gPvsPRes_Eta_P20_P25_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_P20_P25_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_P20_P25_epic->GetYaxis()->SetRangeUser(0.,13.);
  gPvsPRes_Eta_P20_P25_epic->SetLineColor(kRed);
  gPvsPRes_Eta_P20_P25_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_P20_P25_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_P20_P25_epic->Draw("AEP");
  func11->SetParameter(0, 1.8); 
  func11->SetParameter(1, 0.000);
  func11->Draw("SAME"); 
  gPvsPRes_Eta_P20_P25->SetMarkerStyle(21);
  gPvsPRes_Eta_P20_P25->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_P20_P25 = gPvsPRes_Eta_P20_P25->GetFunction("pol1");
  fit_gPvsPRes_Eta_P20_P25->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_P20_P25->SetLineStyle(2);
  gPvsPRes_Eta_P20_P25->Draw("EP, SAME");
  cnv16->SaveAs("./plots/hTrackPRes_P20_P25.png");

  TCanvas *cnv17 = new TCanvas("cnv17", "cnv17");
  gPvsPRes_Eta_P25_P30_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_P25_P30_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_P25_P30_epic->SetTitle("2.5<#eta<3.0");
  gPvsPRes_Eta_P25_P30_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_P25_P30_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_P25_P30_epic->GetYaxis()->SetRangeUser(0.,13.);
  gPvsPRes_Eta_P25_P30_epic->SetLineColor(kRed);
  gPvsPRes_Eta_P25_P30_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_P25_P30_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_P25_P30_epic->Draw("AEP");
  func12->SetParameter(0, 3.0);
  func12->SetParameter(1, 0.037);
  func12->Draw("SAME"); 
  gPvsPRes_Eta_P25_P30->SetMarkerStyle(21);
  gPvsPRes_Eta_P25_P30->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_P25_P30 = gPvsPRes_Eta_P25_P30->GetFunction("pol1");
  fit_gPvsPRes_Eta_P25_P30->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_P25_P30->SetLineStyle(2);
  gPvsPRes_Eta_P25_P30->Draw("EP, SAME");
  cnv17->SaveAs("./plots/hTrackPRes_P25_P30.png");

  TCanvas *cnv18 = new TCanvas("cnv18", "cnv18");
  gPvsPRes_Eta_P30_P35_epic->GetXaxis()->CenterTitle(true);
  gPvsPRes_Eta_P30_P35_epic->GetYaxis()->CenterTitle(true);
  gPvsPRes_Eta_P30_P35_epic->SetTitle("3.0<#eta<3.5");
  gPvsPRes_Eta_P30_P35_epic->GetXaxis()->SetTitle("p [GeV/c]");
  gPvsPRes_Eta_P30_P35_epic->GetYaxis()->SetTitle("#Delta p/p [%]");
  gPvsPRes_Eta_P30_P35_epic->GetYaxis()->SetRangeUser(0.,13.);
  gPvsPRes_Eta_P30_P35_epic->SetLineColor(kRed);
  gPvsPRes_Eta_P30_P35_epic->SetMarkerStyle(21);
  gPvsPRes_Eta_P30_P35_epic->SetMarkerColor(kRed);
  gPvsPRes_Eta_P30_P35_epic->Draw("AEP");
  func13->SetParameter(0, 5.4); 
  func13->SetParameter(1, 0.231);
  func13->Draw("SAME"); 
  gPvsPRes_Eta_P30_P35->SetMarkerStyle(21);
  gPvsPRes_Eta_P30_P35->Fit("pol1","R","", 0.8, 21.);
  TF1 * fit_gPvsPRes_Eta_P30_P35 = gPvsPRes_Eta_P30_P35->GetFunction("pol1");
  fit_gPvsPRes_Eta_P30_P35->SetLineColor(kBlack);
  fit_gPvsPRes_Eta_P30_P35->SetLineStyle(2);
  gPvsPRes_Eta_P30_P35->Draw("EP, SAME");
  cnv18->SaveAs("./plots/hTrackPRes_P30_P35.png");

  cout << "** Done..." << endl;

  delete treeReader;
  delete chain;
}
