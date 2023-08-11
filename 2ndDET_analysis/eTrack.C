//////////////////////////////////////////
// 08/02/2023 Jihee Kim (jkim11@bnl.gov)
// Tracks of particles
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

//-------------------------------------------------------------------------------------------//
//define constants, pt and eta bins
const int GEN=0, REC=1;  //generated particle, reconstructed particle
const int nPBin=10;
const int nPtBin=10;
const int nEtaBin=8;

string tagGenRec[2]={"Gen","Rec"};
const float p[nPBin+1]={0.,0.5,1,1.5,2,3,4,5,10,15,20};
const float pt[nPtBin+1]={0.,0.5,1,1.5,2,3,4,5,10,15,20};
const float eta[nEtaBin+1]={-4,-3,-2,-1,0,1,2,3,4};

//-------------------------------------------------------------------------------------------//
//define functions
void SetStyle(TStyle* kStyle);
void SetHisto(TH1F* h);
int GetPBin(float myP);
int GetPtBin(float myPt);
int GetEtaBin(float myEta);
TF1* GetFit(const char* fname, int ieta, int ip, TH1F* h);

//-------------------------------------------------------------------------------------------//
//Main function
void eTrack(const char *inputFile, const char *outputFile="etrackKin.root")
{
  int i,j;

  // Setting for figures
  TStyle* kStyle = new TStyle("kStyle","Kim's Style");
  SetStyle(kStyle);

  gSystem->Load("libDelphes");

  // Create chain of root trees
  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  // Create object of class ExRootTreeReader
  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);

  // Get pointers to branches used in this analysis
  TClonesArray *branchParticle = treeReader->UseBranch("Particle");
  TClonesArray* branchTrack = treeReader->UseBranch("Track");
  TClonesArray* branchElectron = treeReader->UseBranch("Electron");

  // Get total number of events
  Long64_t allEntries = treeReader->GetEntries();
  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Track *track;  
  //Electron* track;
  
  // Book histograms
  TH2F* h2DKin[2][2]; //[gen,rec][p,pt] 
  for (i=0;i<2;i++) {
    h2DKin[i][0]=new TH2F(Form("h2DKin%sPEvta",tagGenRec[i].c_str()),
			  Form("eta vs p (%s);#eta;p (GeV)",tagGenRec[i].c_str()),
			  90,-4,5,100,0,20);
    h2DKin[i][1]=new TH2F(Form("h2DKin%sPtEvta",tagGenRec[i].c_str()),
                          Form("eta vs pt (%s);#eta;p_{T} (GeV)",tagGenRec[i].c_str()),
                          90,-4,5,100,0,5);

    for (j=0;j<2;j++) h2DKin[i][j]->Sumw2();
  }

  TH1F* hResP[nEtaBin][nPtBin];
  TH1F* hResPt[nEtaBin][nPtBin];
  for (i=0;i<nEtaBin;i++) {
    for (j=0;j<nPBin;j++) {
      hResP[i][j]=new TH1F(Form("hResP_E%dP%d",i,j),
			   Form("res(p) %.1f<#eta<%0.1f, %.1f-%.fGeV;(P^{rec}-P^{gen})/P^{gen};nTrk",eta[i],eta[i+1],p[j],p[j+1]),
			   50,-0.05,0.05);
      hResP[i][j]->Sumw2();
      SetHisto(hResP[i][j]);

      hResPt[i][j]=new TH1F(Form("hResPt_E%dP%d",i,j),
			    Form("res(pt) %.1f<#eta<%0.1f, %.1f-%.fGeV;(P_{T}^{rec}-P^{gen})/P_{T}^{gen};nTrk",eta[i],eta[i+1],p[j],p[j+1]),
			    50,-0.05,0.05);
      hResPt[i][j]->Sumw2();
      SetHisto(hResPt[i][j]);
    } 
  }
  
  Long64_t entry;
  float pGen,pRec,pxRec,pyRec,pzRec,ptGen,ptRec,etaGen, etaRec;

  //---------------------------
  // Loop over all events
  //---------------------------
  for(entry = 0; entry < allEntries; ++entry)
    {
      // Load selected branches with data from specified event
      treeReader->ReadEntry(entry);
      //
      if (branchTrack->GetEntries() > 0)
	//if (branchElectron->GetEntries() > 0)
	{
	  // Loop on all tracks
	  for (i = 0; i < branchTrack->GetEntries(); ++i)
	    //for (i = 0; i < branchElectron->GetEntries(); ++i)
	    {
	      track = (Track*) branchTrack->At(i);
	      //track = (Electron*) branchElectron->At(i); 
	      pRec = track->P;
	      //pxRec = track->Px;
	      //pyRec = track->Py;
	      //pzRec = track->Pz;
	      ptRec = track->PT;
	      etaRec = track->Eta;
	      h2DKin[REC][0]->Fill(etaRec,pRec);
	      h2DKin[REC][1]->Fill(etaRec,ptRec);

	      particle = (GenParticle*) track->Particle.GetObject();
	      pGen = particle->P;
	      ptGen = particle->PT;
	      etaGen = particle->Eta;
	      h2DKin[GEN][0]->Fill(etaGen,pGen);
	      h2DKin[GEN][1]->Fill(etaGen,ptGen);
	      
	      //cout<<"pGen="<<pGen<<", pRec="<<pRec<<", pxRec="<<pxRec<<", pyRec="<<pyRec<<", pzRec="<<pzRec<<", ptRec="<<ptRec<<endl;
	      
	      int iEta=GetEtaBin(etaGen);
	      int iP=GetPBin(pGen);
	      int iPt=GetPtBin(ptGen);
	      if (iEta<0) continue;

	      if (iP>=0) hResP[iEta][iP]->Fill((pRec-pGen)/pGen);
	      if (iPt>=0) hResPt[iEta][iPt]->Fill((ptRec-ptGen)/ptGen);
	    }
	}
    }

  //---------------------------
  // fits for resolutions   
  //---------------------------
  TF1* fResP[nEtaBin][nPBin];
  TF1* fResPt[nEtaBin][nPtBin];
  for (i=0;i<nEtaBin;i++) {
    for(j=0;j<nPBin;j++) fResP[i][j]=GetFit("P",i,j,hResP[i][j]);
    for(j=0;j<nPtBin;j++) fResPt[i][j]=GetFit("Pt",i,j,hResPt[i][j]);
  }


  //---------------------------
  // Plot figures
  //---------------------------
  TCanvas* c2DKin[2][2]; //[gen, rec][p vs eta, pt vs eta]
  for (i=0;i<2;i++) {
    c2DKin[i][0]=new TCanvas(Form("c%s_p_eta",tagGenRec[i].c_str()),Form("p vs eta (%s)",tagGenRec[i].c_str()),600,600);
    c2DKin[i][1]=new TCanvas(Form("c%s_pt_eta",tagGenRec[i].c_str()),Form("pt vs eta (%s)",tagGenRec[i].c_str()),600,600);
    
    for (j=0;j<2;j++) {
      c2DKin[i][j]->Range(0,0,1,1);
      c2DKin[i][j]->cd();
      gPad->SetLogz();

      h2DKin[i][j]->Draw("colz");
    }

    c2DKin[i][0]->SaveAs(Form("p_eta_%s.pdf",tagGenRec[i].c_str()));
    c2DKin[i][1]->SaveAs(Form("pt_eta_%s.pdf",tagGenRec[i].c_str()));
  }

  TCanvas* cRes[2][nEtaBin];  //[p,pt][nEtaBin]
  for (i=0;i<nEtaBin;i++) {
    cRes[0][i]=new TCanvas(Form("cPResEta%d",i),Form("p resolution, %0.1f<#eta<%0.1f",eta[i],eta[i+1]),1000,600);
    cRes[1][i]=new TCanvas(Form("cPtResEta%d",i),Form("pt resolution, %0.1f<#eta<%0.1f",eta[i],eta[i+1]),1000,600);

    cRes[0][i]->Range(0,0,1,1);
    cRes[0][i]->Divide(ceil(nPBin/2),2);
    for (j=0;j<nPBin;j++) {
      cRes[0][i]->cd(j+1);
      hResP[i][j]->Draw();

      if (fResP[i][j]==0) continue;
      cout<<"par0="<<fResP[i][j]->GetParameter(0)<<endl;
      fResP[i][j]->Draw();
      hResP[i][j]->Draw("same");
    }
    cRes[0][i]->SaveAs(Form("pResolutionEta%d.pdf",i));
    
    cRes[1][i]->Range(0,0,1,1);
    cRes[1][i]->Divide(ceil(nPtBin/2),2);
    for (j=0;j<nPtBin;j++) {
      cRes[1][i]->cd(j+1);
      hResPt[i][j]->Draw();
     
      if (fResPt[i][j]==0) continue;
      fResPt[i][j]->Draw("same");
      hResPt[i][j]->Draw("same");
    }  
    cRes[1][i]->SaveAs(Form("ptResolutionEta%d.pdf",i));
  }

  //--------------------------- 
  // Output file to 
  // save all figures
  //--------------------------- 
  TFile* output = new TFile(outputFile,"RECREATE");
  for (i=0;i<2;i++) {
    for (j=0;j<2;j++) {
      c2DKin[i][j]->Write();
      h2DKin[i][j]->Write();
    }
  }

  for (i=0;i<nEtaBin;i++) {
    //for (j=0;j<2;j++) cRes[j][i]->Write();
    for (j=0;j<nPBin;j++) {
      hResP[i][j]->Write();
      //if (fResP[i][j]) fResP[i][j]->Write();
    }
    for (j=0;j<nPtBin;j++) {
      hResPt[i][j]->Write();
      //if (fResPt[i][j]) fResPt[i][j]->Write();
    }
  }
  
  output->Close();
  delete treeReader;
  delete chain;
}
//-------------------------------------------------------------------------------------------//
void SetStyle(TStyle* kStyle) 
{
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

  gStyle->SetPalette(1);
}
//-------------------------------------------------------------------------------------------// 
int GetEtaBin(float myEta)
{
  int i;
  if (myEta<eta[0] || myEta>eta[nEtaBin]) return -1;

  for (i=0;i<nEtaBin;i++) {
    if (myEta>=eta[i] && myEta<eta[i+1]) return i;
  }

  return -1;
}
//-------------------------------------------------------------------------------------------//
int GetPBin(float myP)
{
  int i;
  if (myP<p[0] || myP>p[nPBin]) return -1;

  for (i=0;i<nPBin;i++) {
    if (myP>=p[i] && myP<p[i+1]) return i;
  }

  return -1;
}
//-------------------------------------------------------------------------------------------//
int GetPtBin(float myPt)
{
  int i;
  if (myPt<pt[0] || myPt>pt[nPtBin]) return -1;

  for (i=0;i<nPtBin;i++) {
    if (myPt>=pt[i] && myPt<pt[i+1]) return i;
  }

  return -1;
}
//-------------------------------------------------------------------------------------------//
void SetHisto(TH1F* h)
{
  h->SetStats(0);
  h->SetLineColor(4);
  h->SetMarkerColor(4);
  h->SetMarkerStyle(24);

  h->GetXaxis()->SetNdivisions(505);
  h->GetXaxis()->SetTitleFont(43);
  h->GetXaxis()->SetTitleSize(18);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(14);
  
  h->GetYaxis()->SetNdivisions(505);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(18);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(14);
}
//-------------------------------------------------------------------------------------------//
TF1* GetFit(const char* fname, int ieta, int ip, TH1F* h)
{
  if (h->GetEntries()<20) return 0;

  TF1* f=new TF1(Form("fitRes_Eta%d_%s%d",ieta,fname,ip),"[0]*exp(-((x-[1])^2)/(2*[2]^2))",-0.05,0.05);
  f->SetLineStyle(7);
  f->SetLineColor(2);

  f->SetParameter(0,h->GetMaximum());
  f->SetParameter(1,0);
  f->SetParameter(2,0.001);

  h->Fit(f,"N");
  return f;
}
