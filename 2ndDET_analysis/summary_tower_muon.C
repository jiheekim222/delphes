//////////////////////////////////////////
// 08/17/2023 Jihee Kim (jkim11@bnl.gov)
// Energy resolutions
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

void singleMuon(const char *inputFile)
{
  // Setting for figures
  TStyle* kStyle = new TStyle("kStyle","Kim's Style");
  kStyle->SetOptStat(0);
  kStyle->SetOptTitle(0);
  kStyle->SetOptFit(0);
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
  TClonesArray* branchTower = treeReader->UseBranch("Tower");
  // Get total number of events
  Long64_t allEntries = treeReader->GetEntries();
  cout << "** Chain contains " << allEntries << " events" << endl;

  GenParticle *particle;
  Tower *tower;  
  
  // Book histograms
  TH1* hTowerERes_Eta_N10_P10_P_00_05 = new TH1D("hTowerERes_Eta_N10_P10_P_00_05","0<P[GeV]<=5 & -1.0<=#eta<=1.0;(E^{rec}-E^{gen})/E^{gen};Events",100, -0.2, 0.2); 
  TH1* hTowerERes_Eta_N10_P10_P_05_10 = new TH1D("hTowerERes_Eta_N10_P10_P_05_10","5<P[GeV]<=10 & -1.0<=#eta<=1.0;(E^{rec}-E^{gen})/E^{gen};Events",100, -0.1, 0.1); 
  TH1* hTowerERes_Eta_N10_P10_P_10_15 = new TH1D("hTowerERes_Eta_N10_P10_P_10_15","10<P[GeV]<=15 & -1.0<=#eta<=1.0;(E^{rec}-E^{gen})/E^{gen};Events",100, -0.1, 0.1); 
  TH1* hTowerERes_Eta_N10_P10_P_15_20 = new TH1D("hTowerERes_Eta_N10_P10_P_15_20","15<P[GeV]<=20 & -1.0<=#eta<=1.0;(E^{rec}-E^{gen})/E^{gen};Events",100, -0.1, 0.1); 

  Long64_t entry;
  Int_t i;
  Int_t n_miss_branch = 0;

  // Loop over all events
  for(entry = 0; entry < allEntries; ++entry)
  {
    // Load selected branches with data from specified event
    treeReader->ReadEntry(entry);
    // Check tower branch
    if (branchTower->GetEntries() > 0)
    {
      // Loop on all towers
      for (i = 0; i < branchTower->GetEntries(); ++i)
      {
        tower = (Tower*) branchTower->At(i);
	// Loop over all generated particles which are associated with each tower
        for (int j=0; j< tower->Particles.GetEntriesFast(); ++j)
	{
          particle = (GenParticle*) tower->Particles.At(j);
          // Select tower information when generated particles are electrons
          if (particle->PID == 13)
          {	
 	    // Tower energy resolution divided into energy ranges
            if (particle->E > 0.0 && particle->E <= 5.0)
	    {
              if (particle->Eta >= -1.0 && particle->Eta <= 1.0)
	        hTowerERes_Eta_N10_P10_P_00_05->Fill((tower->E - particle->E)/particle->E);
            }
	    else if (particle->E > 5.0 && particle->E <= 10.0)
            {
              if (particle->Eta >= -1.0 && particle->Eta <= 1.0)
                hTowerERes_Eta_N10_P10_P_05_10->Fill((tower->E - particle->E)/particle->E);
            }
	    else if (particle->E > 10.0 && particle->E <= 15.0)
            {
              if (particle->Eta >= -1.0 && particle->Eta <= 1.0)
                hTowerERes_Eta_N10_P10_P_10_15->Fill((tower->E - particle->E)/particle->E);
            }
	    else if (particle->E > 15.0 && particle->E <= 20.0)
            {
              if (particle->Eta >= -1.0 && particle->Eta <= 1.0)
                hTowerERes_Eta_N10_P10_P_15_20->Fill((tower->E - particle->E)/particle->E);
            }
	  } // if statement for electron PID
          else
	    cout << "** Miss electron PID " << entry << " " << particle->PID << endl; 
        } // for statement for GenParticle loop
      } // for statement for Tower loop
    } // if statement for Tower branch
    else
      n_miss_branch++;
  } // for statement for all events

  // Sanity checks
  int total_nevt = hTowerERes_Eta_N10_P10_P_00_05->GetEntries() + hTowerERes_Eta_N10_P10_P_05_10->GetEntries() + 
		   hTowerERes_Eta_N10_P10_P_10_15->GetEntries() + hTowerERes_Eta_N10_P10_P_15_20->GetEntries();
  cout << "** Nevt proceesed: " << total_nevt << " out of " << allEntries << endl;
  cout << "** Missing tower branch nevt: " << n_miss_branch << " " << endl;
  cout << "** Total events through: " << total_nevt + n_miss_branch << endl;

  // Plot energy resolution each energy region
  TCanvas *cnv1 = new TCanvas("cnv1", "cnv1", 1000, 1000);
  cnv1->Divide(2,2);
  cnv1->cd(1);
  hTowerERes_Eta_N10_P10_P_00_05->GetXaxis()->CenterTitle(true);
  hTowerERes_Eta_N10_P10_P_00_05->GetYaxis()->CenterTitle(true);
  hTowerERes_Eta_N10_P10_P_00_05->GetXaxis()->SetNdivisions(505);
  hTowerERes_Eta_N10_P10_P_00_05->SetLineWidth(2);
  hTowerERes_Eta_N10_P10_P_00_05->Fit("gaus");
  TF1 *gaus_hTowerERes_Eta_N10_P10_P_00_05 = hTowerERes_Eta_N10_P10_P_00_05->GetFunction("gaus");
  gaus_hTowerERes_Eta_N10_P10_P_00_05->SetLineWidth(2);
  gaus_hTowerERes_Eta_N10_P10_P_00_05->SetLineColor(kRed);
  hTowerERes_Eta_N10_P10_P_00_05->Draw();
  TLegend* leg1 = new TLegend(0.17,0.7,0.45,0.85);
  leg1->SetBorderSize(0);
  leg1->SetTextColor(kRed+1);
  leg1->AddEntry((TObject*)0, "0<E^{gen}GeV<=5","");
  leg1->AddEntry((TObject*)0, Form("Mean = %.2f",gaus_hTowerERes_Eta_N10_P10_P_00_05->GetParameter(1)),"");
  leg1->AddEntry((TObject*)0, Form("Sigma = %.2f",gaus_hTowerERes_Eta_N10_P10_P_00_05->GetParameter(2)),"");
  leg1->Draw();
  cnv1->cd(2);
  hTowerERes_Eta_N10_P10_P_05_10->GetXaxis()->CenterTitle(true);
  hTowerERes_Eta_N10_P10_P_05_10->GetYaxis()->CenterTitle(true);
  hTowerERes_Eta_N10_P10_P_05_10->GetXaxis()->SetNdivisions(505);
  hTowerERes_Eta_N10_P10_P_05_10->SetLineWidth(2);
  hTowerERes_Eta_N10_P10_P_05_10->Fit("gaus");
  TF1 *gaus_hTowerERes_Eta_N10_P10_P_05_10 = hTowerERes_Eta_N10_P10_P_05_10->GetFunction("gaus");
  gaus_hTowerERes_Eta_N10_P10_P_05_10->SetLineWidth(2);
  gaus_hTowerERes_Eta_N10_P10_P_05_10->SetLineColor(kRed);
  hTowerERes_Eta_N10_P10_P_05_10->Draw();
  TLegend* leg2 = new TLegend(0.17,0.7,0.45,0.85);
  leg2->SetBorderSize(0);
  leg2->SetTextColor(kRed+1);
  leg2->AddEntry((TObject*)0, "5<E^{gen}GeV<=10","");
  leg2->AddEntry((TObject*)0, Form("Mean = %.2f",gaus_hTowerERes_Eta_N10_P10_P_05_10->GetParameter(1)),"");
  leg2->AddEntry((TObject*)0, Form("Sigma = %.2f",gaus_hTowerERes_Eta_N10_P10_P_05_10->GetParameter(2)),"");
  leg2->Draw();  
  cnv1->cd(3);
  hTowerERes_Eta_N10_P10_P_10_15->GetXaxis()->CenterTitle(true);
  hTowerERes_Eta_N10_P10_P_10_15->GetYaxis()->CenterTitle(true);
  hTowerERes_Eta_N10_P10_P_10_15->GetXaxis()->SetNdivisions(505);
  hTowerERes_Eta_N10_P10_P_10_15->SetLineWidth(2);
  hTowerERes_Eta_N10_P10_P_10_15->Fit("gaus");
  TF1 *gaus_hTowerERes_Eta_N10_P10_P_10_15 = hTowerERes_Eta_N10_P10_P_10_15->GetFunction("gaus");
  gaus_hTowerERes_Eta_N10_P10_P_10_15->SetLineWidth(2);
  gaus_hTowerERes_Eta_N10_P10_P_10_15->SetLineColor(kRed);
  hTowerERes_Eta_N10_P10_P_10_15->Draw();
  TLegend* leg3 = new TLegend(0.17,0.7,0.45,0.85);
  leg3->SetBorderSize(0);
  leg3->SetTextColor(kRed+1);
  leg3->AddEntry((TObject*)0, "10<E^{gen}GeV<=15","");
  leg3->AddEntry((TObject*)0, Form("Mean = %.2f",gaus_hTowerERes_Eta_N10_P10_P_10_15->GetParameter(1)),"");
  leg3->AddEntry((TObject*)0, Form("Sigma = %.2f",gaus_hTowerERes_Eta_N10_P10_P_10_15->GetParameter(2)),"");
  leg3->Draw();
  cnv1->cd(4);
  hTowerERes_Eta_N10_P10_P_15_20->GetXaxis()->CenterTitle(true);
  hTowerERes_Eta_N10_P10_P_15_20->GetYaxis()->CenterTitle(true);
  hTowerERes_Eta_N10_P10_P_15_20->GetXaxis()->SetNdivisions(505);
  hTowerERes_Eta_N10_P10_P_15_20->SetLineWidth(2);
  hTowerERes_Eta_N10_P10_P_15_20->Fit("gaus");
  TF1 *gaus_hTowerERes_Eta_N10_P10_P_15_20 = hTowerERes_Eta_N10_P10_P_15_20->GetFunction("gaus");
  gaus_hTowerERes_Eta_N10_P10_P_15_20->SetLineWidth(2);
  gaus_hTowerERes_Eta_N10_P10_P_15_20->SetLineColor(kRed);
  hTowerERes_Eta_N10_P10_P_15_20->Draw();
  TLegend* leg4 = new TLegend(0.17,0.7,0.45,0.85);
  leg4->SetBorderSize(0);
  leg4->SetTextColor(kRed+1);
  leg4->AddEntry((TObject*)0, "15<E^{gen}GeV<=20","");
  leg4->AddEntry((TObject*)0, Form("Mean = %.2f",gaus_hTowerERes_Eta_N10_P10_P_15_20->GetParameter(1)),"");
  leg4->AddEntry((TObject*)0, Form("Sigma = %.2f",gaus_hTowerERes_Eta_N10_P10_P_15_20->GetParameter(2)),"");
  leg4->Draw();
  cnv1->SaveAs("./plots/hElectronERes.png");
  cnv1->SaveAs("./plots/hElectronERes.pdf");
 
  // Energy resolution data 
  const int n = 4;
  Double_t x[n] = {2.5, 7.5, 12.5, 17.5};
  Double_t x_error[n] = {0., 0., 0., 0.};
  Double_t y_Eta_N10_P10[n] = {gaus_hTowerERes_Eta_N10_P10_P_00_05->GetParameter(2)*100., gaus_hTowerERes_Eta_N10_P10_P_05_10->GetParameter(2)*100.,
			       gaus_hTowerERes_Eta_N10_P10_P_10_15->GetParameter(2)*100., gaus_hTowerERes_Eta_N10_P10_P_15_20->GetParameter(2)*100.};
  Double_t y_Eta_N10_P10_error[n] = {gaus_hTowerERes_Eta_N10_P10_P_00_05->GetParError(2)*100., gaus_hTowerERes_Eta_N10_P10_P_05_10->GetParError(2)*100.,
				     gaus_hTowerERes_Eta_N10_P10_P_10_15->GetParError(2)*100., gaus_hTowerERes_Eta_N10_P10_P_15_20->GetParError(2)*100.};
  TGraphErrors* gEvsERes_Eta_N10_P10 = new TGraphErrors(n, x, y_Eta_N10_P10, x_error, y_Eta_N10_P10_error);
  // ePIC barrel ecal energy resoltuion
  TF1* fepic = new TF1("fepic", "sqrt(pow([0]/sqrt(x),2) + [1]*[1])", 0.05, 25.0);
  fepic->SetParName(0,"S");      // param names
  fepic->SetParName(1,"C");
  fepic->SetParameter(0,4.9); // initial param values
  fepic->SetParameter(1,0.1); // initial param values
  fepic->SetLineColor(kRed);
  fepic->SetLineStyle(2);
  fepic->SetLineWidth(2);  
  // Energy resolution fit function
  TF1* func = new TF1("func", "sqrt(pow([0]/sqrt(x),2) + [1]*[1])", 0.05, 25.0);
  // Plot energy resoltuion curve 
  TCanvas *cnv2 = new TCanvas("cnv2", "cnv2");
  gEvsERes_Eta_N10_P10->GetXaxis()->CenterTitle(true);
  gEvsERes_Eta_N10_P10->GetYaxis()->CenterTitle(true);
  gEvsERes_Eta_N10_P10->SetTitle("-1.0<=#eta<=1.0");
  gEvsERes_Eta_N10_P10->GetXaxis()->SetTitle("E [GeV]");
  gEvsERes_Eta_N10_P10->GetYaxis()->SetTitle("#Delta E/E [%]");
  gEvsERes_Eta_N10_P10->GetYaxis()->SetRangeUser(0.,6.);
  gEvsERes_Eta_N10_P10->SetLineColor(kBlack);
  gEvsERes_Eta_N10_P10->SetMarkerStyle(21);
  gEvsERes_Eta_N10_P10->SetMarkerColor(kBlack);
  gEvsERes_Eta_N10_P10->Fit(func,"R","", 0.1, 25.);
  TF1 * fit_gEvsERes_Eta_N10_P10 = gEvsERes_Eta_N10_P10->GetFunction("func");
  fit_gEvsERes_Eta_N10_P10->SetLineColor(kBlack);
  fit_gEvsERes_Eta_N10_P10->SetLineStyle(2);
  gEvsERes_Eta_N10_P10->Draw("AEP");
  fepic->Draw("SAME"); 
  TLegend* leg = new TLegend(0.2,0.7,0.8,0.85);
  leg->SetBorderSize(0);
  leg->AddEntry("fepic", "ePIC Barrel ECAL 4.9/#surdE #oplus 0.1","l");
  leg->AddEntry(fit_gEvsERes_Eta_N10_P10, Form("Delphes Barrel ECAL %.2f/#surdE #oplus %.2f",fit_gEvsERes_Eta_N10_P10->GetParameter(0),fit_gEvsERes_Eta_N10_P10->GetParameter(1)),"l");
  leg->Draw();
  cnv2->SaveAs("./plots/hElectronEResSummary.png");
  cnv2->SaveAs("./plots/hElectronEResSummary.pdf");

  cout << "** Done..." << endl;

  delete treeReader;
  delete chain;
}
