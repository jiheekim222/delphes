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

TH1F* GetHisto(TFile* f, const char* hname);
TH2F* GetHisto2D(TFile* f, const char* hname);
void SetStyle();
void SetHisto(TH1F* h, int c);
void SetGraph(TGraphErrors* gr,int c,const char* title="");
void SetPad1();
void SetPad2();
TF1* GetFitGaus(TH1F* h, const char* tag, int ieta, int ip,float min, float max);
TF1* GetFitMom(TGraphErrors* gr,const char* tag,int ieta,float min, float max);
TF1* GetFitE(TGraphErrors* gr,const char* tag,int ieta,float min, float max);
TF1* GetFitDCA(TGraphErrors* gr,const char* tag,int ieta,float min, float max);
TGraphErrors* GetGraphRes(TF1** fit,TH1F* hSpec,int ieta,int nBin,
                          const char* tag,float* array,float s=100);
double GetMean(TH1F* h, float* array, int iBin);
void EtaLabel(int ieta);
//------------------------------------------------------------------------------
void singlePiRes(const char* fname)
{
  int i,j;
  
  SetStyle();

  TFile* f=new TFile(fname);

  //-----------------
  // ePIC performance
  //-----------------
  TF1* fePICPRes[nEtaBin];
  for (i=0;i<nEtaBin;i++) {
    fePICPRes[i]=new TF1(Form("fePICPRes%d",i),"sqrt(([0]*x)^2+[1]^2)",0,10);
    fePICPRes[i]->SetLineStyle(7);
    fePICPRes[i]->SetLineColor(2);
  }
  fePICPRes[0]->SetParameter(0,0.048); fePICPRes[0]->SetParameter(1,0.5);
  fePICPRes[1]->SetParameter(0,0.053); fePICPRes[1]->SetParameter(1,0.4);
  fePICPRes[2]->SetParameter(0,0.053); fePICPRes[2]->SetParameter(1,0.4);
  fePICPRes[3]->SetParameter(0,0.051); fePICPRes[3]->SetParameter(1,0.5);

  TF1* fePICDCARes[nEtaBin];
  for (i=0;i<nEtaBin;i++) {
    fePICDCARes[i]=new TF1(Form("fePICDCARes%d",i),"sqrt(([0]/x)^2+[1]^2)",0,10);
    fePICDCARes[i]->SetLineStyle(7);
    fePICDCARes[i]->SetLineColor(2);
  }
  // convert mm to um
  fePICDCARes[0]->SetParameter(0,0.0279e3); fePICDCARes[0]->SetParameter(1,0.0067e3);
  fePICDCARes[1]->SetParameter(0,0.0247e3); fePICDCARes[1]->SetParameter(1,0.0066e3);
  fePICDCARes[2]->SetParameter(0,0.0241e3); fePICDCARes[2]->SetParameter(1,0.0065e3);
  fePICDCARes[3]->SetParameter(0,0.0283e3); fePICDCARes[3]->SetParameter(1,0.0068e3);

  TF1* fePICERes=new TF1("fePICERes","sqrt(([0]/sqrt(x))^2+[1]^2)",0,10);
  fePICERes->SetParameter(0,75);
  fePICERes->SetParameter(1,14.5);
  fePICERes->SetLineStyle(7);
  fePICERes->SetLineColor(2);
  
  //----------------- 
  // declare fits
  //----------------- 
  TF1* fitEResDis[nEtaBin][nEBin];
  TF1* fitPResDis[nEtaBin][nPBin];
  TF1* fitPtResDis[nEtaBin][nPtBin];
  TF1* fitDCAtResDis[nEtaBin][nPtBin];
  TF1* fitDCAzResDis[nEtaBin][nPtBin];

  TF1* fitERes[nEtaBin];
  TF1* fitPRes[nEtaBin];
  TF1* fitPtRes[nEtaBin];
  TF1* fitDCAtRes[nEtaBin];
  TF1* fitDCAzRes[nEtaBin];

  //-----------------
  // declare graphs
  //-----------------
  TGraphErrors* grERes[nEtaBin];
  TGraphErrors* grPRes[nEtaBin];
  TGraphErrors* grPtRes[nEtaBin];
  TGraphErrors* grDCAtRes[nEtaBin];
  TGraphErrors* grDCAzRes[nEtaBin];

  //-----------------
  // get hsito
  //-----------------
  //particle branch
  TH1F *fParticleEta=GetHisto(f,"hParticleEta");
  TH1F *fParticleE[nEtaBin];
  TH1F *fParticleP[nEtaBin];
  TH1F *fParticlePt[nEtaBin];

  for (i=0;i<nEtaBin;i++) {
    fParticleE[i]=GetHisto(f,Form("hParticleE_eta%d",i));
    fParticleP[i]=GetHisto(f,Form("hParticlePt_eta%d",i));
    fParticlePt[i]=GetHisto(f,Form("hParticleP_eta%d",i));
  }

  //track branch
  TH1F *fTrkEta=GetHisto(f,"hEtaTrk");
  TH1F *fTrkP=GetHisto(f,"hPTrk");
  TH1F *fTrkPt=GetHisto(f,"hPTrk");
  TH1F *fTrkOuterPhi=GetHisto(f,"hOuterPhiTrk");
  TH1F *fTrkOuterEta=GetHisto(f,"hOuterEtaTrk");

  TH1F *fTrkDeltaEta[nEtaBin];
  TH1F *fTrkDeltaP[nEtaBin][nPBin];
  TH1F *fTrkDeltaPt[nEtaBin][nPtBin];
  TH1F *fTrkDeltaDCAt[nEtaBin][nPtBin];
  TH1F *fTrkDeltaDCAz[nEtaBin][nPtBin];

  for (i=0;i<nEtaBin;i++) {
    fTrkDeltaEta[i]=GetHisto(f,Form("Trk_delta_eta_eta%d",i));
    for (j=0;j<nPBin;j++) {
      fTrkDeltaP[i][j]=GetHisto(f,Form("Trk_delta_p_eta%d_p%d",i,j));
      SetHisto(fTrkDeltaP[i][j], 1);
      fitPResDis[i][j]=GetFitGaus(fTrkDeltaP[i][j],"P",i,j,-0.03,0.03);
    }
    grPRes[i]=GetGraphRes(fitPResDis[i],fParticleP[i],i,nPBin,"p",p);
    SetGraph(grPRes[i],1,";p^{gen} (GeV);#sigma((p^{gen}-p^{rec})/p^{gen}) (%)");
    grPRes[i]->GetXaxis()->SetLimits(0,p[nPBin]);
    grPRes[i]->SetMinimum(0.15);
    grPRes[i]->SetMaximum(0.8);
    fitPRes[i]=GetFitMom(grPRes[i],"p",i,0,p[nPBin]);

    for (j=0;j<nPtBin;j++) {
      fTrkDeltaPt[i][j]=GetHisto(f,Form("Trk_delta_pt_eta%d_pt%d",i,j));
      SetHisto(fTrkDeltaPt[i][j],1);
      fitPtResDis[i][j]=GetFitGaus(fTrkDeltaPt[i][j],"P",i,j,-0.03,0.03);

      fTrkDeltaDCAt[i][j]=GetHisto(f,Form("Trk_delta_DCAt_eta%d_pt%d",i,j));
      SetHisto(fTrkDeltaDCAt[i][j],1);
      fTrkDeltaDCAt[i][j]->SetTitle(";DCA_{T} (mm);");
      fitDCAtResDis[i][j]=GetFitGaus(fTrkDeltaDCAt[i][j],"P",i,j,-0.03,0.03);

      fTrkDeltaDCAz[i][j]=GetHisto(f,Form("Trk_delta_DCAz_eta%d_pt%d",i,j));
      SetHisto(fTrkDeltaDCAz[i][j],1);
      fTrkDeltaDCAz[i][j]->SetTitle(";DCA_{z} (mm);");
      fitDCAzResDis[i][j]=GetFitGaus(fTrkDeltaDCAz[i][j],"P",i,j,-0.03,0.03);
    }
    
    grPtRes[i]=GetGraphRes(fitPtResDis[i],fParticlePt[i],i,nPtBin,"pt",pt);
    SetGraph(grPtRes[i],1,";p_{T}^{gen} (GeV);#sigma((p_{T}^{gen}-p_{T}^{rec})/p_{T}^{gen}) (%)");
    grPtRes[i]->GetXaxis()->SetLimits(0,pt[nPBin]);
    grPtRes[i]->SetMinimum(0.46);
    grPtRes[i]->SetMaximum(0.75);
    fitPtRes[i]=GetFitMom(grPtRes[i],"pt",i,0,pt[nPBin]);

    grDCAtRes[i]=GetGraphRes(fitDCAtResDis[i],fParticlePt[i],i,nPtBin,"pt",pt,1e3);
    SetGraph(grDCAtRes[i],1,";p_{T}^{gen} (GeV);DCA_{T} (#mum)");
    grDCAtRes[i]->GetXaxis()->SetLimits(0,pt[nPBin]);
    grDCAtRes[i]->SetMinimum(5);
    grDCAtRes[i]->SetMaximum(100);
    fitDCAtRes[i]=GetFitDCA(grDCAtRes[i],"pt",i,0,pt[nPBin]);

    grDCAzRes[i]=GetGraphRes(fitDCAzResDis[i],fParticlePt[i],i,nPtBin,"pt",pt,1e3);
    SetGraph(grDCAzRes[i],1,";p_{T}^{gen} (GeV);DCA_{z} (#mum)");
    grDCAzRes[i]->GetXaxis()->SetLimits(0,pt[nPtBin]);
    grDCAzRes[i]->SetMinimum(0);
    grDCAzRes[i]->SetMaximum(30);
    fitDCAzRes[i]=GetFitDCA(grDCAzRes[i],"pt",i,0,pt[nPtBin]);
  }

  //tower branch
  TH1F *fTowEta=GetHisto(f,"hEtaTow");
  TH1F *fTowPhi=GetHisto(f,"hPhiTow");
  TH1F *fTowE=GetHisto(f,"hETow");
  TH1F *fTowEhad=GetHisto(f,"hEhadTow");
  TH1F *fTowEem=GetHisto(f,"hEemTow");
  TH2F *fTowEfracHad=GetHisto2D(f,"hEfracHadTow");
  TH2F *fTowEfracEM=GetHisto2D(f,"hEfracEMTow");

  TH1F *fTowDeltaEta[nEtaBin];
  TH1F *fTowDeltaE[nEtaBin][nEBin];

  for (i=0;i<nEtaBin;i++) {
    fTowDeltaEta[i]=GetHisto(f,Form("Tow_delta_eta_eta%d",i));
    for (j=0;j<nEBin;j++) {
      fTowDeltaE[i][j]=GetHisto(f,Form("Tow_delta_E_eta%d_E%d",i,j));
      SetHisto(fTowDeltaE[i][j],1);
      fitEResDis[i][j]=GetFitGaus(fTowDeltaE[i][j],"E",i,j,-0.3,0.3);
    }
    grERes[i]=GetGraphRes(fitEResDis[i],fParticleE[i],i,nEBin,"E",E);
    SetGraph(grERes[i],1,";E^{gen} (GeV);#sigma((E^{gen}-E^{rec})/E^{gen}) (%)");
    grERes[i]->GetXaxis()->SetLimits(0,E[nEBin]);
    grERes[i]->SetMinimum(18);
    grERes[i]->SetMaximum(70);
    fitERes[i]=GetFitE(grERes[i],"E",i,0,E[nEBin]);
  }

  f->Close();

  //-----------------
  // set legend
  //-----------------
  TLegend* legP[nEtaBin];
  for (i=0;i<nEtaBin;i++) {
    legP[i]=new TLegend(0.4,0.15,0.98,0.47);
    legP[i]->SetTextFont(43);
    legP[i]->SetTextSize(20);
    legP[i]->SetFillStyle(0);
    legP[i]->AddEntry(fitPRes[i],"Delphes","l");
    legP[i]->AddEntry((TObject*) 0,
		      Form("%.3fp#oplus%.2f",
			   fitPRes[i]->GetParameter(0),
			   fitPRes[i]->GetParameter(1)),
		      "");
    legP[i]->AddEntry(fePICPRes[i],"ePIC","l");
    legP[i]->AddEntry((TObject*) 0,
		      Form("%.3fp#oplus%.2f",
			   fePICPRes[i]->GetParameter(0),
			   fePICPRes[i]->GetParameter(1)),
		      "");
			
  }

  TLegend* legDCA[nEtaBin];
  for (i=0;i<nEtaBin;i++) {
    legDCA[i]=new TLegend(0.4,0.55,0.98,0.9);
    legDCA[i]->SetTextFont(43);
    legDCA[i]->SetTextSize(20);
    legDCA[i]->SetFillStyle(0);
    legDCA[i]->AddEntry(fitDCAtRes[i],"Delphes","l");
    legDCA[i]->AddEntry((TObject*) 0,
			Form("%.3f/p#oplus%.2f",
			     fitDCAtRes[i]->GetParameter(0),
			     fitDCAtRes[i]->GetParameter(1)),
			"");
    legDCA[i]->AddEntry(fePICDCARes[i],"ePIC","l");
    legDCA[i]->AddEntry((TObject*) 0,
			Form("%.3f/p#oplus%.2f",
			     fePICDCARes[i]->GetParameter(0),
			     fePICDCARes[i]->GetParameter(1)),
			"");
  }
  
  TLegend* legE[nEtaBin];
  for (i=0;i<nEtaBin;i++) {
    legE[i]=new TLegend(0.4,0.55,0.98,0.9);
    legE[i]->SetTextFont(43);
    legE[i]->SetTextSize(20);
    legE[i]->SetFillStyle(0);
    legE[i]->AddEntry(fitPRes[i],"Delphes","l");
    legE[i]->AddEntry((TObject*) 0,
		      Form("%.1f/#sqrt{E}#oplus%.1f",
			   fitERes[i]->GetParameter(0),
			   fitERes[i]->GetParameter(1)),
		      "");
    legE[i]->AddEntry(fePICERes,"ePIC","l");
    legE[i]->AddEntry((TObject*) 0,
		      Form("%.1f/#sqrt{E}#oplus%.1f",
			   fePICERes->GetParameter(0),
			   fePICERes->GetParameter(1)),
		      "");

  }

  //----------------- 
  // draw canvas
  //----------------- 
  //TCanvas* cPResDis[nEtaBin];
  //TCanvas* cPtResDis[nEtaBin];
  TCanvas* cDCAtResDis[nEtaBin];
  //TCanvas* cDCAzResDis[nEtaBin];
  //TCanvas* cEResDis[nEtaBin];

 
  for (i=0;i<nEtaBin;i++) {
    /*
    cPResDis[i]=new TCanvas(Form("cPResDis_eta%d",i),
			    Form("P res dist| %.1f<#eta<%.1f",eta[i],eta[i+1]),
			    1200,800);
    cPResDis[i]->Range(0,0,1,1);
    cPResDis[i]->Divide(ceil(nPBin/2.),2);

    cPtResDis[i]=new TCanvas(Form("cPtResDis_eta%d",i),
			     Form("Pt res dist| %.1f<#eta<%.1f",eta[i],eta[i+1]),
			     1200,800);
    cPtResDis[i]->Range(0,0,1,1);
    cPtResDis[i]->Divide(ceil(nPtBin/2.),2);
    */
    
    cDCAtResDis[i]=new TCanvas(Form("cDCAtResDis_eta%d",i),
			       Form("DCAt res dist| %.1f<#eta<%.1f",eta[i],eta[i+1]),
			       1200,1000);
    cDCAtResDis[i]->Range(0,0,1,1);
    cDCAtResDis[i]->Divide(ceil(nPtBin/3.),3);

    /*
    cDCAzResDis[i]=new TCanvas(Form("cDCAzResDis_eta%d",i),
                               Form("DCAz res dist| %.1f<#eta<%.1f",eta[i],eta[i+1]),
                               1200,800);
    cDCAzResDis[i]->Range(0,0,1,1);
    cDCAzResDis[i]->Divide(ceil(nPtBin/2.),2);

    cEResDis[i]=new TCanvas(Form("cEResDis_eta%d",i),
			    Form("E res dist| %.1f<#eta<%.1f",eta[i],eta[i+1]),
			    1200,800);
    cEResDis[i]->Range(0,0,1,1);
    cEResDis[i]->Divide(ceil(nEBin/2.),2);
    
    for (j=0;j<nPBin;j++) {
      cPResDis[i]->cd(j+1);
      SetPad1();
      fTrkDeltaP[i][j]->Draw();
      if (!fitPResDis[i][j]) continue;
      fitPResDis[i][j]->Draw("same");
      fTrkDeltaP[i][j]->Draw("same");
    }

    for (j=0;j<nPtBin;j++) {
      cPtResDis[i]->cd(j+1);
      SetPad1();
      fTrkDeltaPt[i][j]->Draw();
      if (!fitPtResDis[i][j]) continue;
      fitPtResDis[i][j]->Draw("same");
      fTrkDeltaPt[i][j]->Draw("same");
    }
    */

    for (j=0;j<nPtBin;j++) {
      cDCAtResDis[i]->cd(j+1);
      SetPad1();
      fTrkDeltaDCAt[i][j]->Draw();
      if (!fitDCAtResDis[i][j]) continue;
      fitDCAtResDis[i][j]->Draw("same"); 
      fTrkDeltaDCAt[i][j]->Draw("same"); 
    } 
    
    /*
    for (j=0;j<nPtBin;j++) {
      cDCAzResDis[i]->cd(j+1);
      SetPad1();
      fTrkDeltaDCAz[i][j]->Draw();
      if (!fitDCAzResDis[i][j]) continue;
      fitDCAzResDis[i][j]->Draw("same");
      fTrkDeltaDCAz[i][j]->Draw("same");
    }
    
    for (j=0;j<nEBin;j++) {
      cEResDis[i]->cd(j+1);
      SetPad1();
      fTowDeltaE[i][j]->Draw();
      if (! fitEResDis[i][j]) continue;
      fitEResDis[i][j]->Draw("same");
      fTowDeltaE[i][j]->Draw("same");
      cout<<"mean = "<<fitEResDis[i][j]->GetParameter(1)<<endl;
    }
    */
    
    //cPResDis[i]->SaveAs(Form("pResDis_eta%d.png",i));
    //cPtResDis[i]->SaveAs(Form("ptResDis_eta%d.png",i));
    cDCAtResDis[i]->SaveAs(Form("DCAtResDis_eta%d.png",i));
    //cDCAzResDis[i]->SaveAs(Form("DCAzResDis_eta%d.png",i));
    //cEResDis[i]->SaveAs(Form("EResDis_eta%d.png",i));
  }

  //TCanvas* cPRes=new TCanvas("cPRes","P res",1200,350);
  //TCanvas* cPtRes=new TCanvas("cPtRes","Pt res",1200,350);
  TCanvas* cDCAtRes=new TCanvas("cDCAtRes","DCAt res",1200,350);
  //TCanvas* cERes=new TCanvas("cERes","E res",1200,350);
  
  
  //cPRes->Range(0,0,1,1);
  //cPtRes->Range(0,0,1,1);
  cDCAtRes->Range(0,0,1,1);
  //cERes->Range(0,0,1,1);

  //cPRes->Divide(nEtaBin);
  //cPtRes->Divide(nEtaBin);
  cDCAtRes->Divide(nEtaBin);
  //cERes->Divide(nEtaBin);

  for (i=0;i<nEtaBin;i++) {
    /*
    cPRes->cd(i+1);
    SetPad2();
    grPRes[i]->Draw("ap");
    fePICPRes[i]->Draw("same");
    fitPRes[i]->Draw("same");
    grPRes[i]->Draw("psame");
    EtaLabel(i);
    legP[i]->Draw();
    
    cPtRes->cd(i+1);
    SetPad2();
    grPtRes[i]->Draw("ap");
    fitPtRes[i]->Draw("same");
    grPtRes[i]->Draw("psame");
    EtaLabel(i);
    */

    cDCAtRes->cd(i+1);
    SetPad2();
    grDCAtRes[i]->Draw("ap");
    fitDCAtRes[i]->Draw("same");
    fePICDCARes[i]->Draw("same");
    grDCAtRes[i]->Draw("psame");
    EtaLabel(i);
    legDCA[i]->Draw();
    
    /*
    cERes->cd(i+1);
    SetPad2();
    grERes[i]->Draw("ap");
    fitERes[i]->Draw("same");
    fePICERes->Draw("same");
    //fitERes[i]->Draw("same");
    grERes[i]->Draw("PSAME");
    EtaLabel(i);
    legE[i]->Draw();
    */
  }

  //cPRes->SaveAs("singlePi_pRes.png");
  //cPtRes->SaveAs("singlePi_ptRes.png");
  cDCAtRes->SaveAs("singlePi_DCAtRes.png");
  //cERes->SaveAs("singlePi_ERes.png");
}
//------------------------------------------------------------------------------
TH1F* GetHisto(TFile* f, const char* hname)
{
  TH1F* tmp=(TH1F*) f->Get(hname);
  if (!tmp) {
    cout<<"can't find "<<hname<<endl;
    return 0;
  }

  gROOT->cd();
  TH1F* h=(TH1F*) tmp->Clone(Form("%s_clone",hname));

  delete tmp;
  return h;
}
//------------------------------------------------------------------------------
TH2F* GetHisto2D(TFile* f, const char* hname)
{
  TH2F* tmp=(TH2F*) f->Get(hname);
  gROOT->cd();
  TH2F* h=(TH2F*) tmp->Clone(Form("%s_clone",hname));

  delete tmp;
  return h;
}
//------------------------------------------------------------------------------
void SetStyle()
{
  gStyle->SetPalette(1);
  gStyle->SetLegendBorderSize(0);
  TGaxis().SetMaxDigits(4);
}
//------------------------------------------------------------------------------
TF1* GetFitGaus(TH1F* h, const char* tag, int ieta, int ip,float min, float max)
{
  if (h->GetEntries()==0) return 0;


  TF1* fit=new TF1(Form("fit%s%d%d",tag,ieta,ip),
		   "[0]*exp(-((x-[1])^2/(2*[2]^2)))",
		   h->GetBinCenter(1),
		   h->GetBinCenter(h->GetNbinsX()));
  fit->SetParameter(0,h->GetMaximum());
  fit->SetParameter(1,0);
  fit->SetParameter(2,h->GetStdDev());
  fit->SetLineStyle(7);
  fit->SetLineColor(2);

  h->Fit(fit,"NR","",min,max);
  
  return fit;
}
//------------------------------------------------------------------------------
TF1* GetFitMom(TGraphErrors* gr,const char* tag,int ieta,float min, float max)
{
  TF1* fit=new TF1(Form("fit%s%d",tag,ieta),"sqrt(([0]*x)^2+[1]^2)",min,max);
  fit->SetParameter(0,0.48);
  fit->SetParameter(1,0.5);
  fit->SetLineStyle(7);
  fit->SetLineColor(1);

  gr->Fit(fit,"N");//,"",min,max);

  return fit;

}
//------------------------------------------------------------------------------  
TF1* GetFitE(TGraphErrors* gr,const char* tag,int ieta,float min, float max)
{
  TF1* fit=new TF1(Form("fit%s%d",tag,ieta),"sqrt(([0]/sqrt(x))^2+[1]^2)",min,max);
  fit->SetParameter(0,75);
  fit->SetParameter(1,14);
  fit->SetLineStyle(7);
  fit->SetLineColor(1);

  gr->Fit(fit,"NR","",min,max);

  return fit;

}
//------------------------------------------------------------------------------ 
TF1* GetFitDCA(TGraphErrors* gr,const char* tag,int ieta,float min, float max)
{
  TF1* fit=new TF1(Form("fit%s%d",tag,ieta),"sqrt(([0]/x)^2+[1]^2)",min,max);
  fit->SetParameter(0,28);
  fit->SetParameter(1,7);
  fit->SetLineStyle(7);
  fit->SetLineColor(1);

  gr->Fit(fit,"N");//,"",min,max);                                                                                                                                              

  return fit;

}
//------------------------------------------------------------------------------
void SetHisto(TH1F* h,int c)
{
  h->SetStats(0);
  h->SetLineColor(c);
  h->SetMarkerColor(c);
  h->SetMarkerStyle(23+c);

  h->GetXaxis()->SetNdivisions(510);
  h->GetXaxis()->SetTitleFont(43);
  h->GetXaxis()->SetTitleSize(24);
  h->GetXaxis()->SetTitleOffset(0.9);
  h->GetXaxis()->SetLabelFont(43);
  h->GetXaxis()->SetLabelSize(20);

  h->GetYaxis()->SetNdivisions(510);
  h->GetYaxis()->SetTitleFont(43);
  h->GetYaxis()->SetTitleSize(24);
  h->GetYaxis()->SetLabelFont(43);
  h->GetYaxis()->SetLabelSize(20);
}
//------------------------------------------------------------------------------
void SetGraph(TGraphErrors* gr,int c,const char* title="")
{
  gr->SetTitle(title);
  gr->SetLineColor(c);
  gr->SetMarkerColor(c);
  gr->SetMarkerStyle(23+c);

  gr->GetXaxis()->SetNdivisions(510);
  gr->GetXaxis()->SetTitleFont(43);
  gr->GetXaxis()->SetTitleSize(24);
  gr->GetXaxis()->SetTitleOffset(0.9);
  gr->GetXaxis()->SetLabelFont(43);
  gr->GetXaxis()->SetLabelSize(20);

  gr->GetYaxis()->SetNdivisions(510);
  gr->GetYaxis()->SetTitleFont(43);
  gr->GetYaxis()->SetTitleSize(24);
  gr->GetYaxis()->SetLabelFont(43);
  gr->GetYaxis()->SetLabelSize(20);
}
//------------------------------------------------------------------------------
void SetPad1()
{
  gPad->SetTopMargin(0.1);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.005);
  gPad->SetLeftMargin(0.18);
}
//------------------------------------------------------------------------------ 
void SetPad2()
{
  gPad->SetTopMargin(0.001);
  gPad->SetBottomMargin(0.15);
  gPad->SetRightMargin(0.005);
  gPad->SetLeftMargin(0.32);
}
//------------------------------------------------------------------------------
TGraphErrors* GetGraphRes(TF1** fit,TH1F* hSpec,int ieta,int nBin,
			  const char* tag,float* array,float s=100)  //s: scale
{
  int i,cnt;
  double mean,res,err;

  cnt=0;
  for (i=0;i<nBin;i++) {
    if (fit[i]!=0 && GetMean(hSpec,array,i)>=0) cnt++;
  }

  TGraphErrors* gr=new TGraphErrors(cnt);
  gr->SetName(Form("gr%sRes%d",tag,ieta));

  cnt=0;
  for (i=0;i<nBin;i++) {
    if (!fit[i]) continue;
    
    mean=GetMean(hSpec,array,i);
    if (mean<0) continue;

    res=fabs(fit[i]->GetParameter(2))*s;
    err=fabs(fit[i]->GetParError(2))*s;
    //cout<<"mean="<<mean<<" , res="<<res<<endl;
    gr->SetPoint(cnt,mean, res);
    gr->SetPointError(cnt,0,err);
    cnt++;
  }

  return gr;
}
//------------------------------------------------------------------------------
double GetMean(TH1F* h, float* array, int iBin)
{
  int i;
  double mean=0;
  double cnt=0;
  
  for (i=0;i<h->GetNbinsX();i++) {
    if (h->GetBinCenter(i)>=array[iBin] && h->GetBinCenter(i)<array[iBin+1]) {
      mean=mean+h->GetBinContent(i)*h->GetBinCenter(i);
      cnt=cnt+h->GetBinContent(i);
    }
  }
  
  if (cnt) mean=mean/cnt;
  else mean=-1;

  return mean;
}
//------------------------------------------------------------------------------
void EtaLabel(int ieta)
{
  float xMin, xRange, yMax, yRange;

  gPad->Update();
  xMin=gPad->GetUxmin();
  xRange=gPad->GetUxmax()-xMin;
  yMax=gPad->GetUymax();
  yRange=yMax-gPad->GetUymin();

  TLatex* txt=new TLatex(xMin+0.2*xRange,yMax-0.08*yRange,
			 Form("%.1f #leq #eta < %.1f",eta[ieta],eta[ieta+1]));

  txt->SetTextFont(43);
  txt->SetTextSize(20);
  txt->Draw();
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------

