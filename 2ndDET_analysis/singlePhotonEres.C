//------------------------------------------------------------------------------  
// constant and binning

const int nEtaBin=4;
float eta[nEtaBin+1]={-1,-0.5,0,0.5,1};

const int nEBin=5;
float E[nEBin+1]={5,6,7,8,9,10};

//------------------------------------------------------------------------------ 

TH1F* GetHisto(TFile* f, const char* hname);
int GetEtaBin(float myEta);
int GetEBin(float myE);
void SetMax(TH1F** h, int n);
TF1* GetFit(TH1F* h,int ieta,int iE);
TF1* GetFit(TGraphErrors* gr, int ieta);
void SetHisto(TH1F* h, int c);
void SetGraph(TGraphErrors* gr, int c,const char* title);
void SetStyle();
void SetStyle(TStyle* kStyle);
double GetMean(TH1F* h, int iE);
TGraphErrors* GetGraph(TH1F* hE, TF1** fit, int ieta, int nBin);
void EtaLabel(int ieta);
//------------------------------------------------------------------------------
void singlePhotonEres(const char* fname)
{
  int i,j;
  
  // Setting for figures
  TStyle* kStyle = new TStyle("kStyle","Kim's Style");
  //SetStyle(kStyle);
  SetStyle();

  TFile* f=new TFile(fname);

  //------------------- 
  // book histo
  //------------------- 
  TH1F* hE[nEtaBin];
  TH1F* hResDis[nEtaBin][nEBin];
  TF1* fitRes[nEtaBin][nEBin];
  TF1* fitResFunc[nEtaBin];

  TGraphErrors* grRes[nEtaBin];

  //------------------- 
  // ePIC numbers
  //------------------- 
  TF1* fePIC=new TF1("fePIC","sqrt((4.9/sqrt(x))^2+0.1^2)",5,10);
  fePIC->SetLineColor(2);
  fePIC->SetLineStyle(7);

  //-------------------
  // get histo
  //-------------------
  for (i=0;i<nEtaBin;i++) {
    hE[i]=GetHisto(f,Form("hParticleE_eta%d",i));
    SetHisto(hE[i],i+1);
    hE[i]->SetTitle(";E_{gen} (GeV);");

    for (j=0;j<nEBin;j++) {
      hResDis[i][j]=GetHisto(f,Form("photon_delta_energy_eta%dE%d",i,j));
      SetHisto(hResDis[i][j],1);
      hResDis[i][j]->GetXaxis()->SetRangeUser(-0.1,0.1);
      fitRes[i][j]=GetFit(hResDis[i][j],i,j);
    }

    grRes[i]=GetGraph(hE[i],fitRes[i],i,nEBin);
    SetGraph(grRes[i],1,";E_{gen} (GeV);#sigma(#DeltaE/E_{gen}) (%)");
    //SetGraph(grRes[i],1,";E_{gen} (GeV);#sigma(#DeltaE/E_{gen}) (%%)");
    
    fitResFunc[i]=GetFit(grRes[i],i);
  }
  f->Close();
  
  SetMax(hE,nEtaBin);
  
  TLegend* legE=new TLegend(0.1,0.8,0.9,0.95);
  legE->SetTextFont(43);
  legE->SetTextSize(20);
  legE->SetFillStyle(0);
  legE->SetNColumns(4);
  for (i=0;i<nEtaBin;i++) legE->AddEntry(hE[i],Form("%.1f #leq #eta < %.1f",eta[i],eta[i+1]),"p");
    
  TLegend* legFit[i];
  for (i=0;i<nEtaBin;i++) {
    legFit[i]=new TLegend(0.4,0.6,0.98,0.9);
    legFit[i]->SetTextFont(43);
    legFit[i]->SetTextSize(20);
    legFit[i]->SetFillStyle(0);
    //legE->SetNColumns(4);
    legFit[i]->AddEntry(fitResFunc[i],"Delphes","l");
    legFit[i]->AddEntry((TObject*) 0,
                        Form("%.2f/#sqrt{E}#oplus%.2f",
                             fitResFunc[i]->GetParameter(0),
                             fitResFunc[i]->GetParameter(1)),
                        "");
    legFit[i]->AddEntry(fePIC,"ePIC","l");
    legFit[i]->AddEntry((TObject*) 0,"4.9/#sqrt{E}#oplus0.1","");
  }

  //-------------------
  // draw canvas
  //-------------------
  TCanvas* cEdis=new TCanvas("c energy","energy distribution",800,600);
  cEdis->Range(0,0,1,1);
  cEdis->cd();
  gPad->SetTopMargin(0.05);
  gPad->SetBottomMargin(0.15);
  hE[0]->Draw();
  for (i=0;i<nEtaBin;i++) {
    hE[i]->Draw("same");
  }
  legE->Draw();
  cEdis->SaveAs("singlePhotonEdis.png");
 
  TCanvas* cRes[nEtaBin];
  for (i=0;i<nEtaBin;i++) {
    cRes[i]=new TCanvas(Form("cRes%d",i),Form("%.1f < eta < %.1f",eta[i],eta[i+1]),1000,750);
    cRes[i]->Range(0,0,1,1);
    cRes[i]->Divide(ceil(nEBin/2.),2);
    
    for (j=0;j<nEBin;j++) {
      cRes[i]->cd(j+1);
      gPad->SetBottomMargin(0.15);
      hResDis[i][j]->Draw();
      fitRes[i][j]->Draw("same");
      hResDis[i][j]->Draw("same");
    }
    cRes[i]->SaveAs(Form("singlePhotonResDisEta%d.png",i));
  }

  TCanvas* cResE=new TCanvas("cResFunc","resolution",1200,400);
  cResE->Range(0,0,1,1);
  cResE->Divide(nEtaBin);
  for (i=0;i<nEtaBin;i++) {
    cResE->cd(i+1);
    gPad->SetTopMargin(0.01);
    gPad->SetBottomMargin(0.15);
    gPad->SetRightMargin(0.005);
    gPad->SetLeftMargin(0.25);
    grRes[i]->Draw("pa");
    fePIC->Draw("same");
    fitResFunc[i]->Draw("same");
    legFit[i]->Draw();
    EtaLabel(i);
  }
  cResE->SaveAs("singlePhotonResEfunc.png");
}

//------------------------------------------------------------------------------

TH1F* GetHisto(TFile* f, const char* hname)
{
  TH1F* tmp=(TH1F*) f->Get(hname);

  gROOT->cd();
  TH1F* h=(TH1F*) tmp->Clone(Form("%sClone",hname));
  
  delete tmp;
  return h;
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

void SetMax(TH1F** h, int n)
{
  int i;
  double max=-999;

  for (i=0;i<n;i++) {
    if (h[i]->GetMaximum()>max) max=h[i]->GetMaximum();
  }
  
  for (i=0;i<n;i++) {
    h[i]->SetMaximum(max*1.2);
  }
}
//------------------------------------------------------------------------------ 
TF1* GetFit(TH1F* h,int ieta,int iE)
{
  TF1* f=new TF1(Form("fitRes_eta%d_E%d",ieta,iE),"[0]*exp(-((x-[1])^2)/(2*[2]^2))",-0.05,0.05);
  f->SetParameter(0,h->GetMaximum());
  f->SetParameter(1,0);
  f->SetParameter(2,h->GetStdDev());
  
  f->SetLineStyle(7);
  
  h->Fit(f,"N");
  
  return f;
}
//------------------------------------------------------------------------------  
TF1* GetFit(TGraphErrors* gr, int ieta)
{
  TF1* f=new TF1(Form("fitResE_eta%d",ieta),"sqrt(([0]/sqrt(x))^2+[1]^2)",5,10);
  f->SetParameter(1,35e-3);
  f->SetParameter(0,2e-3);

  f->SetLineStyle(7);
  f->SetLineColor(1);

  gr->Fit(f,"N");

  return f;
}
//------------------------------------------------------------------------------
void SetHisto(TH1F* h, int c)
{
  h->SetStats(0);
  h->SetLineColor(c);
  h->SetMarkerColor(c);
  h->SetMarkerStyle(23+c);
  
  h->SetTitleFont(43);
  h->SetTitleSize(24);
  
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
//------------------------------------------------------------------------------
void SetStyle()
{
  gStyle->SetLegendBorderSize(0);
  TGaxis().SetMaxDigits(3);
}

//------------------------------------------------------------------------------
TGraphErrors* GetGraph(TH1F* hE, TF1** fit, int ieta, int nBin)
{
  int i;
  double meanE, res, err;

  TGraphErrors* gr=new TGraphErrors(nBin);
  gr->SetName(Form("grRes_eta%d",ieta));

  for (i=0;i<nBin;i++) {
    meanE=GetMean(hE,i);
    cout<<fit[i]->GetName()<<endl;
    res=fabs(fit[i]->GetParameter(2))*100;
    err=fabs(fit[i]->GetParError(2))*100;
    cout<<"E="<<meanE<<" , res="<<res<<endl;
    gr->SetPoint(i,meanE, res);
    gr->SetPointError(i,0,err);
  }

  return gr;
}
//------------------------------------------------------------------------------ 
double GetMean(TH1F* h, int iE)
{
  int i;
  double mean=0;
  double cnt=0;
  
  for (i=0;i<h->GetNbinsX();i++) {
    if (h->GetBinCenter(i)>=E[iE] && h->GetBinCenter(i)<E[iE+1]) {
      mean=mean+h->GetBinContent(i)*h->GetBinCenter(i);
      cnt=cnt+h->GetBinContent(i);
    }
  }
  mean=mean/cnt;
  return mean;
}
//------------------------------------------------------------------------------  
void SetGraph(TGraphErrors* gr, int c,const char* title)
{
  gr->SetTitle(title);
  gr->SetLineColor(c);
  gr->SetMarkerColor(c);
  gr->SetMarkerStyle(24+c);
  gr->SetMinimum(1.5);
  gr->SetMaximum(2.3);

  gr->GetXaxis()->SetNdivisions(505);
  gr->GetXaxis()->SetTitleFont(43);
  gr->GetXaxis()->SetTitleSize(24);
  gr->GetXaxis()->SetLabelFont(43);
  gr->GetXaxis()->SetLabelSize(20);

  gr->GetYaxis()->SetNdivisions(505);
  gr->GetYaxis()->SetTitleFont(43);
  gr->GetYaxis()->SetTitleSize(24);
  gr->GetYaxis()->SetLabelFont(43);
  gr->GetYaxis()->SetLabelSize(20);
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
  
  TLatex* txt=new TLatex(xMin+0.2*xRange,yMax-0.06*yRange,
			 Form("%.1f #leq #eta < %1.f",eta[ieta],eta[ieta+1]));

  //cout<<xMin+0.2*xRange<<", "<<yMax-0.1*yRange<<endl;

  txt->SetTextFont(43);
  txt->SetTextSize(22);
  txt->Draw();
}
