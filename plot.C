#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "THStack.h"
#include "TLegend.h"

void plot(){
  //gROOT->LoadMacro("tdrstyle.C");
  //setTDRStyle();

  TFile *f = new TFile("multiplehiggssamplefitold.root","READ");

  TH2D *h_chi2comparehfzf = (TH2D*) f->Get("h_chi2comparehfzf");
  TH2D *h_chi2comparehflikelihoodzflikelihood = (TH2D*) f->Get("h_chi2comparehflikelihoodzflikelihood");
  TH2D *h_chi2comparehfhflikelihood = (TH2D*) f->Get("h_chi2comparehfhflikelihood");
  TH2D *h_chi2comparehfzflikelihood = (TH2D*) f->Get("h_chi2comparehfzflikelihood");


  TH2D *h_hfLhflikelihoodL = (TH2D*) f->Get("h_hfLhflikelihoodL");
  TH2D *h_zfLzflikelihoodL = (TH2D*) f->Get("h_zfLzflikelihoodL");
  TH2D *h_hfLkdkincompareLkddyn = (TH2D*) f->Get("h_hfLkdkincompareLkddyn");
  TH2D *h_zfLkdkincompareLkddyn = (TH2D*) f->Get("h_zfLkdkincompareLkddyn");
  TH2D *h_EkcompareEkdtau1h = (TH2D*) f->Get("h_EkcompareEkdtau1h");
  TH2D *h_EkcompareEkdtau2h = (TH2D*) f->Get("h_EkcompareEkdtau2h");

  TH1D *h_Eresulutionkhtau1 = (TH1D*) f->Get("h_Eresulutionkhtau1");
  TH1D *h_Eresulutionkdhtau1 = (TH1D*) f->Get("h_Eresulutionkdhtau1");


  THStack *Eresulutionkhcomparekdhtau1= new THStack("Eresulutionkhcomparekdhtau1","energyresulution from tau1 in the k-fit(red) and kd-fit(blue)");
  Eresulutionkhcomparekdhtau1->Add(h_Eresulutionkhtau1);
  h_Eresulutionkhtau1->SetLineColor(kRed);
  Eresulutionkhcomparekdhtau1->Add(h_Eresulutionkdhtau1);




  /*TCanvas *c_hfLhflikelihoodL = new TCanvas("h_hfLhflikelihoodL","h_hfLhflikelihoodL",1000,1000);
  h_hfLhflikelihoodL->Draw("colz");
  c_hfLhflikelihoodL->Print("histos.pdf[");
  c_hfLhflikelihoodL->Print("histos.pdf");
  c_hfLhflikelihoodL->Close();

  TCanvas *c_zfLzflikelihoodL = new TCanvas("zfLzflikelihoodL","zfLzflikelihoodL",1000,1000);
  h_zfLzflikelihoodL->Draw("colz");
  c_zfLzflikelihoodL->Print("histos.pdf");
  c_zfLzflikelihoodL->Close();*/


  TCanvas *c_Eresulutionkhcomparekdhtau1 = new TCanvas("Eresulutionkhtau1","Eresulutionkhtau1",2000,2000);
  Eresulutionkhcomparekdhtau1->Draw("nostack");
  Eresulutionkhcomparekdhtau1->GetXaxis()->SetTitle("(E-EFit)/E");
  Eresulutionkhcomparekdhtau1->GetYaxis()->SetTitle("Counts");
  c_Eresulutionkhcomparekdhtau1->Print("energyresulution.jpg");
  c_Eresulutionkhcomparekdhtau1->Print("histos.pdf[");
  c_Eresulutionkhcomparekdhtau1->Print("histos.pdf");
  c_Eresulutionkhcomparekdhtau1->Close();

  TCanvas *c_hfLkdkincompareLkddyn = new TCanvas("hfLkdkincompareLkddyn","hfLkdkincompareLkddyn",2400,2000);
  h_hfLkdkincompareLkddyn->Draw("colz");
  h_hfLkdkincompareLkddyn->GetYaxis()->SetTitle("PDF1(z1)*PDF2(z2)");
  h_hfLkdkincompareLkddyn->GetXaxis()->SetTitle("N*exp(-#chi^2/2)");
  gStyle->SetOptStat(0);
  h_hfLkdkincompareLkddyn->Draw("colz");
  h_hfLkdkincompareLkddyn->SetMaximum(600);
  c_hfLkdkincompareLkddyn->SetLogz();
  c_hfLkdkincompareLkddyn->Print("higgslkvslkd.jpg");
  c_hfLkdkincompareLkddyn->Print("histos.pdf");
  c_hfLkdkincompareLkddyn->Close();

  TCanvas *c_zfLkdkincompareLkddyn = new TCanvas("zfLkdkincompareLkddyn","zfLkdkincompareLkddyn",2400,2000);
  h_zfLkdkincompareLkddyn->GetYaxis()->SetTitle("PDF1(z1)*PDF2(z2)");
  h_zfLkdkincompareLkddyn->GetXaxis()->SetTitle("N*exp(-#chi^2/2)");
  h_zfLkdkincompareLkddyn->Draw("colz");
  h_zfLkdkincompareLkddyn->SetMaximum(600);
  h_zfLkdkincompareLkddyn->SetMinimum(0);
  c_zfLkdkincompareLkddyn->SetLogz();
  c_zfLkdkincompareLkddyn->Print("z.jpg");
  c_zfLkdkincompareLkddyn->Print("histos.pdf");
  c_zfLkdkincompareLkddyn->Close();

  TCanvas *c_hfLkdkincompareLkddyn2 = new TCanvas("hfLkdkincompareLkddyn","hfLkdkincompareLkddyn",1000,1000);
  h_hfLkdkincompareLkddyn->GetXaxis()->SetRangeUser(0.001,0.00125);
  h_hfLkdkincompareLkddyn->Draw("colz");
  h_hfLkdkincompareLkddyn->SetMaximum(600);
  h_hfLkdkincompareLkddyn->SetMinimum(0);
  c_hfLkdkincompareLkddyn2->SetLogz();
  c_hfLkdkincompareLkddyn2->Print("histos.pdf");
  c_hfLkdkincompareLkddyn2->Close();

  TCanvas *c_zfLkdkincompareLkddyn2 = new TCanvas("zfLkdkincompareLkddyn","zfLkdkincompareLkddyn",1000,1000);
  h_zfLkdkincompareLkddyn->GetXaxis()->SetRangeUser(0.001,0.00125);
  h_zfLkdkincompareLkddyn->Draw("colz");
  h_zfLkdkincompareLkddyn->SetMaximum(600);
  c_zfLkdkincompareLkddyn2->SetLogz();
  c_zfLkdkincompareLkddyn2->Print("histos.pdf");
  c_zfLkdkincompareLkddyn2->Close();

  TCanvas *c_EkcompareEkdtau1h = new TCanvas("EkcompareEkdtau1h","EkcompareEkdtau1h",1000,1000);
  h_EkcompareEkdtau1h->Draw("colz");
  c_EkcompareEkdtau1h->Print("histos.pdf");
  c_EkcompareEkdtau1h->Close();

  TCanvas *c_EkcompareEkdtau2h = new TCanvas("EkcompareEkdtau2h","EkcompareEkdtau2h",1000,1000);
  h_EkcompareEkdtau2h->Draw("colz");
  c_EkcompareEkdtau2h->Print("histos.pdf");
  c_EkcompareEkdtau2h->Close();

  TCanvas *c_EkcompareEkd = new TCanvas("EkcompareEkd","EkcompareEkd",2000,1000);
 // h_EkcompareEkdtau1h->Draw("colz");
 // h_EkcompareEkdtau2h->Draw("colz,same");
  TH2D *h_EkcompareEkd = (TH2D*) f->Get("h_EkcompareEkdtau1h");
  h_EkcompareEkd->Add(h_EkcompareEkdtau2h);
  h_EkcompareEkd->Draw("colz");
  h_EkcompareEkd->SetTitle("Final Tau-Energys from k-fit compared with kd-fit");
  h_EkcompareEkd->GetXaxis()->SetTitle("E_k[GeV]");
  h_EkcompareEkd->GetYaxis()->SetTitle("E_kd[GeV]");
  h_EkcompareEkd->Draw("colz");
  c_EkcompareEkd->Print("histos.pdf");
  c_EkcompareEkd->Print("EkcompareEkd.jpg");
  c_EkcompareEkd->Close();



  TCanvas *c_chi2comparehfhflikelihood = new TCanvas("chi2comparehfhflikelihood","chi2comparehfhflikelihood",2400,2000);
  h_chi2comparehfhflikelihood->Draw("colz");
  h_chi2comparehfhflikelihood->GetXaxis()->SetTitle("#chi^2-k-fit");
  h_chi2comparehfhflikelihood->GetYaxis()->SetTitle("#chi^2-kd-fit");
  h_chi2comparehfhflikelihood->SetTitle("compare chi2 from the k-fit and kd-fit");
  c_chi2comparehfhflikelihood->Print("chi2compare.jpg");
  c_chi2comparehfhflikelihood->Print("histos.pdf");
  c_chi2comparehfhflikelihood->Print("histos.pdf]");
  c_chi2comparehfhflikelihood->Close();















}
