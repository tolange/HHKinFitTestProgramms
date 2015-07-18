#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TCanvas.h"
#include "TPaveStats.h"
#include "THStack.h"
#include "TLegend.h"
#include "TTree.h"



void plotnew(){

TFile *f1 = new TFile("newToyMctest.root","READ");

TTree *tree1 = (TTree*) f1->Get("HHTauTau_tree");

TH1D* eresulutiontau1k= new TH1D("eresulutiontau1k","eresulutiontau1k",100,-1,1);
tree1->Draw("(Energy_fromTau1_truth-Energy_fromTau1_kfit)/Energy_fromTau1_truth>>eresulutiontau1k","abs((Energy_fromTau1_truth-Energy_fromTau1_kfit)/Energy_fromTau1_truth)<1","");
TH1D* eresulutiontau2k= new TH1D("eresulutiontau2k","eresulutiontau2k",100,-1,1);
tree1->Draw("(Energy_fromTau2_truth-Energy_fromTau2_kfit)/Energy_fromTau2_truth>>eresulutiontau2k","abs((Energy_fromTau2_truth-Energy_fromTau2_kfit)/Energy_fromTau2_truth)<1","");
TH1D* eresulutiontauk= new TH1D("eresulutiontauk","eresulutiontauk",100,-1,1);
eresulutiontauk->Add(eresulutiontau1k);
eresulutiontauk->Add(eresulutiontau2k);

TH1D* eresulutiontau1kd= new TH1D("eresulutiontau1kd","eresulutiontau1kd",100,-1,1);
tree1->Draw("(Energy_fromTau1_truth-Energy_fromTau1_kdfit)/Energy_fromTau1_truth>>eresulutiontau1kd","abs((Energy_fromTau1_truth-Energy_fromTau1_kdfit)/Energy_fromTau1_truth)<1","");
TH1D* eresulutiontau2kd= new TH1D("eresulutiontau2kd","eresulutiontau2kd",100,-1,1);
tree1->Draw("(Energy_fromTau2_truth-Energy_fromTau2_kdfit)/Energy_fromTau2_truth>>eresulutiontau2kd","abs((Energy_fromTau2_truth-Energy_fromTau2_kdfit)/Energy_fromTau2_truth)<1","");
TH1D* eresulutiontaukd= new TH1D("eresulutiontaukd","eresulutiontaukd",100,-1,1);
eresulutiontaukd->Add(eresulutiontau1kd);
eresulutiontaukd->Add(eresulutiontau2kd);

TH1D* deltachikddivdeltachikplus = new TH1D("deltachikddivdeltachikplus","deltachikddivdeltachikplus",100,-1,0);
tree1->Draw("(deltachi2kddynamicplusz/deltachi2kdkinematicplusz)>>deltachikddivdeltachikplus","","");
TH1D* deltachikddivdeltachikminus = new TH1D("deltachikddivdeltachikminus","deltachikddivdeltachikminus",100,-5,5);
tree1->Draw("deltachi2kddynamicminusz/deltachi2kdkinematicminusz>>deltachikddivdeltachikminus","abs(deltachi2kddynamicminusz/deltachi2kdkinematicminusz)<5","");


TH1D* z1divAkd =new TH1D("z1divAkd","z1divAkd",100,0,50);
tree1->Draw("(visible_Energy_fromTau1_vis/Energy_fromTau1_kdfit)/constantz1dotz2_kdfit>>z1divAkd","(visible_Energy_fromTau1_vis/Energy_fromTau1_kdfit)/constantz1dotz2_kdfit<50","");

TH1D* deltazkvskddivsigma=new TH1D("deltazkvskddivsigma","deltazkvskddivsigma",100,-1,0);
tree1->Draw("deltazkvskddivsigma>>deltazkvskddivsigma","a!=0","");

TH1D* deltazkvskddivsigmahighzcut=new TH1D("deltazkvskddivsigmahighzcut","deltazkvskddivsigmahighzcut",100,-1,0);
tree1->Draw("deltazkvskddivsigma>>deltazkvskddivsigmahighzcut","ztau1truth>0.85&&a!=0","");

TH1D* deltazkvskddivsigmalowzcut=new TH1D("deltazkvskddivsigmalowzcut","deltazkvskddivsigmalowzcut",100,-1,0);
tree1->Draw("deltazkvskddivsigma>>deltazkvskddivsigmalowzcut","ztau1truth<0.3&&a!=0","");

TH1D* deltazkvskddivsigmacutchi2k0_03=new TH1D("deltazkvskddivsigmacutchi2k0_03","deltazkvskddivsigmacutchi2k0_03",100,-1,0);
tree1->Draw("deltazkvskddivsigma>>deltazkvskddivsigmacutchi2k0_03","a!=0&&kinematic_part_of_the_final_chi2_kdfit<0.3","");

TH1D* deltazkvskddivsigmacutchi2k03_07=new TH1D("deltazkvskddivsigmacutchi2k03_07","deltazkvskddivsigmacutchi2k03_07",100,-1,0);
tree1->Draw("deltazkvskddivsigma>>deltazkvskddivsigmacutchi2k03_07","a!=0&&kinematic_part_of_the_final_chi2_kdfit>0.3&&kinematic_part_of_the_final_chi2_kdfit<0.7","");

TH1D* deltazkvskddivsigmacutchi2k07_1=new TH1D("deltazkvskddivsigmacutchi2k07_1","deltazkvskddivsigmacutchi2k07_1",100,-1,0);
tree1->Draw("deltazkvskddivsigma>>deltazkvskddivsigmacutchi2k07_1","a!=0&&kinematic_part_of_the_final_chi2_kdfit>0.7&&kinematic_part_of_the_final_chi2_kdfit<1","");

TH1D* deltazkvskddivsigmacutchi2k1_15=new TH1D("deltazkvskddivsigmacutchi2k1_15","deltazkvskddivsigmacutchi2k1_15",100,-1,0);
tree1->Draw("deltazkvskddivsigma>>deltazkvskddivsigmacutchi2k1_15","a!=0&&kinematic_part_of_the_final_chi2_kdfit>1&&kinematic_part_of_the_final_chi2_kdfit<1.5","");

/////////////////////////////////////////////////

TH2D* sigmazvsz =new TH2D("sigmazvsz","sigmazvsz",100,0,1,100,0,1);
tree1->Draw("deltaz:ztau1truth>>sigmazvsz","a!=0","colz");

TH2D* deltazvsz=new TH2D("deltazvsz","deltazvsz",100,-1,1,100,0,1);
tree1->Draw("(deltazkvskddivsigma*deltaz):ztau1truth>>deltazvsz","a!=0","colz");

TH2D* deltazvssimaz=new TH2D("deltazvssimaz","deltazvssimaz",100,-1,1,100,0,1);
tree1->Draw("(deltazkvskddivsigma*deltaz):deltaz>>deltazvssimaz","a!=0","colz");

TH1D* zkminztruth=new TH1D("zkminztruth","zkminztruth",100,-1,1);
tree1->Draw("(visible_Energy_fromTau1_vis/Energy_fromTau1_kfit-ztau1truth)>>zkminztruth","","");

TH1D* zkminztruthdivsigmaz=new TH1D("zkminztruthdivsigmaz","zkminztruthdivsigmaz",100,-1,1);
tree1->Draw("(visible_Energy_fromTau1_vis/Energy_fromTau1_kfit-ztau1truth)/deltaz>>zkminztruthdivsigmaz","a!=0","");

TH1D* zkdminztruth=new TH1D("zkdminztruth","zkdminztruth",100,-1,1);
tree1->Draw("(visible_Energy_fromTau1_vis/Energy_fromTau1_kdfit-ztau1truth)>>zkdminztruth","","");

TH1D* zkdminztruthdivsigmaz=new TH1D("zkdminztruthdivsigmaz","zkdminztruthdivsigmaz",100,-1,1);
tree1->Draw("(visible_Energy_fromTau1_vis/Energy_fromTau1_kdfit-ztau1truth)/deltaz>>zkdminztruthdivsigmaz","a!=0","");



TCanvas* c1 = new TCanvas("newToyMcStudies1","newToyMcStudies1",1500,1000);
eresulutiontau1k->SetTitle("Energyresulution from Tau 1 (kinematic Fit)");
eresulutiontau1k->GetXaxis()->SetTitle("#frac{E-E_{F}}{E}");
eresulutiontau1k->GetYaxis()->SetTitleOffset(1.2);
eresulutiontau1k->GetYaxis()->SetTitle("Counts");
eresulutiontau1k->Draw();
c1->Print("newToyMcStudies.pdf[");
c1->Print("newToyMcStudies.pdf");
c1->Close();
TCanvas* c2 = new TCanvas("newToyMcStudies2","newToyMcStudies2",1500,1000);
eresulutiontau2k->SetTitle("Energyresulution from Tau 1(blue) and Tau 2(red)  (kinematic Fit)");
eresulutiontau2k->GetXaxis()->SetTitle("#frac{E-E_{F}}{E}");
eresulutiontau2k->GetYaxis()->SetTitleOffset(1.2);
eresulutiontau2k->GetYaxis()->SetTitle("Counts");
eresulutiontau2k->SetLineColor(kRed);
eresulutiontau2k->Draw();
eresulutiontau1k->Draw("same");
c2->Print("newToyMcStudies.pdf");
c2->Close();
TCanvas* c3 = new TCanvas("newToyMcStudie3","newToyMcStudie3",1500,1000);
eresulutiontauk->SetTitle("Energyresulution (kinematic Fit)");
eresulutiontauk->GetXaxis()->SetTitle("#frac{E-E_{F}}{E}");
eresulutiontauk->GetYaxis()->SetTitleOffset(1.2);
eresulutiontauk->GetYaxis()->SetTitleOffset(1.5);
eresulutiontauk->GetYaxis()->SetTitle("Counts");
eresulutiontauk->Draw();
c3->Print("newToyMcStudies.pdf");
c3->Close();
TCanvas* c4 = new TCanvas("newToyMcStudies4","newToyMcStudies4",1500,1000);
eresulutiontaukd->SetTitle("Energyresulution from the kinematic Fit(blue) and dynamic Fit (red)");
eresulutiontaukd->GetXaxis()->SetTitle("#frac{E-E_{F}}{E}");
eresulutiontaukd->GetYaxis()->SetTitleOffset(1.2);
eresulutiontaukd->GetYaxis()->SetTitleOffset(1.2);
eresulutiontaukd->GetYaxis()->SetTitle("Counts");
eresulutiontaukd->SetLineColor(kRed);
eresulutiontaukd->Draw();
eresulutiontauk->Draw("same");
c4->Print("newToyMcStudies.pdf");
c4->Close();
TF1* restauk =new TF1("restauk","gaus",-0.2,0.2);
eresulutiontauk->Fit(restauk,"R");
TF1* restaukd =new TF1("restaukd","gaus",-0.2,0.2);
eresulutiontaukd->Fit(restaukd,"R");
TCanvas* c5 = new TCanvas("newToyMcStudies5","Energyresulution from the kinematic Fit(blue) and dynamic Fit (red)",1500,1000);
TLegend* leg =new TLegend (0,0,0.3,0.3);
leg->AddEntry(restauk,"kinematic-fit","l");
leg->AddEntry((TObject*)0,"Mean: (2.071#pm0.271)*10^{-3}","");
leg->AddEntry((TObject*)0,"#sigma: (1.00545#pm0.00275)*10^{-1}","");
leg->AddEntry(restaukd,"dynamic-fit","l");
leg->AddEntry((TObject*)0,"Mean: (1.970#pm0.263)*10^{-3}","");
leg->AddEntry((TObject*)0,"#sigma: (0.98891#pm0.00264)*10^{-1}","");
leg->AddEntry((TObject*)0,"#frac{#Delta #sigma}{#sigma}: (1.645#pm 0.368)%","");
restauk->SetLineWidth(0.5);
restauk->SetTitle("Energyresulution from kinematic-(blue) and dynamic-fit(red) after an gaussian-fit ");
restauk->GetXaxis()->SetTitle("#frac{E-E_{F}}{E}");
restauk->GetYaxis()->SetTitleOffset(1.2);
restauk->GetXaxis()->SetTitleOffset(1.2);
restauk->GetYaxis()->SetTitle("Counts");
restauk->SetLineColor(kBlue);
restauk->Draw();
restaukd->SetLineWidth(0.5);
restaukd->Draw("same");
leg->Draw();
c5->Print("newToyMcStudies.pdf");
c5->Close();
TCanvas* c6 = new TCanvas("newToyMcStudies6","newToyMcStudies6",1500,1000);
z1divAkd->SetTitle("#frac{z_{1}}{A} ,z= #frac{E_{vis}}{E}, A=z_{1}*z_{2}   (#chi_{d}^{2}=-2*ln(4*(z_{1}-A)))");
z1divAkd->GetXaxis()->SetTitleOffset(1.15);
z1divAkd->GetXaxis()->SetTitle("#frac{z_{1}}{A}");
z1divAkd->GetYaxis()->SetTitleOffset(1.25);
z1divAkd->GetYaxis()->SetTitle("Counts");
z1divAkd->Draw();
c6->Print("newToyMcStudies.pdf");
c6->Close();
TCanvas* c7 = new TCanvas("newToyMcStudies7","newToyMcStudies7",1350,1000);
c7->SetMargin(0.325,0.05,0.3,0.2);
deltachikddivdeltachikplus->SetTitleOffset(2);
deltachikddivdeltachikplus->SetTitle("#frac{#Delta #chi_{d}^{2}}{#Delta #chi_{k}^{2}}for #frac{#Delta z}{z}=+10%");
deltachikddivdeltachikplus->GetXaxis()->SetTitleOffset(1.5);
deltachikddivdeltachikplus->GetXaxis()->SetTitle("#frac{#Delta #chi_{d}^{2}}{#Delta #chi_{k}^{2}}");
deltachikddivdeltachikplus->GetYaxis()->SetTitleOffset(1.25);
deltachikddivdeltachikplus->GetYaxis()->SetTitle("Counts");
deltachikddivdeltachikplus->Draw();
c7->Print("newToyMcStudies.pdf");
c7->Close();
TCanvas* c8 = new TCanvas("newToyMcStudie8","newToyMcStudies8",1350,1000);
c8->SetMargin(0.325,0.05,0.3,0.2);
deltachikddivdeltachikminus->SetTitleOffset(2);
deltachikddivdeltachikminus->SetTitle("#frac{#Delta #chi_{d}^{2}}{#Delta #chi_{k}^{2}} for #frac{#Delta z}{z}=-10%");
deltachikddivdeltachikminus->GetXaxis()->SetTitleOffset(1.5);
deltachikddivdeltachikminus->GetXaxis()->SetTitle("#frac{#Delta #chi_{d}^{2}}{#Delta #chi_{k}^{2}}");
deltachikddivdeltachikminus->GetYaxis()->SetTitleOffset(1.3);
deltachikddivdeltachikminus->GetYaxis()->SetTitle("Counts");
deltachikddivdeltachikminus->Draw();
c8->Print("newToyMcStudies.pdf");
c8->Close();
TCanvas* c9= new TCanvas("newToyMcStudie9","newToyMcStudies9",1500,1000);
c9->SetBottomMargin(0.2);
deltazkvskddivsigma->SetTitle("Impact of the dynamic constraint on the final energyfraction");
deltazkvskddivsigma->GetXaxis()->SetTitle("#frac{z_{k}-z_{kd}}{#sigma}");
deltazkvskddivsigma->GetXaxis()->SetTitleOffset(1.2);
deltazkvskddivsigma->GetYaxis()->SetTitle("Counts");
deltazkvskddivsigma->Draw();
c9->Print("newToyMcStudies.pdf");
c9->Close();
TCanvas* c10= new TCanvas("newToyMcStudie10","newToyMcStudies10",1500,1000);
c10->SetBottomMargin(0.2);
deltazkvskddivsigmahighzcut->SetTitle("Impact of the dynamic constraint on the final energyfraction (z>0.85)");
deltazkvskddivsigmahighzcut->GetXaxis()->SetTitle("#frac{z_{k}-z_{kd}}{#sigma}");
deltazkvskddivsigmahighzcut->GetXaxis()->SetTitleOffset(1.2);
deltazkvskddivsigmahighzcut->GetYaxis()->SetTitle("Counts");
deltazkvskddivsigmahighzcut->Draw();
c10->Print("newToyMcStudies.pdf");
c10->Close();
TCanvas* c11= new TCanvas("newToyMcStudie11","newToyMcStudies11",1500,1000);
c11->SetBottomMargin(0.2);
deltazkvskddivsigmalowzcut->SetTitle("Impact of the dynamic constraint on the final energyfraction (z<0.3)");
deltazkvskddivsigmalowzcut->GetXaxis()->SetTitle("#frac{z_{k}-z_{kd}}{#sigma}");
deltazkvskddivsigmalowzcut->GetXaxis()->SetTitleOffset(1.2);
deltazkvskddivsigmalowzcut->GetYaxis()->SetTitle("Counts");
deltazkvskddivsigmalowzcut->Draw();
c11->Print("newToyMcStudies.pdf");
c11->Close();


TCanvas* c12= new TCanvas("newToyMcStudie12","newToyMcStudies12",1500,1000);
c12->SetBottomMargin(0.2);
deltazkvskddivsigmacutchi2k0_03->SetTitle("Impact of the dynamic constraint on the final energyfraction (kinematic part of #chi_{kd}^{2}<0.3)");
deltazkvskddivsigmacutchi2k0_03->GetXaxis()->SetTitle("#frac{z_{k}-z_{kd}}{#sigma}");
deltazkvskddivsigmacutchi2k0_03->GetXaxis()->SetTitleOffset(1.2);
deltazkvskddivsigmacutchi2k0_03->GetYaxis()->SetTitle("Counts");
deltazkvskddivsigmacutchi2k0_03->Draw();
c12->Print("newToyMcStudies.pdf");
c12->Close();

TCanvas* c13= new TCanvas("newToyMcStudie13","newToyMcStudies13",1500,1000);
c13->SetBottomMargin(0.2);
deltazkvskddivsigmacutchi2k03_07->SetTitle("Impact of the dynamic constraint on the final energyfraction (kinematic part of #chi_{kd}^{2}>0.3&&<0.7)");
deltazkvskddivsigmacutchi2k03_07->GetXaxis()->SetTitle("#frac{z_{k}-z_{kd}}{#sigma}");
deltazkvskddivsigmacutchi2k03_07->GetXaxis()->SetTitleOffset(1.2);
deltazkvskddivsigmacutchi2k03_07->GetYaxis()->SetTitle("Counts");
deltazkvskddivsigmacutchi2k03_07->Draw();
c13->Print("newToyMcStudies.pdf");
c13->Close();

TCanvas* c14= new TCanvas("newToyMcStudie14","newToyMcStudies14",1500,1000);
c14->SetBottomMargin(0.2);
deltazkvskddivsigmacutchi2k07_1->SetTitle("Impact of the dynamic constraint on the final energyfraction (kinematic part of #chi_{kd}^{2}>0.7&&<1)");
deltazkvskddivsigmacutchi2k07_1->GetXaxis()->SetTitle("#frac{z_{k}-z_{kd}}{#sigma}");
deltazkvskddivsigmacutchi2k07_1->GetXaxis()->SetTitleOffset(1.2);
deltazkvskddivsigmacutchi2k07_1->GetYaxis()->SetTitle("Counts");
deltazkvskddivsigmacutchi2k07_1->Draw();
c14->Print("newToyMcStudies.pdf");
c14->Close();

TCanvas* c15= new TCanvas("newToyMcStudie15","newToyMcStudies15",1500,1000);
c15->SetBottomMargin(0.2);
deltazkvskddivsigmacutchi2k1_15->SetTitle("Impact of the dynamic constraint on the final energyfraction (kinematic part of #chi_{kd}^{2}>1&&<1.5)");
deltazkvskddivsigmacutchi2k1_15->GetXaxis()->SetTitle("#frac{z_{k}-z_{kd}}{#sigma}");
deltazkvskddivsigmacutchi2k1_15->GetXaxis()->SetTitleOffset(1.2);
deltazkvskddivsigmacutchi2k1_15->GetYaxis()->SetTitle("Counts");
deltazkvskddivsigmacutchi2k1_15->Draw();
c15->Print("newToyMcStudies.pdf");
c15->Close();

/////////////////////////////////////////////////////

TCanvas* c16= new TCanvas("newToyMcStudie16","newToyMcStudies16",1500,1000);
c16->SetBottomMargin(0.2);
sigmazvsz->SetTitle("The Error from z_{1} against z");
sigmazvsz->GetXaxis()->SetTitle("#sigma_{z}");
sigmazvsz->GetXaxis()->SetTitleOffset(1.2);
sigmazvsz->GetYaxis()->SetTitle("#z=frac{#E_{#tau,vis}}{#E_{#tau,Fit}}");
sigmazvsz->Draw("colz");
c16->Print("newToyMcStudies.pdf");
c16->Close();

TCanvas* c17= new TCanvas("newToyMcStudie17","newToyMcStudies17",1500,1000);
c17->SetBottomMargin(0.2);
deltazvsz->SetTitle("Impact of the dynamic constraint on the final energyfraction z against z");
deltazvsz->GetXaxis()->SetTitle("z_{k}-z_{kd}");
deltazvsz->GetXaxis()->SetTitleOffset(1.2);
deltazvsz->GetYaxis()->SetTitle("#z=frac{#E_{#tau,vis}}{#E_{#tau,Fit}}");
deltazvsz->Draw("colz");
c17->Print("newToyMcStudies.pdf");
c17->Close();

TCanvas* c18= new TCanvas("newToyMcStudie18","newToyMcStudies18",1500,1000);
c18->SetBottomMargin(0.2);
deltazvssimaz->SetTitle("Impact of the dynamic constraint on the final energyfraction z against the errot on z");
deltazvssimaz->GetXaxis()->SetTitle("z_{k}-z_{kd}");
deltazvssimaz->GetXaxis()->SetTitleOffset(1.2);
deltazvssimaz->GetYaxis()->SetTitle("#sigma_{z}");
deltazvssimaz->Draw("colz");
c18->Print("newToyMcStudies.pdf");
c18->Close();

TCanvas* c19= new TCanvas("newToyMcStudie19","newToyMcStudies19",1500,1000);
c19->SetBottomMargin(0.2);
zkminztruth->SetTitle("Difference of the energyfraction z_{k} out of the k-fit and z_{truth}");
zkminztruth->GetXaxis()->SetTitle("z_{k}-z_{truth}");
zkminztruth->GetXaxis()->SetTitleOffset(1.2);
zkminztruth->GetYaxis()->SetTitle("Counts");
zkminztruth->Draw("");
c19->Print("newToyMcStudies.pdf");
c19->Close();

TCanvas* c20= new TCanvas("newToyMcStudie20","newToyMcStudies20",1500,1000);
c20->SetBottomMargin(0.2);
zkminztruthdivsigmaz->SetTitle("Difference of the energyfraction z_{k} out of the k-fit and z_{truth}");
zkminztruthdivsigmaz->GetXaxis()->SetTitle("frac{z_{k}-z_{truth}}{#sigma_{z}}");
zkminztruthdivsigmaz->GetXaxis()->SetTitleOffset(1.2);
zkminztruthdivsigmaz->GetYaxis()->SetTitle("Counts");
zkminztruthdivsigmaz->Draw("");
c20->Print("newToyMcStudies.pdf");
c20->Close();

TCanvas* c21= new TCanvas("newToyMcStudie21","newToyMcStudies21",1500,1000);
c21->SetBottomMargin(0.2);
zkdminztruth->SetTitle("Difference of the energyfraction z_{kd} out of the k-fit and z_{truth}");
zkdminztruth->GetXaxis()->SetTitle("z_{kd}-z_{truth}");
zkdminztruth->GetXaxis()->SetTitleOffset(1.2);
zkdminztruth->GetYaxis()->SetTitle("Counts");
zkdminztruth->Draw("");
c21->Print("newToyMcStudies.pdf");
c21->Close();

TCanvas* c22= new TCanvas("newToyMcStudie20","newToyMcStudies20",1500,1000);
c22->SetBottomMargin(0.2);
zkdminztruthdivsigmaz->SetTitle("Difference of the energyfraction z_{kd} out of the k-fit and z_{truth}");
zkdminztruthdivsigmaz->GetXaxis()->SetTitle("frac{z_{kd}-z_{truth}}{#sigma_{z}}");
zkdminztruthdivsigmaz->GetXaxis()->SetTitleOffset(1.2);
zkdminztruthdivsigmaz->GetYaxis()->SetTitle("Counts");
zkdminztruthdivsigmaz->Draw("");
c22->Print("newToyMcStudies.pdf");
c22->Print("newToyMcStudies.pdf]");
c22->Close();










}





