#include "HHTauTauEventGenerator.h"
#include "TFile.h"
#include "TH1D.h"
#include "TDirectory.h"
#include <iostream>
#include "HHKinFit.h"
#include "HHFitConstraint4Vector.h"
#include "HHFitConstraintEHardM.h"
#include "HHFitConstraint.h"
#include "HHFitConstraintLikelihood.h"
#include "HHFitObjectEConstM.h"
#include "HHFitObjectE.h"
#include "HHFitObjectMET.h"
#include "HHFitObject.h"
#include "HHFitObjectComposite.h"
#include "HHLorentzVector.h"
#include "exceptions/HHCovarianceMatrixException.h"
#include "exceptions/HHEnergyRangeException.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TFile.h"
#include "HHFitConstraintLikelihood.h"
#include "TMath.h"
#include "TH2D.h"
#include "TGraph.h"
#include "TRandom3.h"
#include <cmath>
#include "TTree.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLine.h"
#include"TLegend.h"


using namespace HHKinFit2;

int main(int argc, char* argv[])
{
  // creating outfile
  TFile* outfile=new TFile("newToyMctest.root","RECREATE");

  //Define Tree
  TTree * HHTauTau_tree=new TTree("HHTauTau_tree","HHTauTau_tree");
  HHTauTau_tree->SetDirectory(0);
  HHTauTau_tree->SetAutoSave(100000000);


  unsigned int j;
  HHTauTau_tree->Branch("eventnumber",&j);
  double etau1truth;
  HHTauTau_tree->Branch("Energy_fromTau1_truth",&etau1truth);
  double etau2truth;
  HHTauTau_tree->Branch("Energy_fromTau2_truth",&etau2truth);
  double etau1kfit;
  HHTauTau_tree->Branch("Energy_fromTau1_kfit",&etau1kfit);
  double etau2kfit;
  HHTauTau_tree->Branch("Energy_fromTau2_kfit",&etau2kfit);
  double etau1kdfit;
  HHTauTau_tree->Branch("Energy_fromTau1_kdfit",&etau1kdfit);
  double etau2kdfit;
  HHTauTau_tree->Branch("Energy_fromTau2_kdfit",&etau2kdfit);
  double etau1vis;
  HHTauTau_tree->Branch("visible_Energy_fromTau1_vis",&etau1vis);
  double etau2vis;
  HHTauTau_tree->Branch("visible_Energy_fromTau2_vis",&etau2vis);
  double phitau1;
  HHTauTau_tree->Branch("The_angle_phi_from_tau1",&phitau1);
  double phitau2;
  HHTauTau_tree->Branch("The_angle_phi_from_tau2",&phitau2);
  double etatau1;
  HHTauTau_tree->Branch("Eta_from_tau1",&etatau1);
  double etatau2;
  HHTauTau_tree->Branch("Eta_phi_from_tau2",&etatau2);
  double METx;
  HHTauTau_tree->Branch("Missing_energy_in_x_direction_METx",&METx);
  double METy;
  HHTauTau_tree->Branch("Missing_energy_in_y_direction_METy",&METy);
  double chi2k;
  HHTauTau_tree->Branch("final_chi2_kfit",&chi2k);
  double chi2kd;
  HHTauTau_tree->Branch("final_chi2_kdfit",&chi2kd);
  double chi2kinematick;
  HHTauTau_tree->Branch("kinematic_part_of_the_final_chi2_kfit",&chi2kinematick);
  double chi2kinematickd;
  HHTauTau_tree->Branch("kinematic_part_of_the_final_chi2_kdfit",&chi2kinematickd);
  double chi2dynamick;
  HHTauTau_tree->Branch("dynamic_part_of_the_final_chi2_kfit",&chi2dynamick);
  double chi2dynamickd;
  HHTauTau_tree->Branch("dynamic_part_of_the_final_chi2_kdfit",&chi2dynamickd);
  double Ak;
  HHTauTau_tree->Branch("constantz1dotz2_kfit",&Ak);
  double Akd;
  HHTauTau_tree->Branch("constantz1dotz2_kdfit",&Akd);
  double ztau1truth;
  HHTauTau_tree->Branch("ztau1truth",&ztau1truth);
  double deltachi2kdkinematicplusz;
  HHTauTau_tree->Branch("deltachi2kdkinematicplusz",&deltachi2kdkinematicplusz);
  double deltachi2kddynamicplusz;
  HHTauTau_tree->Branch("deltachi2kddynamicplusz",&deltachi2kddynamicplusz);
  double deltachi2kdkinematicminusz;
  HHTauTau_tree->Branch("deltachi2kdkinematicminusz",&deltachi2kdkinematicminusz);
  double deltachi2kddynamicminusz;
  HHTauTau_tree->Branch("deltachi2kddynamicminusz",&deltachi2kddynamicminusz);


  //chi2 parabel a*(x-b)+c

  double a;
  HHTauTau_tree->Branch("a",&a);
  double b;
  HHTauTau_tree->Branch("b",&b);
  double c;
  HHTauTau_tree->Branch("c",&c);
  double deltaz;
  HHTauTau_tree->Branch("deltaz",&deltaz);
  double deltazkvskddivsigma;
  HHTauTau_tree->Branch("deltazkvskddivsigma",&deltazkvskddivsigma);

  TF1 *PDF1 = new TF1("PDF1","2*x",0,1);
  PDF1->SetNpx(100000);
  TF1 *PDF2= new TF1("PDF2","2-2*x",0,1);
  PDF2->SetNpx(100000);
  TMatrixD covarmatrix(2,2);
  covarmatrix[0][0]=130;
  covarmatrix[0][1]=0;
  covarmatrix[1][0]=0;
  covarmatrix[1][1]=130;
  HHTauTauEventGenerator Higgsgenerator(PDF1,PDF2,covarmatrix);
  //Higgsgenerator.setMtau1(0);
  //Higgsgenerator.setMtau2(0);
unsigned int maxnumber=100000;

  for(unsigned int i=0; i<maxnumber; i++){
	  double progress =(i*1.0/maxnumber)*100;
	  if(progress/10==0.5||progress/10==1||progress/10==1.5||progress/10==2||progress/10==2.5||progress/10==3||progress/10==3.5||progress/10==4||progress/10==4.5||progress/10==5||progress/10==5.5||progress/10==6||progress/10==6.5||progress/10==7||progress/10==7.5||progress/10==8||progress/10==8.5||progress/10==9||progress/10==9.5){
		  std::cout << "progress: " << progress <<" %" <<std::endl;
	  }


    try{
	    Higgsgenerator.generateEvent();
	   }
	catch(const HHEnergyRangeException& e){
		/*std::cout << Higgsgenerator.getvisfrac1() << std::endl;
		std::cout << Higgsgenerator.getTau1().E() << std::endl;
		std::cout << Higgsgenerator.getTau1Vis().E() <<std::endl;
		std::cout << e.what() << std::endl;
		std::cout << "---------------------------------------" << std::endl;*/
	  	  i--;
	  	  continue;
	      }


	//K-Fit -Objects

	double mass=125.7;

	HHFitObjectE* tau1k = new HHFitObjectEConstM(Higgsgenerator.getTau1Vis());
	HHFitObjectE* tau2k = new HHFitObjectEConstM(Higgsgenerator.getTau2Vis());

	tau1k->setLowerFitLimitE(tau1k->getInitial4Vector());
	tau1k->setUpperFitLimitE(mass,tau2k->getInitial4Vector());
	tau2k->setLowerFitLimitE(tau2k->getInitial4Vector());
	tau2k->setUpperFitLimitE(mass,tau1k->getInitial4Vector());

	HHFitObjectMET* metk = new HHFitObjectMET(TVector2(Higgsgenerator.getMETwithsigma()[0],Higgsgenerator.getMETwithsigma()[1]));
    metk->setCovMatrix(Higgsgenerator.getCovarmatrix()[0][0],Higgsgenerator.getCovarmatrix()[1][1],Higgsgenerator.getCovarmatrix()[1][0]);

    HHFitObject* higgsk  = new HHFitObjectComposite(tau1k, tau2k, metk);

	//KD-Fit -Objects

	HHFitObjectE* tau1kd = new HHFitObjectEConstM(Higgsgenerator.getTau1Vis());
	HHFitObjectE* tau2kd = new HHFitObjectEConstM(Higgsgenerator.getTau2Vis());


    tau1kd->setLowerFitLimitE(tau1kd->getInitial4Vector());
	tau1kd->setUpperFitLimitE(mass,tau2kd->getInitial4Vector());
	tau2kd->setLowerFitLimitE(tau2kd->getInitial4Vector());
	tau2kd->setUpperFitLimitE(mass,tau1kd->getInitial4Vector());



	HHFitObjectMET* metkd = new HHFitObjectMET(TVector2(Higgsgenerator.getMETwithsigma()[0],Higgsgenerator.getMETwithsigma()[1]));
	metkd->setCovMatrix(Higgsgenerator.getCovarmatrix()[0][0],Higgsgenerator.getCovarmatrix()[1][1],Higgsgenerator.getCovarmatrix()[1][0]);

	HHFitObject* higgskd  = new HHFitObjectComposite(tau1kd, tau2kd, metkd);


	//K-Fit -Constraints

	HHFitConstraint* invmk = new HHFitConstraintEHardM(tau1k, tau2k, mass);
	HHFitConstraint* balancek = new HHFitConstraint4Vector(higgsk, true, true, false, false);

	//KD-Fit -Constraints

	HHFitConstraint* invmkd = new HHFitConstraintEHardM(tau1kd, tau2kd, mass);
	HHFitConstraint* balancekd = new HHFitConstraint4Vector(higgskd, true, true, false, false);
	HHFitConstraint* Likelihoodkd = new HHFitConstraintLikelihood(tau1kd,tau2kd,PDF1,PDF2);


	//K-Fit -Fit

	HHKinFit* singlefitk = new HHKinFit();
	tau1k->setInitStart((tau1k->getUpperFitLimitE()+tau1k->getLowerFitLimitE())/2);
	tau1k->setInitPrecision(0.1);
	tau1k->setInitStepWidth(0.1*(tau1k->getUpperFitLimitE() - tau1k->getLowerFitLimitE()));
	tau1k->setInitDirection(1.0);

	singlefitk->addFitObjectE(tau1k);
    singlefitk->addConstraint(invmk);
    singlefitk->addConstraint(balancek);

    //KD-Fit -Fit

    HHKinFit* singlefitkd = new HHKinFit();
    tau1kd->setInitStart((tau1kd->getUpperFitLimitE()+tau1kd->getLowerFitLimitE())/2);
    tau1kd->setInitPrecision(0.1);
    tau1kd->setInitStepWidth(0.1*(tau1kd->getUpperFitLimitE() - tau1kd->getLowerFitLimitE()));
    tau1kd->setInitDirection(1.0);

    singlefitkd->addFitObjectE(tau1kd);
    singlefitkd->addConstraint(invmkd);
    singlefitkd->addConstraint(balancekd);
    singlefitkd->addConstraint(Likelihoodkd);


    //K-Fit -Fitting


    try {
  		 singlefitk->fit();
  		 /*if (!((singlefitk->getConvergence()==1)||(singlefitk->getConvergence()==2))) {
  		     i--;
  		     continue;
  		    }*/
  	   }
  	catch(HHEnergyRangeException const& e){
  		  i--;
  		  continue;
  	     }

  	//KD-Fit -Fitting

    try {
  	  	 singlefitkd->fit();
  	    /* if (!((singlefitkd->getConvergence()==1)||(singlefitkd->getConvergence()==2))) {
  	  		 i--;
  	  		 continue;
  	  		}*/
  		}
  	catch(HHEnergyRangeException const& e){
  	  	  i--;
  	  	  continue;
  	  	 }

  //Filling Tree


  j=i;
  etau1truth=Higgsgenerator.getTau1boosted().E();
  etau2truth=Higgsgenerator.getTau2boosted().E();
  etau1kfit=tau1k->getFit4Vector().E();
  etau2kfit=tau2k->getFit4Vector().E();
  etau1kdfit=tau1kd->getFit4Vector().E();
  etau2kdfit=tau2kd->getFit4Vector().E();
  etau1vis=tau1k->getInitial4Vector().E();
  etau2vis=tau2k->getInitial4Vector().E();
  phitau1=tau1k->getInitial4Vector().Phi();
  phitau2=tau2k->getInitial4Vector().Phi();
  etatau1=tau1k->getInitial4Vector().Eta();
  etatau2=tau2k->getInitial4Vector().Eta();
  METx=Higgsgenerator.getMETwithsigma()[0];
  METy=Higgsgenerator.getMETwithsigma()[0];
  chi2k=singlefitk->getChi2();
  chi2kd=singlefitkd->getChi2();
  chi2kinematick=singlefitk->getChi2();
  chi2kinematickd=balancekd->getChi2();
  chi2dynamick=(-2*log(PDF1->Eval((etau1truth-etau1kfit)/etau1truth))-2*log(PDF2->Eval((etau2truth-etau2kfit)/etau2truth)));
  chi2dynamickd=Likelihoodkd->getChi2();
  Ak=(etau1vis/etau1kfit)-exp(-(chi2dynamick+4*log(2))/2);
  Akd=(etau1vis/etau1kdfit)-exp(-(chi2dynamickd+4*log(2))/2);
  ztau1truth=Higgsgenerator.getvisfrac1();

  //std::cout  << "---------------------"<< i <<"----------------------------" << std::endl;
  //testing delta chik and delta chid
  double chikdkinematicold=balancekd->getChi2();
  double chikddynamicold=Likelihoodkd->getChi2();
 // std::cout << "kinematicchi2old: " << chikdkinematicold  << " dynamicold: " <<chikddynamicold << " sum: " <<chikddynamicold+chikdkinematicold<<std::endl;
  double etemp=tau1kd->getFit4Vector().E();
  double frac1oldd=tau1kd->getInitial4Vector().E()/tau1kd->getFit4Vector().E();
  if(frac1oldd<0.9){
	  tau1kd->changeEandSave((1/1.1)*tau1kd->getFit4Vector().E());
	  invmkd->getChi2();
	  deltachi2kdkinematicplusz=balancekd->getChi2()-chikdkinematicold;
	  deltachi2kddynamicplusz=Likelihoodkd->getChi2()-chikddynamicold;
	//  std::cout << "kinematic+: " << balancekd->getChi2()<< " dynamic+: " <<Likelihoodkd->getChi2() <<" sum: " << balancekd->getChi2()+Likelihoodkd->getChi2()<< std::endl;
	  tau1kd->changeEandSave(etemp);
	  invmk->getChi2();
  }
    if((1/0.9)*tau1kd->getFit4Vector().E()<tau1kd->getUpperFitLimitE()){
  	tau1kd->changeEandSave((1/0.9)*tau1kd->getFit4Vector().E());
  	invmkd->getChi2();
  	deltachi2kdkinematicminusz=balancekd->getChi2()-chikdkinematicold;
  	deltachi2kddynamicminusz=Likelihoodkd->getChi2()-chikddynamicold;
  	//std::cout << "kinematic-: " << balancekd->getChi2()<< " dynamic-: " <<Likelihoodkd->getChi2() <<" sum: " << balancekd->getChi2()+Likelihoodkd->getChi2() <<std::endl;
  	//std::cout  << "-------------------------------------------------" << std::endl;
  	tau1kd->changeEandSave(etemp);
  	invmk->getChi2();
    }


  double chi2kmin=singlefitk->getChi2();
  double zmin=tau1k->getInitial4Vector().E()/tau1k->getFit4Vector().E();
  double chi2plusten;
  double zplusten;
  double chi2minusten;
  double zminusten;
  double etemp2=tau1k->getFit4Vector().E();
  a=0;
  b=0;
  c=0;
  deltaz=0;
  deltazkvskddivsigma=0;


  if(zmin<0.9&&(1/0.9)*tau1k->getFit4Vector().E()<tau1k->getUpperFitLimitE()){
	  tau1k->changeEandSave((1/1.1)*tau1k->getFit4Vector().E());
	  invmk->getChi2();
	  chi2minusten=balancek->getChi2();
	  zminusten=tau1k->getInitial4Vector().E()/tau1k->getFit4Vector().E();
	  tau1k->changeEandSave(etemp2);
	  invmk->getChi2();
	  tau1k->changeEandSave((1/0.9)*tau1k->getFit4Vector().E());
	  invmk->getChi2();
	  chi2plusten=balancek->getChi2();
	  zplusten=tau1k->getInitial4Vector().E()/tau1k->getFit4Vector().E();
	  tau1k->changeEandSave(etemp2);
	  invmk->getChi2();
	  double x1=zmin;
	  double y1=chi2kmin;
	  double x2=zplusten;
	  double y2=chi2plusten;
	  double x3=zminusten;
	  double y3=chi2minusten;
	  /*b=(chi2plusten*(zminusten*zminusten-zmin*zmin)-chi2minusten*(zplusten*zplusten-zmin*zmin))/(2*(chi2plusten*(zminusten-zmin)+(zmin-zplusten)*chi2minusten));
	  a=(chi2plusten-chi2kmin)/(zminusten*zminusten-zmin*zmin-2*b*(zmin-zminusten));
	  c=chi2kmin-a*(zmin-b)*a*(zmin-b);/*
	  double tempA =(x2-x3)/(x2-x1);
	  c=(tempA*(y2-y1)-(y2-y3))/(tempA*(x2*x2-x1*x1)-(x2*x2-x3*x3));
	  b=((y2-y1)-c*(x2*x2-x1*x1))/(x2-x1);
	  a=y3-b*x3-c*x3*x3;*/
	  b=(x2*x2-x1*x1-(y2-y1)/(y3-y1)*(x3*x3-x1*x1))/(2*(y2-y1)/(y3-y1)*(x1-x3)-2*x1+2*x2);
	  a=(y2-y1)/(pow(x2-b,2)-pow(x1-b,2));
	 // std::cout <<"senitycheck"<<a-(y3-y1)/(pow(x3-b,2)-pow(x3-b,2))<< std::endl;
	  c=y1-a*pow(x1-b,2);
	  //deltaz=(chi2kmin+1-c)/a+b-zmin;
	  //deltaz=sqrt(pow(x1-b,2)+1/a)+b-x1 ;
	  deltaz=sqrt(1/a);
	  deltazkvskddivsigma=(zmin-frac1oldd)/deltaz;
	  tau1k->changeEandSave(etemp2);
	  invmk->getChi2();

  }

if(i<25){
double graphlowerlimit;
double graphupperlimit;
TGraph* gr=new TGraph();
gr->SetName(Form("graph_%d",i));
TGraph* aprox=new TGraph();
aprox->SetName(Form("aprox_%d",i));
TGraph* chidyn=new TGraph();
chidyn->SetName(Form("chidyn_%d",i));
TGraph* chikin=new TGraph();
chikin->SetName(Form("chikin_%d",i));
if((zmin-2*deltaz)<tau1k->getInitial4Vector().E()/tau1k->getUpperFitLimitE()){
   graphlowerlimit=tau1k->getInitial4Vector().E()/tau1k->getUpperFitLimitE();
}
else{
	graphlowerlimit=zmin-2*deltaz;
}
if((zmin+2*deltaz)>tau1k->getInitial4Vector().E()/tau1k->getLowerFitLimitE()){
   graphupperlimit=tau1k->getInitial4Vector().E()/tau1k->getLowerFitLimitE();
}
else{
	graphupperlimit=zmin+2*deltaz;
}
double stepsize=(graphupperlimit-graphlowerlimit)/100;
double etemp3=tau1k->getFit4Vector().E();
std::cout << "upper/lower z limit: " << graphupperlimit << "/ " << graphlowerlimit <<std::endl;
std::cout << "visible tau energy" <<tau1k->getInitial4Vector().E() << "upper/lower tau limit" <<  tau1k->getUpperFitLimitE() << "/" << tau1k->getLowerFitLimitE() << std::endl;
std::cout << "a/b/c" << a <<"/"<<b<<"/"<<c<<std::endl;
std::cout << "delta z" <<deltaz << std::endl;

for(unsigned int k=0; k<98;k++){
	tau1k->changeEandSave(tau1k->getInitial4Vector().E()/(graphlowerlimit+(k+1)*stepsize));
	invmk->getChi2();
	gr->SetPoint(k,tau1k->getInitial4Vector().E()/tau1k->getFit4Vector().E(),balancek->getChi2());
	aprox->SetPoint(k,tau1k->getInitial4Vector().E()/tau1k->getFit4Vector().E(),a*pow((tau1k->getInitial4Vector().E()/tau1k->getFit4Vector().E()-b),2)+c);
	tau1k->changeEandSave(etemp3);
	invmk->getChi2();

}

if(tau1kd->getInitial4Vector().E()/(zmin-deltaz)<tau1kd->getUpperFitLimitE()&&tau1kd->getInitial4Vector().E()/(zmin+deltaz)>tau1kd->getLowerFitLimitE()){
double etempkdbla=tau1kd->getFit4Vector().E();
for(unsigned int k=0; k<99;k++){
	double zchange=zmin-deltaz+0.02*k*deltaz;
	tau1kd->changeEandSave(tau1kd->getInitial4Vector().E()/zchange);
	invmkd->getChi2();
	chidyn->SetPoint(k,(tau1kd->getInitial4Vector().E()/tau1kd->getFit4Vector().E()-zmin)/deltaz,Likelihoodkd->getChi2());
	chikin->SetPoint(k,(tau1kd->getInitial4Vector().E()/tau1kd->getFit4Vector().E()-zmin)/deltaz,balancekd->getChi2());
	tau1kd->changeEandSave(etempkdbla);
	invmkd->getChi2();
}
}
else{
	for(unsigned int k=0; k<99;k++){
		chidyn->SetPoint(k,0,0);
	}
}

TGraph* onesigma=new TGraph();
//onesigma->SetPoint(1,zmin-deltaz,0);
onesigma->SetPoint(1,b+deltaz,0.9*chi2kmin);  // use b instead of zmin (minimum of fitkurve)
onesigma->SetPoint(2,b+deltaz,chi2kmin+3);
onesigma->SetPoint(3,b-deltaz,chi2kmin+3);
onesigma->SetPoint(4,b-deltaz,0.9*chi2kmin);
onesigma->SetFillColor(kGreen);
//onesigma->SetFillStyle(1);
onesigma->SetLineColor(kGreen);
TCanvas* c = new TCanvas(Form("c_%d",i),Form("c_%d",i),1500,1500);
c->GetPad(0)->SetBottomMargin(0.15);
TLine* zkdline=new TLine();
TLine* zkline=new TLine();
TLegend* leg =new TLegend(0,0,0.35,0.3);
leg->AddEntry(gr,"#chi_{k}^{2}(z)","l");
leg->AddEntry(aprox,"#chi_{k}^{2}(z)-fit","l");
leg->AddEntry(onesigma,"z_{min}#pm 1#sigma_{z}","F");
leg->AddEntry(zkdline,"Position of #chi^{2}-minimum in k-Fit","l");
leg->AddEntry(zkdline,"Position of #chi^{2}-minimum in kd-Fit","l");
gr->GetYaxis()->SetRangeUser(0.9*chi2kmin,chi2kmin+3);
gr->GetXaxis()->SetTitle("z=#frac{E_{Vis}}{E_{Fit}}");
gr->GetXaxis()->SetTitleOffset(1.5);
gr->GetYaxis()->SetTitle("#chi^{2}");
gr->SetTitle("Impact of the dynamic constraint on the final energyfraction");
gr->Draw("al");
onesigma->Draw("F");
gr->Draw("l");
aprox->SetLineColor(kBlue);
aprox->Draw("l");
zkdline->Draw();
zkdline->SetLineStyle(3);
zkdline->SetLineColor(kRed);
zkdline->DrawLine(frac1oldd,0.9*chi2kmin,frac1oldd,chi2kmin+3);
zkline->Draw();
zkline->SetLineStyle(0);
zkline->SetLineColor(kBlack);
zkline->DrawLine(zmin,0.9*chi2kmin,zmin,chi2kmin+3);
leg->Draw();

/*c->Print(Form("chi2funktion_%d.pdf[",i));
c->Print(Form("chi2funktion_%d.pdf",i));
c->Print(Form("chi2funktion_%d.pdf]",i));*/
if(i==0)c->Print("chi2funktion.pdf[");
c->Print("chi2funktion.pdf");
if(i==24)c->Print("chi2funktion.pdf]");

c->Close();
c->Delete();
gr->Delete();
aprox->Delete();
onesigma->Delete();
zkdline->Delete();
zkline->Delete();


TCanvas* c2 = new TCanvas(Form("c2_%d",i),Form("c2_%d",i),1500,1500);
c2 ->SetBottomMargin(0.18);
chikin->SetLineColor(kRed);
chikin->GetYaxis()->SetRangeUser(-5,5);
chikin->SetTitle("the #chi^{2}-contributions from the k(red) and d(black) constraint in the kd-fit");
chikin->GetXaxis()->SetTitle("#frac{z-z_{min,kfit}}{#sigma_{z}}");
chikin->GetXaxis()->SetTitleOffset(1.25);
chikin->GetYaxis()->SetTitle("#chi^{2}(z)");
chikin->Draw("al");
chidyn->Draw("l");
TLine* zkdline2=new TLine();
TLine* zkline2=new TLine();
zkdline2->Draw();
zkdline2->SetLineStyle(3);
zkline2->SetLineColor(kRed);
zkdline2->DrawLine((frac1oldd-zmin)/deltaz,-5,(frac1oldd-zmin)/deltaz,5);
zkline2->Draw();
zkline2->SetLineStyle(0);
zkline2->SetLineColor(kRed);
zkline2->DrawLine(0,-5,0,5);
TLegend *leg2=new TLegend(0.015,0.015,0.255,0.135);
leg2->AddEntry(zkdline2,"Position of #chi^{2}-minimum in k-Fit","l");
leg2->AddEntry(zkdline2,"Position of #chi^{2}-minimum in kd-Fit","l");
leg2->Draw();
if(i==0)c2->Print("dynamicconstraintoverzminpmsigmaz.pdf[");
c2->Print("dynamicconstraintoverzminpmsigmaz.pdf");
if(i==24)c2->Print("dynamicconstraintoverzminpmsigmaz.pdf]");
c2->Close();
chikin->Delete();
chidyn->Delete();
zkdline2->Delete();
zkline2->Delete();
c2->Delete();



}







  HHTauTau_tree->Fill();

  //clear heap

  delete(tau1k);
  delete(tau2k);
  delete(tau1kd);
  delete(tau2kd);
  delete(metk);
  delete(metkd);
  delete(higgsk);
  delete(higgskd);
  delete(invmk);
  delete(invmkd);
  delete(balancek);
  delete(balancekd);
  delete(Likelihoodkd);
  delete(singlefitk);
  delete(singlefitkd);

  }

  std::cout << "finished analysing, now writing outfile" <<std::endl;
  HHTauTau_tree->Write();


  delete(outfile);
  delete(HHTauTau_tree);
  delete(PDF1);
  delete(PDF2);










}
