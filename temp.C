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


using namespace HHKinFit2;

int main(int argc, char* argv[])
{
  //Define Tree
  TTree * HHTauTau_tree=new TTree("HHTauTau_tree","HHTauTau_tree");
  HHTauTau_tree->SetDirectory(0);
  HHTauTau_tree->SetAutoSave(100000000);

  unsigned int j;
  HHTauTau_tree->Branch("eventnumber",&j);
  double etau1truth;
  HHTauTau_tree->Branch("Energy fromTau1 (truth)",&etau1truth);
  double etau2truth;
  HHTauTau_tree->Branch("Energy fromTau2 (truth)",&etau2truth);
  double etau1kfit;
  HHTauTau_tree->Branch("Energy fromTau1 (kfit)",&etau1kfit);
  double etau2kfit;
  HHTauTau_tree->Branch("Energy fromTau2 (kfit)",&etau2kfit);
  double etau1kdfit;
  HHTauTau_tree->Branch("Energy fromTau1 (kdfit)",&etau1kdfit);
  double etau2kdfit;
  HHTauTau_tree->Branch("Energy fromTau2 (kdfit)",&etau2kdfit);
  double etau1vis;
  HHTauTau_tree->Branch("visible Energy fromTau1 (vis)",&etau1vis);
  double etau2vis;
  HHTauTau_tree->Branch("visible Energy fromTau2 (vis)",&etau2vis);
  double phitau1;
  HHTauTau_tree->Branch("The angle phi from tau1",&phitau1);
  double phitau2;
  HHTauTau_tree->Branch("The angle phi from tau2",&phitau2);
  double etatau1;
  HHTauTau_tree->Branch("Eta from tau1",&etatau1);
  double etatau2;
  HHTauTau_tree->Branch("Eta phi from tau2",&etatau2);
  double METx;
  HHTauTau_tree->Branch("Missing energy in x-direction (METx)",&METx);
  double METy;
  HHTauTau_tree->Branch("Missing energy in y-direction (METy)",&METy);
  double chi2k;
  HHTauTau_tree->Branch("final chi^2 (kfit)",&chi2k);
  double chi2kd;
  HHTauTau_tree->Branch("final chi^2 (kdfit)",&chi2kd);
  double chi2kinematick;
  HHTauTau_tree->Branch("kinematic part of the final chi^2 (kfit)",&chi2kinematick);
  double chi2kinematickd;
  HHTauTau_tree->Branch("kinematic part of the final chi^2 (kdfit)",&chi2kinematickd);
  double chi2dynamick;
  HHTauTau_tree->Branch("dynamic part of the final chi^2 (kfit)",&chi2dynamick);
  double chi2dynamickd;
  HHTauTau_tree->Branch("dynamic part of the final chi^2 (kdfit)",&chi2dynamickd);


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

  for(unsigned int i=0; i<100000; i++){

    try{
	    Higgsgenerator.generateEvent();
	   }
	catch(const HHEnergyRangeException& e){
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
	singlefitk->addFitObjectE(tau1k);
    singlefitk->addConstraint(invmk);
    singlefitk->addConstraint(balancek);

    //KD-Fit -Fit

    HHKinFit* singlefitkd = new HHKinFit();
    singlefitkd->addFitObjectE(tau1kd);
    singlefitkd->addConstraint(invmkd);
    singlefitkd->addConstraint(balancekd);
    singlefitkd->addConstraint(Likelihoodkd);

    //K-Fit -Fitting

    try {
  		 singlefitk->fit();
  		 if (!((singlefitk->getConvergence()==1)||(singlefitk->getConvergence()==2))) {
  		     i--;
  		     continue;
  		    }
  	   }
  	catch(HHEnergyRangeException const& e){
  		  i--;
  		  continue;
  	     }

  	//KD-Fit -Fitting

    try {
  	  	 singlefitkd->fit();
  	     if (!((singlefitkd->getConvergence()==1)||(singlefitkd->getConvergence()==2))) {
  	  		 i--;
  	  		 continue;
  	  		}
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

  TFile* outfile=new TFile("newToyMctest.root","RECREATE");
  HHTauTau_tree->Write();













}
