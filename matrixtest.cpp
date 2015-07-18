#include "TMatrixD.h"
#include "TVectorD.h"
#include "TArrayD.h"

int main(){

	TVectorD testvector(2);
	testvector[0]=25;
	testvector[1]=92;

/*
	double m_tau1vis=2;
	double m_tau2vis=3;
	double* b=&m_tau2vis;
	double* a=&m_tau1vis;
	testvector.SetElements(a);
	testvector.SetElements(b);
*/
	testvector.Print();

	TMatrixD testmatrix(2,2);
	testmatrix[0][0]=1;
	testmatrix[0][1]=2;
	testmatrix[1][0]=3;
	testmatrix[1][1]=4;
/*

	TArrayD fill(4);
	fill[0]=1;
	fill[1]=2;
	fill[2]=3;
	fill[3]=4;
	testmatrix.SetMatrixArray(fill.GetArray());
*/
	testmatrix.Print();


	TVectorD testergebnis(2);
	testergebnis=testmatrix*testvector;
	testergebnis.Print();




	return(0);
}
