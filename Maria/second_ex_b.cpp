// #include <TMatrixD.h>
// #include <TVectorD.h>
// #include <iostream>
// #include <Tmath.h>
// #include "TROOT.h"
// #include "TLegend.h"

Double_t func(TVector A, TVector X, int i, Double_t x) {

	return (A(i) * (X(i + 1) - x) + A(i + 1) * (x - X(i))) / (X(i + 1) - X(i));

}

Double_t BB(TVector a, TVector X, Double_t h, Double_t x, int i) {
	cout << x << " " << X(i + 3) << endl;
	return (a(i) * pow((2 - (x - X(i)) / h), 3) / 6 + a(i + 1) * (0.6667 - pow(((x - X(i + 1)) / h), 2) + pow(((x - X(i + 1)) / h), 3) / 2) + a(i + 2) * (0.6667 - pow(((x - X(i + 2)) / h), 2) - pow(((x - X(i + 2)) / h), 3) / 2) + a(i + 3) * pow((2 + (x - X(i + 3)) / h), 3) / 6);
}

void seconda_ex() {

	TCanvas* c1 = new TCanvas();

	int n = 10;
	Double_t a1 = 0;
	Double_t b1 = 3;
	TVectorD A(n);
	TVectorD X(n);
	TVectorD X2(n+2);
	TMatrixD A2(n+2, n+2);
	TVectorD F(n + 2);
	TVectorD a(n + 2);
	Double_t xxx[n];
	Double_t yyy[n];
	for (int i = 0; i < n; i++) {
		X(i) = a1 + (b1 - a1) * i / (n - 1);
		X2(i+1) = X(i);
		A(i) = sin(X(i)) * exp(-X(i));
		xxx[i] = X(i);
		yyy[i] = A(i);
	}
	X2(0) = (a1 - b1) / (n - 1);
	X2(n + 1) = X2(n) + (b1 - a1) / (n - 1);
	

	int p = 10;
	int z = (n - 1) * p;
	Double_t x[z];
	Double_t y[z];
	Double_t er[z];

	for (int i = 0; i < n - 1; i++) {
		for (int j = p * i; j < p * (i + 1); j++) {
			x[j] = a1 + (b1 - a1) * j / (z - 1);
			y[j] = func(A, X, i, x[j]);
			er[j] = abs(y[j] - sin(x[j]) * exp(-x[j]));
		}
	}

	int k = 0;
	Double_t h = (b1 - a1) / (n - 1);
	Double_t s = 0.66667;

	for (int i = 0; i < n+2; i++) {

		if ((i != 0) && (i != (n + 1))) {

			F(i) = sin(X2(i))*exp(-X2(i));
			
			for (int j = 0; j < n+2; j++) {

				A2(i, j) = 0;

				if (j == k) {

					A2(i, j) = pow((2 - (X2(i) - X2(i - 1)) / h), 3) / 6;
					//cout << i << " " << j << " A2 " << A2(i, j) << " ";
				}
				if (j == k + 1) {
				
					A2(i, j) = s;
					//cout << i << " " << j << " A2 " << A2(i, j) << " ";
				}
				if (j == k + 2) {
					
					A2(i, j) = pow((2 + (X2(i) - X2(i + 1)) / h), 3) / 6;
					//cout << i << " " << j << " A2 " << A2(i, j) << " ";
				}
			}
			k++;
		}

		if ((i == 0) || (i == (n + 1))) {
			F(i) = 0;
			for (int j = 0; j < n+2; j++) {

				A2(i, j) = 0;
				//cout << i << " " << j << " A2 " << A2(i, j) << " ";

			}
		}
		//cout << " " << endl;
	}

	A2(0, 0) = (2 - (X2(1) - X2(0)) / h);
	A2(0, 1) = -2;
	A2(0, 2) = (2 + (X2(1) - X2(2)) / h);
	A2(n + 1, n - 1) = (2 - (X2(n) - X2(n-1)) / h);
	A2(n + 1, n) = -2;
	A2(n + 1, n + 1) = (2 + (X2(n) - X2(n+1)) / h);



	for (int i = 0; i < n + 2; i++) {
		for (int j = 0; j < n + 2; j++) {
			cout << A2(i, j) << " " ;
		}
		cout << " " << endl;
	}
	a = A2.Invert() * F;


	for (int i = 0; i < n + 2; i++) {
		cout << a(i) << " " << i << endl;
	}
	cout << "fuck" << endl;
	for (int i = 0; i < n + 2; i++) {
		cout << X2(i) << " " << i << endl;
	}
	int p2 = 10;
	int z2 = (n - 1) * p2;
	Double_t x2[z2];
	Double_t y2[z2];	
	Double_t er2[z2];	
	for (int i = 0; i < n - 1; i++) {
		for (int j = p2 * i; j < p2 * (i + 1); j++) {
			x2[j] = a1 + (b1 - a1) * j / (z2 - 1);
			y2[j] = BB(a, X2, h, x2[j], i);
			er2[j] = abs(y2[j] - sin(x[j]) * exp(-x[j]));
		}
	}

	Double_t xx[n];
	for (int i = 0; i < n; i++) {
		xx[i] = X2(i + 1);
		cout << "BB " << BB(a, X2, h, xx[i], 0) << endl;
	}

	TGraph* gr1 = new TGraph(z, x, y); 	TGraph* grer1 = new TGraph(z, x, er);
	TGraph* gr2 = new TGraph(z2, x2, y2); 	TGraph* grer2 = new TGraph(z2, x2, er2);
	TMultiGraph* mg = new TMultiGraph(); 	TMultiGraph* mger = new TMultiGraph();
	TGraph* grr = new TGraph(n, xxx, yyy);

	c1->Divide(1, 2);
	c1->cd(1);
	mg->Add(gr1);
	gr1->SetLineColor(kRed);
	mg->Add(gr2);
	mg->Add(grr);
	grr->SetLineWidth(0);
	grr->SetMarkerStyle(24);
	mg->Draw("ACP");

	auto legend1 = new TLegend(0.6, 0.7, 0.8, 0.9);
	legend1->SetHeader("Legend");
	legend1->AddEntry(gr1, "first-order B-spline");
	legend1->AddEntry(gr2, "third-order B-spline");
	legend1->Draw();


	c1->cd(2);
	mger->Add(grer1);
	mger->Add(grer2);
	grer1->SetLineColor(kRed);
	mger->Draw("ACP");
	auto legend12 = new TLegend(0.6, 0.7, 0.8, 0.9);
	legend12->SetHeader("Legend");
	legend12->AddEntry(grer1, "Error for first-order B-spline");
	legend12->AddEntry(grer2, "Error for third-order B-spline");
	legend12->Draw();
	cout << h << endl;
}