// #include <TMatrixD.h>
// #include <TVectorD.h>
// #include <iostream>
// #include <Tmath.h>
// #include "TROOT.h"
// #include "TLegend.h"


void second_ex() {

	TCanvas* c1 = new TCanvas();

	int n = 10;
	Double_t a1 = 0;
	Double_t b1 = 3;
	int k = 0;
	TMatrixD A(n, n);	TMatrixD A2(n, n);
	TVectorD B(n);
	TVectorD X(n);
	TVectorD Y(n);
	TVectorD C(n);	TVectorD C2(n);
	TVectorD F(n);
	TVectorD b(n-1);	TVectorD b2(n-1);
	TVectorD d(n-1);	TVectorD d2(n-1);
	TVectorD a(n-1);	TVectorD a2(n-1);
	Double_t grrx[n];
	Double_t grry[n];

	for (int i = 0; i < n; i++) {

		X(i) = (b1 - a1) * i / (n - 1);
		Y(i) = sin(X(i)) * exp(-X(i));
		grrx[i] = X(i);
		grry[i] = Y(i);
	}


	for (int i = 0; i < n; i++) {

		if ((i != 0) && (i != (n-1))) {
			int i1 = i + 1;
			int i2 = i - 1;
			F(i) = 6 * (Y(i1) - Y(i)) / (X(i1) - X(i)) - 6 * (Y(i) - Y(i2)) / (X(i) - X(i2));
			cout << "i= " << i << " "<< Y(i1) <<" "<<Y(i)<<" "<<Y(i2)<< endl;

			for (int j = 0; j < n; j++) {
				A(i, j) = 0;
				if (j == k) {
					int p = k + 1;
					A(i, j) = X(p) - X(j);
				}
				if (j == k + 1) {
					int p = k + 2;
					A(i, j) = 2 * (X(p) - X(k));
				}
				if (j == k + 2) {
					int p = k + 1;
					A(i, j) = X(j) - X(p);
				}
			}
			k++;
		}

		if ((i == 0) || (i == (n-1))) {
			F(i) = 0;
			for (int j = 0; j < n; j++) {

				A(i, j) = 0;

				if ((i == 0) && (j == 0)) {
					A(i, j) = 1;
				}
				if ((i == (n-1)) && (j == (n-1))) {
					A(i, j) = 1;
				}
			}
		}
	}
	A2 = A;

	A2(0, 0) = 1 / (X(1) - X(0));
	A2(0, 1) = -1 / (X(1) - X(0)) - 1 / (X(2) - X(1));
	A2(0, 2) = 1 / (X(2) - X(1));
	A2(n - 1, n - 1) = 1 / (X(n - 1) - X(n - 2));
	A2(n-1, n-2) = -1 / (X(n-2) - X(n-3)) - 1 / (X(n-1) - X(n-2));
	A2(n - 1, n - 3) = 1 / (X(n - 2) - X(n - 3));

	C2 = A2.Invert() * F;
	C = A.Invert() * F;

	for (int i = 0; i < n-1; i++) {
		int i1 = i + 1;
		a(i) = Y(i1);
		d(i) = (C(i1) - C(i)) / (X(i1) - X(i));
		b(i) = (Y(i1) - Y(i)) / (X(i1) - X(i)) + C(i1) * (X(i1) - X(i)) / 3 + C(i) * (X(i1) - X(i)) / 6;
		a2(i) = Y(i1);
		d2(i) = (C2(i1) - C2(i)) / (X(i1) - X(i));
		b2(i) = (Y(i1) - Y(i)) / (X(i1) - X(i)) + C2(i1) * (X(i1) - X(i)) / 3 + C2(i) * (X(i1) - X(i)) / 6;
	}
	int p = 10;
	int z = (n - 1) * p;
	Double_t x[z];	
	Double_t y[z];	Double_t y2[z];
	Double_t er[z];	Double_t er2[z];
	for (int i = 0; i < n - 1; i++) {
		int i1 = i + 1;
		for (int j = p * i; j < p * (i + 1); j++) {
			x[j] = a1 + (b1 - a1) * j / (z - 1);
			y[j] = a(i) + b(i) * (x[j] - X(i1)) + C(i1) * pow((x[j] - X(i1)), 2) / 2 + d(i) * pow((x[j] - X(i1)), 3) / 6;
			er[j] = abs(y[j] - sin(x[j]) * exp(-x[j]));
			y2[j] = a2(i) + b2(i) * (x[j] - X(i1)) + C2(i1) * pow((x[j] - X(i1)), 2) / 2 + d2(i) * pow((x[j] - X(i1)), 3) / 6;
			er2[j] = abs(y2[j] - sin(x[j]) * exp(-x[j]));
		}
	}
	TGraph* gr1 = new TGraph(z, x, y); 	TGraph* grer1 = new TGraph(z, x, er);	
	TGraph* gr2 = new TGraph(z, x, y2); 	TGraph* grer2 = new TGraph(z, x, er2);
	TMultiGraph* mg = new TMultiGraph();	TMultiGraph* mger = new TMultiGraph();
	TGraph* grr = new TGraph(n, grrx, grry);

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
	legend1->AddEntry(gr1, "Second derivative");
	legend1->AddEntry(gr2, "Third derivative");
	legend1->Draw();

	c1->cd(2);
	mger->Add(grer1);
	grer1->SetLineColor(kRed);
	mger->Add(grer2);
	mger->Draw("ACP");
	auto legend12 = new TLegend(0.6, 0.7, 0.8, 0.9);
	legend12->SetHeader("Legend");
	legend12->AddEntry(grer1, "Error for second derivative");
	legend12->AddEntry(grer2, "Error for third derivative");
	legend12->Draw();

}