// #include <iostream>

void third_ex(){

TCanvas* c1 = new TCanvas();

double a = 1.0;
double b = 10.0;
int n = 12;
Double_t x = 1.0;

Double_t f1a[n]; Double_t f1b[n]; Double_t f1c[n];
Double_t f2a[n]; Double_t f2b[n];
Double_t h[n];

  for (int i = 0; i < n; i++) {
		h[i] = a + i * (b - a) / (n - 1);
		f1a[i] = -exp(-x) - (exp(-(x + h[i])) - exp(-x)) / h[i];
		f1b[i] = -exp(-x) - (exp(-(x + h[i])) - exp(-(x - h[i]))) / (2 * h[i]);
		f1c[i] = -exp(-x) - (-3 * exp(-(x)) + 4 * exp(-(x + h[i])) - exp(-(x + 2 * h[i]))) / (2 * h[i]);
		f2a[i] = exp(-x) - (exp(-(x)) - 2 * exp(-(x + h[i])) + exp(-(x + 2 * h[i]))) / pow(h[i], 2);
		f2b[i] = exp(-x) - (exp(-(x + h[i])) - 2 * exp(-x) + exp(-(x - h[i]))) / pow(h[i], 2);
		f1a[i] = abs(f1a[i]);
		f1b[i] = abs(f1b[i]);
		f1c[i] = abs(f1c[i]);
		f2a[i] = abs(f2a[i]);
		f2b[i] = abs(f2b[i]);
		cout << f1c[i] << endl;
}

TGraph* gr1 = new TGraph(n, h, f1a);TGraph* gr2 = new TGraph(n, h, f1b);
TGraph* gr3 = new TGraph(n, h, f1c);TGraph* gr4 = new TGraph(n, h, f2a);
TGraph* gr5 = new TGraph(n, h, f2b);

TMultiGraph* mg = new TMultiGraph(); 
c1->cd();
mg->Add(gr1);mg->Add(gr2);mg->Add(gr3);mg->Add(gr4);mg->Add(gr5);
gr1->SetLineColor(kRed);
gr2->SetLineColor(kBlue);
gr3->SetLineColor(kGreen);
gr5->SetLineColor(kOrange);
gr1->SetMarkerStyle(24);
gr2->SetMarkerStyle(24);
gr3->SetMarkerStyle(24);
gr4->SetMarkerStyle(24);
gr5->SetMarkerStyle(24);
mg->Draw("ACP");

auto legend = new TLegend(0.6, 0.7, 0.8, 0.9);
legend->SetHeader("Legend");
legend->AddEntry(gr1, "Error for 1a");
legend->AddEntry(gr2, "Error for 1b");
legend->AddEntry(gr3, "Error for 1c");
legend->AddEntry(gr4, "Error for 2a");
legend->AddEntry(gr5, "Error for 2b");
legend->Draw();

}
