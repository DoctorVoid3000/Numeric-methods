#include <math.h>

void Fillh(double*, const int);
void FillR(const int, double, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*, double*);

void Printh(double*, const int);

#define LOG 

int MainLab3_1(void){

	const int NSTEPS = 10;
	double x = 1.0;

	double h[NSTEPS];
	Fillh(h, NSTEPS);
	Printh(h, NSTEPS);

	double R1a[NSTEPS], R1b[NSTEPS], R1c[NSTEPS], R2a[NSTEPS], R2b[NSTEPS];
	double R1aTeor[NSTEPS], R1bTeor[NSTEPS], R1cTeor[NSTEPS], R2aTeor[NSTEPS], R2bTeor[NSTEPS];

	FillR(NSTEPS, x, h, R1a, R1b, R1c, R2a, R2b, R1aTeor, R1bTeor, R1cTeor, R2aTeor, R2bTeor);

	//***********Оформление**********//
	TCanvas *c1 = new TCanvas("c1", "c1");
	c1->SetTitle("Numeric Derive");
	c1->SetGrid();

	#ifdef LOG
	c1->SetLogy();
	c1->SetLogx();
	#endif

	TGraph *grR1a = new TGraph(NSTEPS, h, R1a);
	grR1a->SetTitle("Numeric Derive");
	grR1a->GetXaxis()->SetTitle("h");
	grR1a->GetYaxis()->SetTitle("R");
	grR1a->SetLineWidth(4);
	grR1a->SetMarkerColor(kGreen);
	grR1a->SetMarkerSize(3);
	grR1a->SetMarkerStyle(45);
	grR1a->SetLineColor(kGreen);
	
	TGraph *grR1b = new TGraph(NSTEPS, h, R1b);
	grR1b->SetLineWidth(4);
	grR1b->SetMarkerColor(kBlue);
	grR1b->SetMarkerSize(3);
	grR1b->SetMarkerStyle(45);
	grR1b->SetLineColor(kBlue);

	TGraph *grR1c = new TGraph(NSTEPS, h, R1c);
	grR1c->SetLineWidth(4);
	grR1c->SetMarkerColor(kRed);
	grR1c->SetMarkerSize(3);
	grR1c->SetMarkerStyle(45);
	grR1c->SetLineColor(kRed);
	//grR1b->GetYaxis()->SetRangeUser(0, 10);

	TGraph *grR2a = new TGraph(NSTEPS, h, R2a);
	grR2a->SetLineWidth(4);
	grR2a->SetMarkerColor(kOrange);
	grR2a->SetMarkerSize(3);
	grR2a->SetMarkerStyle(45);
	grR2a->SetLineColor(kOrange);

	TGraph *grR2b = new TGraph(NSTEPS, h, R2b);
	grR2b->SetLineWidth(4);
	grR2b->SetMarkerColor(kBlack);
	grR2b->SetMarkerSize(3);
	grR2b->SetMarkerStyle(45);
	grR2b->SetLineColor(kBlack);




	TGraph *grR1aTeor = new TGraph(NSTEPS, h, R1aTeor);
	grR1aTeor->SetLineWidth(6);
	grR1aTeor->SetLineStyle(9);
	grR1aTeor->SetMarkerColor(kGreen);
	grR1aTeor->SetMarkerSize(3);
	grR1aTeor->SetMarkerStyle(45);
	grR1aTeor->SetLineColor(kGreen);
	
	TGraph *grR1bTeor = new TGraph(NSTEPS, h, R1bTeor);
	grR1bTeor->SetLineWidth(6);
	grR1bTeor->SetLineStyle(9);
	grR1bTeor->SetMarkerColor(kBlue);
	grR1bTeor->SetMarkerSize(3);
	grR1bTeor->SetMarkerStyle(45);
	grR1bTeor->SetLineColor(kBlue);

	TGraph *grR1cTeor = new TGraph(NSTEPS, h, R1cTeor);
	grR1cTeor->SetLineWidth(6);
	grR1cTeor->SetLineStyle(9);
	grR1cTeor->SetMarkerColor(kRed);
	grR1cTeor->SetMarkerSize(3);
	grR1cTeor->SetMarkerStyle(45);
	grR1cTeor->SetLineColor(kRed);

	TGraph *grR2aTeor = new TGraph(NSTEPS, h, R2aTeor);
	grR2aTeor->SetLineWidth(6);
	grR2aTeor->SetLineStyle(9);
	grR2aTeor->SetMarkerColor(kOrange);
	grR2aTeor->SetMarkerSize(3);
	grR2aTeor->SetMarkerStyle(45);
	grR2aTeor->SetLineColor(kOrange);

	TGraph *grR2bTeor = new TGraph(NSTEPS, h, R2bTeor);
	grR2bTeor->SetLineWidth(6);
	grR2bTeor->SetLineStyle(9);
	grR2bTeor->SetMarkerColor(kBlack);
	grR2bTeor->SetMarkerSize(3);
	grR2bTeor->SetMarkerStyle(45);
	grR2bTeor->SetLineColor(kBlack);

	TMultiGraph* mg = new TMultiGraph(); 
	mg->Add(grR1a);
	mg->Add(grR1b);
	mg->Add(grR1c);
	mg->Add(grR2a);
	mg->Add(grR2b);
	mg->Add(grR1aTeor);
	mg->Add(grR1bTeor);
	mg->Add(grR1cTeor);
	mg->Add(grR2aTeor);
	mg->Add(grR2bTeor);
	mg->Draw("ACP");

	TLegend *legend = new TLegend(0.50, 0.1, 0.9, 0.40);
	legend->SetHeader("Numeric Derive","C"); 
    legend->AddEntry(grR1a, "R^{(1a)}", "lp");
    legend->AddEntry(grR1b, "R^{(1b)}", "lp");
    legend->AddEntry(grR1c, "R^{(1c)}", "lp");
    legend->AddEntry(grR2a, "R^{(2a)}", "lp");
    legend->AddEntry(grR2b, "R^{(2b)}", "lp");

    legend->AddEntry(grR1aTeor, "R^{(1a)}Teor", "lp");
    legend->AddEntry(grR1bTeor, "R^{(1b)}Teor", "lp");
    legend->AddEntry(grR1cTeor, "R^{(1c)}Teor", "lp");
    legend->AddEntry(grR2aTeor, "R^{(2a)}Teor", "lp");
    legend->AddEntry(grR2bTeor, "R^{(2b)}Teor", "lp");
    legend->Draw();

    // TCanvas *c2 = new TCanvas("c2", "c2");
    // c2->SetGrid();
    // c2->cd();

    // #ifdef LOG
	// c2->SetLogy();
	// #endif

    // TMultiGraph* mg2 = new TMultiGraph(); 
    // mg2->Add(grR1aTeor);
	// mg2->Draw("ACP");

	// TLegend *legend2 = new TLegend(0.50, 0.1, 0.9, 0.4);
	// legend2->SetHeader("Numeric Derive","C"); 
    // legend2->AddEntry(grR1aTeor, "R^{(2a)Teor}", "lp");
    // legend2->Draw();

    c1->SaveAs("NumericDerive.pdf");
    // c2->SaveAs("NumericDerive2.pdf");


	return 0;
}

//***********Функции заполнения**********//

void Fillh(double* h, const int NSTEPS)
{
	for (int i = 0; i < NSTEPS; ++i)
	{
		h[i] = (i+1)/10000.0;
	}
}

void FillR(const int NSTEPS, double x, double* h, double* R1a, double* R1b, double* R1c, double* R2a, double* R2b, double* R1aTeor, double* R1bTeor, double* R1cTeor, double* R2aTeor, double* R2bTeor)
{
	for (int i = 0; i < NSTEPS; ++i)
	{
		#ifndef LOG
		R1a[i] = -exp(-x) - ( exp(-(x+h[i])) - exp(-x) )/h[i];
		R1b[i] = -exp(-x) - ( exp(-(x+h[i])) - exp(-(x-h[i])) )/(2*h[i]);
		R1c[i] = -exp(-x) - ( -3*exp(-x) + 4*exp(-(x+h[i])) - exp(-(x+2*h[i])) )/(2*h[i]);
		R2a[i] = exp(-x) - ( exp(-x) - 2*exp(-(x+h[i])) + exp(-(x+2*h[i])) )/(h[i]*h[i]);
		R2b[i] = exp(-x) - ( exp(-(x+h[i])) -2*exp(-x) + exp(-(x-h[i])) )/(h[i]*h[i]);
		cout << R2a[i] << "\t" << R2b[i] << endl;
		#endif

		#ifdef LOG
		R1a[i] = fabs(-exp(-x) - ( exp(-(x+h[i])) - exp(-x) )/h[i]);
		R1b[i] = fabs(-exp(-x) - ( exp(-(x+h[i])) - exp(-(x-h[i])) )/(2*h[i]));	
		R1c[i] = fabs(-exp(-x) - ( -3*exp(-x) + 4*exp(-(x+h[i])) - exp(-(x+2*h[i])) )/(2*h[i]));
		R2a[i] = fabs(exp(-x) - ( exp(-x) - 2*exp(-(x+h[i])) + exp(-(x+2*h[i])) )/(h[i]*h[i]));
		R2b[i] = fabs(exp(-x) - ( exp(-(x+h[i])) -2*exp(-x) + exp(-(x-h[i])) )/(h[i]*h[i]));

		//teor
		R1aTeor[i] = fabs(-0.5*exp(-x)*h[i]);
		R1bTeor[i] = fabs(-(1.0/6.0)*(-exp(-x))*h[i]*h[i]);
		R1cTeor[i] = fabs((-1.0/3.0)*(-exp(-(x+h[i])))*h[i]*h[i]);
		R2aTeor[i] = fabs(-exp(-x)*h[i]);
		R2bTeor[i] = fabs(-(1.0/12.0)*exp(-x)*h[i]*h[i]);

		cout << R1a[i] << "\t" << R1aTeor[i] << endl;
		#endif
	}
}

//***********Конец функций заполнения**********//



//***********Функции печати**********//

void Printh(double* h, const int NSTEPS)
{
	for (int i = 0; i < NSTEPS; ++i)
	{
		printf("h[%d] = %f\n", i, h[i]);
	}
	printf("\n");
}
//***********Конец функций печати**********//