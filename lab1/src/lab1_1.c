#include <stdio.h>
#include <math.h>

#define NUM_POINTS  7
#define NUM_POLINOMS_POINTS  500

#define LEFT_BOARD 0.0
#define RIGHT_BOARD 3.0
#define LEFT_LIMIT -0.5
#define RIGHT_LIMIT 3.5

#define POLINOM_LEFT_LIMIT -10.0
#define POLINOM_RIGHT_LIMIT 10.0


void FillX(double*);
void FillY(double*, double*);
void FillChebX(double*);
void FillChebY(double*, double*);
void FillMatrixF(TMatrixD&, double*);
void FillVectorY(TVectorD&, double*);
void FillPolinomsPoints(double*);

TF1* PolinomMethod(double*, double*);
TF1* PolinomMethodCheb(double*, double*);

void PrintX(double*);
void PrintY(double*);
void PrintMatrixF(TMatrixD);
void PrintVectorY(TVectorD);
void PrintVectorA(TVectorD);
//void PrintPolinomsPoints(double*);

int Factorial(int);

int MainLab2()
{
	double x[NUM_POINTS];
	double x_cheb[NUM_POINTS];

	double y[NUM_POINTS];
	double y_cheb[NUM_POINTS];

	FillX(x);
	FillChebX(x_cheb);

	FillY(y, x);
	FillChebY(y_cheb, x_cheb);
	
	PrintX(x_cheb);

	double PolinomsPoints[NUM_POLINOMS_POINTS];
	FillPolinomsPoints(PolinomsPoints);

	TF1 *PolinomFunc = PolinomMethod(x, y);
	TF1 *PolinomFuncCheb = PolinomMethodCheb(x_cheb, y_cheb);

	TCanvas *canv = new TCanvas("canv", "canv");
	canv->SetTitle("Chebyshev");
	canv->SetGrid();
	canv->Divide(1, 2);

	canv->cd(1);
    gPad->SetGrid();

	TF1 *OriginalFunc = new TF1("OriginalFunc", "sin(x)*exp(-x)", LEFT_LIMIT, RIGHT_LIMIT);
	TGraph *OriginalFuncPoints = new TGraph(NUM_POINTS, x, y);
	TGraph *OriginalFuncPointsCheb = new TGraph(NUM_POINTS, x_cheb, y_cheb);

	OriginalFunc->SetTitle("Result of interpolation of function y(x) = sin(x) #upoint e^{-x}");
	OriginalFunc->GetXaxis()->SetTitle("X");
	OriginalFunc->GetYaxis()->SetTitle("Y");
	OriginalFunc->SetLineWidth(10);
	OriginalFunc->Draw();
	OriginalFunc->GetYaxis()->SetRangeUser(-0.5, 0.4);

	OriginalFuncPoints->SetMarkerColor(kBlue);
	OriginalFuncPoints->SetMarkerSize(3);
	OriginalFuncPoints->SetMarkerStyle(45);
	OriginalFuncPoints->Draw("P same");

	OriginalFuncPointsCheb->SetMarkerColor(kTeal);
	OriginalFuncPointsCheb->SetMarkerSize(3);
	OriginalFuncPointsCheb->SetMarkerStyle(45);
	OriginalFuncPointsCheb->Draw("P same");

	PolinomFunc->SetLineColor(kBlue);
	PolinomFunc->SetLineStyle(1);
	PolinomFunc->SetLineWidth(6);
	PolinomFunc->Draw("same");

	PolinomFuncCheb->SetLineColor(kTeal);
	PolinomFuncCheb->SetLineStyle(1);
	PolinomFuncCheb->SetLineWidth(2);
	PolinomFuncCheb->Draw("same");

	TLegend *legend = new TLegend(0.50, 0.1, 0.9, 0.45);
    legend->SetHeader("Interpolation Functions ","C"); 
    legend->AddEntry(OriginalFunc,"OriginalFunc y(x) = sin(x) #upoint e^{-x}", "lp");
    legend->AddEntry(PolinomFunc,"PolinomFunc #Phi_{n}(x) = #sum_{j = 0}^{n} a_{j} #upoint x^{j}, const", "lp");
    legend->AddEntry(PolinomFuncCheb,"PolinomFunc #Phi_{n}(x) = #sum_{j = 0}^{n} a_{j} #upoint x^{j}, Chebyshev", "lp");
    legend->Draw();


    TF1 *UncertanityFunc = new TF1("Uncetarnities", "abs(OriginalFunc - PolinomFunc)", LEFT_LIMIT, RIGHT_LIMIT);
    TF1 *UncertanityFuncCheb = new TF1("UncetarnitiesCheb", "abs(OriginalFunc - PolinomFuncCheb)", LEFT_LIMIT, RIGHT_LIMIT);



    canv->cd(2);
    gPad->SetGrid();

    TLegend *legend2 = new TLegend(0.30, 0.55, 0.70, 0.90);

    UncertanityFunc->SetTitle("R(x) = |y(x) - #Phi(x)|");
    UncertanityFunc->SetLineColor(kBlue);
	UncertanityFunc->SetLineStyle(1);
	UncertanityFunc->SetLineWidth(2);
	UncertanityFunc->GetYaxis()->SetRangeUser(0.0, 0.0009);
    UncertanityFunc->Draw();

    UncertanityFuncCheb->SetLineColor(kTeal);
	UncertanityFuncCheb->SetLineStyle(1);
	UncertanityFuncCheb->SetLineWidth(2);
    UncertanityFuncCheb->Draw("same");

    legend2->AddEntry(UncertanityFunc, Form("UncertanityFunc R(x) = |y(x) - #Phi(x)|, n = %d, const", NUM_POINTS), "lp");
    legend2->AddEntry(UncertanityFuncCheb, Form("UncertanityFuncCheb R(x) = |y(x) - #Phi(x)|, n = %d, Chebyshev", NUM_POINTS), "lp");

    legend2->Draw();

    int k = NUM_POINTS;

    double omega_NUM_POINTS[NUM_POLINOMS_POINTS];
    double R_max_NUM_POINTS[NUM_POLINOMS_POINTS];

    for (int i = 0; i < NUM_POLINOMS_POINTS; ++i)
    {
    	omega_NUM_POINTS[i] = 1;
    	for (int j = 0; j < NUM_POINTS; ++j)
    	{
    		omega_NUM_POINTS[i] *= (PolinomsPoints[i] - x_cheb[j]);
    	}
    	R_max_NUM_POINTS[i] = 16*sin(M_PI/4)*exp(-M_PI/4)*fabs(omega_NUM_POINTS[i])/Factorial(k);
    }

    TGraph *R_MAX = new TGraph(NUM_POLINOMS_POINTS, PolinomsPoints, R_max_NUM_POINTS);

    R_MAX->SetLineColor(kTeal);
	R_MAX->SetLineStyle(9);
	R_MAX->SetLineWidth(3);
	//R_MAX->Draw("same C");

	int k2 = NUM_POINTS;

    double omega_NUM_POINTS2[NUM_POLINOMS_POINTS];
    double R_max_NUM_POINTS2[NUM_POLINOMS_POINTS];

    for (int i = 0; i < NUM_POLINOMS_POINTS; ++i)
    {
    	omega_NUM_POINTS2[i] = 1;
    	for (int j = 0; j < NUM_POINTS; ++j)
    	{
    		omega_NUM_POINTS2[i] *= (PolinomsPoints[i] - x[j]);
    	}
    	R_max_NUM_POINTS2[i] = 16*sin(M_PI/4)*exp(-M_PI/4)*fabs(omega_NUM_POINTS2[i])/Factorial(k2);
    	//R_max_NUM_POINTS2[i] = -8*exp(-3*M_PI/4)*(sin(3*M_PI/4)+cos(3*M_PI/4))*fabs(omega_NUM_POINTS[i])/Factorial(k);
    }

    TGraph *R_MAX2 = new TGraph(NUM_POLINOMS_POINTS, PolinomsPoints, R_max_NUM_POINTS2);

    R_MAX2->SetLineColor(kBlue);
	R_MAX2->SetLineStyle(9);
	R_MAX2->SetLineWidth(3);
	//R_MAX2->Draw("same C");

	//legend2->AddEntry(R_MAX, Form("R(x) #leq #frac{max|y^{n+1}(x)|}{(n+1)!}#upoint|#omega_{n+1}(x)|, Chebyshev n = %d", NUM_POINTS), "lp");
	//legend2->AddEntry(R_MAX2, Form("R(x) #leq #frac{max|y^{n+1}(x)|}{(n+1)!}#upoint|#omega_{n+1}(x)|, function n = %d", NUM_POINTS), "lp");

	canv->SaveAs("Chebyshev.pdf");

	return 0;
}



void FillX(double* x)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		x[i] = LEFT_BOARD + (RIGHT_BOARD - LEFT_BOARD)/(NUM_POINTS - 1) * i;
	}
}


void FillY(double* y, double* x)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		y[i] = sin(x[i])*exp(-x[i]);
	}
}

void FillChebX(double* x_cheb)
{
	double t[NUM_POINTS];
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		t[i] = cos(1./NUM_POINTS*(M_PI/2. + M_PI*i));
		x_cheb[i] = (RIGHT_BOARD - LEFT_BOARD)/2.0*t[i] + (LEFT_BOARD + RIGHT_BOARD)/2.0;
	}
}


void FillChebY(double* y_cheb, double* x_cheb)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		y_cheb[i] = sin(x_cheb[i])*exp(-x_cheb[i]);
	}
}

void FillPolinomsPoints(double* PolinomsPoints)
{
	for (int i = 0; i < NUM_POLINOMS_POINTS; ++i)
	{

		PolinomsPoints[i] = POLINOM_LEFT_LIMIT + (POLINOM_RIGHT_LIMIT - POLINOM_LEFT_LIMIT)/(NUM_POLINOMS_POINTS - 1)*i;
	}
}

void FillMatrixF(TMatrixD& F, double* arr)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		for (int j = 0; j < NUM_POINTS; ++j)
		{
			F(i, j) = pow(arr[i], (double)j);
		}
	}
}

void FillVectorY(TVectorD& Y, double* arr)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		Y(i) = arr[i];
	}
}

//***********Конец функций заполнения**********//


TF1* PolinomMethod(double* x, double* y)
{
	TMatrixD F(NUM_POINTS, NUM_POINTS);
	FillMatrixF(F, x);

	//PrintMatrixF(F);

	if(F.Determinant() == 0)
	{
		printf("Determinant of matrix is 0\n");
		exit(1);
		return 0;
	}

	TVectorD Y(NUM_POINTS);
	FillVectorY(Y, y);

	TVectorD A = F.Invert() * Y;


	//PrintVectorA(A);

	string formula = "pol" + to_string(NUM_POINTS - 1);

	TF1 *PolinomFunc = new TF1("PolinomFunc", formula.c_str(), LEFT_LIMIT, RIGHT_LIMIT); 

	for (int i = 0; i < NUM_POINTS; ++i)
	{
		PolinomFunc->SetParameter(i, A(i));
	}

	return PolinomFunc;
}

TF1* PolinomMethodCheb(double* x_cheb, double* y_cheb)
{
	TMatrixD F(NUM_POINTS, NUM_POINTS);
	FillMatrixF(F, x_cheb);

	//PrintMatrixF(F);

	if(F.Determinant() == 0)
	{
		printf("Determinant of matrix is 0\n");
		exit(1);
		return 0;
	}

	TVectorD Y(NUM_POINTS);
	FillVectorY(Y, y_cheb);

	TVectorD A = F.Invert() * Y;


	//PrintVectorA(A);

	string formula = "pol" + to_string(NUM_POINTS - 1);

	TF1 *PolinomFuncCheb = new TF1("PolinomFuncCheb", formula.c_str(), LEFT_LIMIT, RIGHT_LIMIT); 

	for (int i = 0; i < NUM_POINTS; ++i)
	{
		PolinomFuncCheb->SetParameter(i, A(i));
	}

	return PolinomFuncCheb;
}

//***********Функции печати**********//

void PrintX(double* x)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		printf("x[%d] = %f\n", i, x[i]);
	}
	printf("\n");
}

void PrintY(double* y)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		printf("y[%d] = %f\n", i, y[i]);
	}
	printf("\n");
}

void PrintMatrixF(TMatrixD F)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		for (int j = 0; j < NUM_POINTS; ++j)
		{
			printf("A(%d, %d) = %f \t", i, j, F(i, j));
		}
		printf("\n");
	}
	printf("\n");
}

void PrintVectorY(TVectorD Y)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		printf("Y(%d) = %f\n", i, Y(i));
	}
	printf("\n");
}

void PrintVectorA(TVectorD A)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		printf("A(%d) = %f\n", i, A(i));
	}
	printf("\n");
}

// void PrintPolinomsPoints(double* PolinomsPoints)
// {
// 	printf("\n");
// 	for (int i = 0; i < NUM_POLINOMS_POINTS; ++i)
// 	{
// 		printf("PolinomsPoints[%d] = %f \n", i, PolinomsPoints[i]);
// 	}
// 	printf("\n");
// }

int Factorial(int k)
{
	int result = 1;
	for (int i = 1; i <= k; ++i)
	{
		result *= i;
	}

	return result;
}





