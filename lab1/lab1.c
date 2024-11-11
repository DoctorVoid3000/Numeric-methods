#include <stdio.h>
#include <math.h>

#define NUM_POINTS  7
#define NUM_POLINOMS_POINTS  1000

#define LEFT_BOARD 0.0
#define RIGHT_BOARD 3.0

#define LEFT_LIMIT -0.5
#define RIGHT_LIMIT 3.5

#define POLINOM_LEFT_LIMIT -10.0
#define POLINOM_RIGHT_LIMIT 10.0

#define NUM_POINTS_2 8

void FillX(double*);
void FillY(double*, double*);
void FillMatrixF(TMatrixD&, double*);
void FillVectorY(TVectorD&, double*);
void FillPolinomsPoints(double*);

TF1* PolinomMethod(double*, double*);
TGraph* LagrangeMethod(double*, double*, double*);
void DevDif(double*, double*, double*);
TGraph* NewtonMethod(double*, double*, double*, double*);

void PrintX(double*);
void PrintY(double*);
void PrintMatrixF(TMatrixD);
void PrintVectorY(TVectorD);
void PrintVectorA(TVectorD);
void PrintPolinomsPoints(double*);

void FillX2(double*);
void FillY2(double*, double*);
void FillMatrixF2(TMatrixD&, double*);
void FillVectorY2(TVectorD&, double*);
TF1* PolMethod2(double*, double*);
void PrintVectorA2(TVectorD);
void PrintX2(double*);
void PrintY2(double*);
void PrintMatrixF2(TMatrixD);
void PrintVectorY2(TVectorD Y);

int Factorial(int);	

int MainLab()
{
	double x[NUM_POINTS];
	double y[NUM_POINTS];

	FillX(x);
	FillY(y, x);

	//PrintY(y);

	//**********Метод1 - степенной полином**********//
	TF1 *PolinomFunc = PolinomMethod(x, y);
	//**********Конец метода1**********//

	double PolinomsPoints[NUM_POLINOMS_POINTS];
	FillPolinomsPoints(PolinomsPoints);

	//**********Метод2 - полином Лагранжа**********//
	TGraph *LagrangeGraph = LagrangeMethod(PolinomsPoints, x, y);
	//**********Конец метода2**********//

	//**********Метод3 - полином Ньютона**********//
	double Differences[NUM_POINTS];
	DevDif(Differences, x, y);
	TGraph *NewtonGraph = NewtonMethod(PolinomsPoints, Differences, x, y);
	//**********Конец метода3**********//

	//***********Оформление**********//
	TCanvas *canv = new TCanvas("Result of interpolation", "Result of interpolation");
	canv->SetTitle("Result of interpolation");
	canv->SetGrid();

	TF1 *OriginalFunc = new TF1("OriginalFunc", "sin(x)*exp(-x)", LEFT_LIMIT, RIGHT_LIMIT);
	TGraph *OriginalFuncPoints = new TGraph(NUM_POINTS, x, y);

	OriginalFunc->SetTitle("Result of interpolation of function y(x) = sin(x) #upoint e^{-x}");
	OriginalFunc->GetXaxis()->SetTitle("X");
	OriginalFunc->GetYaxis()->SetTitle("Y");
	OriginalFunc->SetLineWidth(10);
	OriginalFunc->Draw();
	//OriginalFunc->GetYaxis()->SetRangeUser(-0.5, 0.4);

	OriginalFuncPoints->SetMarkerColor(kViolet);
	OriginalFuncPoints->SetMarkerSize(4);
	OriginalFuncPoints->SetMarkerStyle(45);
	OriginalFuncPoints->Draw("P same");

	PolinomFunc->SetLineColor(kBlue);
	PolinomFunc->SetLineStyle(1);
	PolinomFunc->SetLineWidth(7);
	PolinomFunc->Draw("same");

	LagrangeGraph->SetLineColor(kBlack);
	LagrangeGraph->SetLineStyle(1);
	LagrangeGraph->SetLineWidth(4);
	LagrangeGraph->Draw("same C");

	NewtonGraph->SetLineColor(kGreen);
	NewtonGraph->SetLineStyle(1);
	NewtonGraph->SetLineWidth(1);
	NewtonGraph->Draw("same C");

	TLegend *legend = new TLegend(0.30, 0.1, 0.9, 0.60);
    legend->SetHeader("Interpolation Functions ","C"); 
    legend->AddEntry(OriginalFunc,"OriginalFunc y(x) = sin(x) #upoint e^{-x}", "lp");
    legend->AddEntry(PolinomFunc,"PolinomFunc #Phi_{n}(x) = #sum_{j = 0}^{n} a_{j} #upoint x^{j}", "lp");
    legend->AddEntry(LagrangeGraph,"LagrangeGraph L_{n}(x) = #sum_{j = 0}^{n} y_{j} #upoint l_{j}", "lp");
    legend->AddEntry(NewtonGraph,"NewtonGraph N_{n}(x) = y(x_{0})+#sum_{j = 1}^{n} f(x_{0}, ..., x_{j})#prod_{i = 0}^{j-1}(x - x_{i})", "lp");
    legend->Draw();

    //***********Конец оформления**********//



    //***********Часть2***********//
    TCanvas *canv2 = new TCanvas("Uncetarnities", "Uncetarnities");
	canv2->SetTitle("Uncetarnities");
	canv2->Divide(1, 2);
	canv2->cd();

    TF1 *UncertanityFunc = new TF1("Uncetarnities", "abs(OriginalFunc - PolinomFunc)", LEFT_LIMIT, RIGHT_LIMIT);
  
    canv2->cd(1);
    gPad->SetGrid();

    TF1 *OriginalFuncCopy = new TF1("OriginalFuncCopy", "OriginalFunc", LEFT_LIMIT, RIGHT_LIMIT);
    TF1 *PolinomFuncCopy = new TF1("PolinomFuncCopy", "PolinomFunc", LEFT_LIMIT, RIGHT_LIMIT);

    //***********Другое количество точек**********//
    double x2[NUM_POINTS_2];
	double y2[NUM_POINTS_2];

	FillX2(x2);
	FillY2(y2, x2);
    TF1 *Pol = PolMethod2(x2, y2); 
    TF1 *UncertanityFunc2 = new TF1("Uncetarnities2", "abs(OriginalFunc - Polinom)", LEFT_LIMIT, RIGHT_LIMIT);
    TGraph *OriginalFuncPoints2 = new TGraph(NUM_POINTS_2, x2, y2);
    //***********Другое количество точек**********//

    OriginalFuncCopy->SetTitle("Uncetarnities R(x) = |y(x) - #Phi(x)|");
	OriginalFuncCopy->GetXaxis()->SetTitle("X");
	OriginalFuncCopy->GetYaxis()->SetTitle("Y");
	OriginalFuncCopy->SetLineWidth(8);
	OriginalFuncCopy->Draw();

	OriginalFuncPoints->SetMarkerColor(kViolet);
	OriginalFuncPoints->SetMarkerSize(3);
	OriginalFuncPoints->SetMarkerStyle(45);
	OriginalFuncPoints->Draw("P same");

	OriginalFuncPoints2->SetMarkerColor(kTeal);
	OriginalFuncPoints2->SetMarkerSize(3);
	OriginalFuncPoints2->SetMarkerStyle(45);
	OriginalFuncPoints2->Draw("P same");

	PolinomFuncCopy->SetLineColor(kBlack);
	PolinomFuncCopy->SetLineStyle(1);
	PolinomFuncCopy->SetLineWidth(5);
	PolinomFuncCopy->Draw("same");

	Pol->SetLineColor(kGreen);
	Pol->SetLineStyle(1);
	Pol->SetLineWidth(2);
	Pol->Draw("same");

    TLegend *legend2 = new TLegend(0.30, 0.1, 0.60, 0.70);
    legend2->SetHeader("Uncetarnities","C"); 
    legend2->AddEntry(OriginalFuncCopy,"OriginalFunc y(x) = sin(x) #upoint e^{-x}", "lp");
    legend2->AddEntry(PolinomFuncCopy, Form("PolinomFunc #Phi_{n}(x) = #sum_{j = 0}^{n} a_{j} #upoint x^{j}, n = %d", NUM_POINTS), "lp");
    legend2->AddEntry(Pol, Form("PolinomFunc #Phi_{n}(x) = #sum_{j = 0}^{n} a_{j} #upoint x^{j}, n = %d", NUM_POINTS_2) , "lp");
    legend2->Draw();

    canv2->cd(2);
    gPad->SetGrid();

    TLegend *legend3 = new TLegend(0.35, 0.35, 0.65, 0.80);

    UncertanityFunc->SetTitle("R(x) = |y(x) - #Phi(x)|");
    UncertanityFunc->SetLineColor(kGreen);
	UncertanityFunc->SetLineStyle(1);
	UncertanityFunc->SetLineWidth(2);
	UncertanityFunc->GetYaxis()->SetRangeUser(0.0, 0.001);
    UncertanityFunc->Draw();

    UncertanityFunc2->SetLineColor(kBlue);
	UncertanityFunc2->SetLineStyle(1);
	UncertanityFunc2->SetLineWidth(2);
    UncertanityFunc2->Draw("same");

    legend3->AddEntry(UncertanityFunc, Form("UncertanityFunc R(x) = |y(x) - #Phi(x)|, n = %d", NUM_POINTS), "lp");
    legend3->AddEntry(UncertanityFunc2, Form("UncertanityFunc R(x) = |y(x) - #Phi(x)|, n = %d", NUM_POINTS_2), "lp");
    legend3->Draw();

    //***********Поиск максимума***********//

    int k = NUM_POINTS;

    double omega_NUM_POINTS[NUM_POLINOMS_POINTS];
    double R_max_NUM_POINTS[NUM_POLINOMS_POINTS];

    for (int i = 0; i < NUM_POLINOMS_POINTS; ++i)
    {
    	omega_NUM_POINTS[i] = 1;
    	for (int j = 0; j < NUM_POINTS; ++j)
    	{
    		omega_NUM_POINTS[i] *= (PolinomsPoints[i] - x[j]);
    	}
    	R_max_NUM_POINTS[i] = 16*sin(M_PI/4)*exp(-M_PI/4)*fabs(omega_NUM_POINTS[i])/Factorial(k);
    }

    TGraph *R_MAX = new TGraph(NUM_POLINOMS_POINTS, PolinomsPoints, R_max_NUM_POINTS);

    R_MAX->SetLineColor(kGreen);
	R_MAX->SetLineStyle(9);
	R_MAX->SetLineWidth(3);
	R_MAX->Draw("same C");

	legend3->AddEntry(R_MAX, Form("R(x) #leq #frac{max|y^{n+1}(x)|}{(n+1)!}#upoint|#omega_{n+1}(x)|, n = %d", NUM_POINTS), "lp");

	int k2 = NUM_POINTS_2;

    double omega_NUM_POINTS2[NUM_POLINOMS_POINTS];
    double R_max_NUM_POINTS2[NUM_POLINOMS_POINTS];

    for (int i = 0; i < NUM_POLINOMS_POINTS; ++i)
    {
    	omega_NUM_POINTS2[i] = 1;
    	for (int j = 0; j < NUM_POINTS_2; ++j)
    	{
    		omega_NUM_POINTS2[i] *= (PolinomsPoints[i] - x2[j]);
    	}
    	R_max_NUM_POINTS2[i] = 16*sin(M_PI/4)*exp(-M_PI/4)*fabs(omega_NUM_POINTS2[i])/Factorial(k2);
    }

    TGraph *R_MAX2 = new TGraph(NUM_POLINOMS_POINTS, PolinomsPoints, R_max_NUM_POINTS2);

    R_MAX2->SetLineColor(kBlue);
	R_MAX2->SetLineStyle(9);
	R_MAX2->SetLineWidth(3);
	R_MAX2->Draw("same C");

	legend3->AddEntry(R_MAX2, Form("R(x) #leq #frac{max|y^{n+1}(x)|}{(n+1)!}#upoint|#omega_{n+1}(x)|, n = %d", NUM_POINTS_2), "lp");

	return 0;
}

//***********Функции заполнения**********//

void FillX(double* x)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		x[i] = (RIGHT_BOARD - LEFT_BOARD)/(NUM_POINTS - 1) * i;
	}
}

void FillY(double* y, double* x)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		y[i] = sin(x[i])*exp(-x[i]);
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
			F(i, j) = pow(arr[i], j);
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

TGraph* LagrangeMethod(double* PolinomsPoints, double* x, double* y)
{
	double Lagrange[NUM_POLINOMS_POINTS];
	double L = 0;
	double y_l = 0;

	for (int k = 0; k < NUM_POLINOMS_POINTS; ++k)
	{
		L = 0;
		for (int i = 0; i < NUM_POINTS; ++i)
		{
			y_l = y[i];
			for (int j = 0; j < NUM_POINTS; ++j)
			{
				if (j != i)
					y_l *= (PolinomsPoints[k] - x[j])/(x[i] - x[j]);
			}
			L += y_l;
		}
		Lagrange[k] = L;
	}

	TGraph *LagrangeGraph = new TGraph(NUM_POLINOMS_POINTS, PolinomsPoints, Lagrange);

	return LagrangeGraph;
}


void DevDif(double* Differences, double* x, double* y)
{
	double Composition = 1;
	Differences[0] = 1;

	for (int k = 1; k < NUM_POINTS; ++k)
	{
		for (int j = 0; j <= k; ++j)
		{
			Composition = 1;
			for (int i = 0; i <= k; ++i)
			{
				if (i != j)
					Composition *= (x[j] - x[i]);
			}
			Differences[k] += y[j]/Composition;
		}
	}

	// for (int i = 0; i < NUM_POINTS; ++i)
	// {
	// 	printf("Differences[%d] = %f\n", i, Differences[i]);
	// }
}


TGraph* NewtonMethod(double* PolinomsPoints, double* Differences, double* x, double* y)
{
	double Newton[NUM_POLINOMS_POINTS];
	double Composition = 1;

	for (int k = 0; k < NUM_POINTS; ++k)
	{
		for (int j = 0; j < NUM_POINTS; ++j)
		{
			Composition = 1;
			for (int i = 0; i <= j - 1; ++i)
			{
				Composition *= PolinomsPoints[k] - x[i];
			}
			Newton[k] += Composition*Differences[j];
		}
	}

	TGraph *NewtonGraph = new TGraph(NUM_POLINOMS_POINTS, PolinomsPoints, Newton);

	return NewtonGraph;
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

void PrintPolinomsPoints(double* PolinomsPoints)
{
	printf("\n");
	for (int i = 0; i < NUM_POLINOMS_POINTS; ++i)
	{
		printf("PolinomsPoints[%d] = %f \n", i, PolinomsPoints[i]);
	}
	printf("\n");
}





void FillX2(double* x2)
{
	for (int i = 0; i < NUM_POINTS_2; ++i)
	{
		x2[i] = (RIGHT_BOARD - LEFT_BOARD)/(NUM_POINTS_2 - 1) * i;
	}

	// PrintX2(x2);
}

void FillY2(double* y2, double* x2)
{
	for (int i = 0; i < NUM_POINTS_2; ++i)
	{
		y2[i] = sin(x2[i])*exp(-x2[i]);
	}

	// PrintY2(y2);
}

void FillMatrixF2(TMatrixD& F, double* arr)
{
	for (int i = 0; i < NUM_POINTS_2; ++i)
	{
		for (int j = 0; j < NUM_POINTS_2; ++j)
		{
			F(i, j) = pow(arr[i], (double)j);
		}
	}
}

void FillVectorY2(TVectorD& Y, double* arr)
{
	for (int i = 0; i < NUM_POINTS_2; ++i)
	{
		Y(i) = arr[i];
	}
}

void PrintX2(double* x2)
{
	for (int i = 0; i < NUM_POINTS_2; ++i)
	{
		printf("x[%d] = %f\n", i, x2[i]);
	}
	printf("\n");
}

void PrintY2(double* y2)
{
	for (int i = 0; i < NUM_POINTS_2; ++i)
	{
		printf("y[%d] = %f\n", i, y2[i]);
	}
	printf("\n");
}

void PrintMatrixF2(TMatrixD F)
{
	for (int i = 0; i < NUM_POINTS_2; ++i)
	{
		for (int j = 0; j < NUM_POINTS_2; ++j)
		{
			printf("A(%d, %d) = %f \t", i, j, F(i, j));
		}
		printf("\n");
	}
	printf("\n");
}

void PrintVectorY2(TVectorD Y)
{
	for (int i = 0; i < NUM_POINTS_2; ++i)
	{
		printf("Y(%d) = %f\n", i, Y(i));
	}
	printf("\n");
}

void PrintVectorA2(TVectorD A)
{
	for (int i = 0; i < NUM_POINTS_2; ++i)
	{
		printf("A(%d) = %f\n", i, A(i));
	}
	printf("\n");
}

TF1* PolMethod2(double* x2, double* y2)
{
	TMatrixD F(NUM_POINTS_2, NUM_POINTS_2);
	FillMatrixF2(F, x2);

	//PrintMatrixF2(F);

	if(F.Determinant() == 0)
	{
		printf("Determinant of matrix is 0\n");
		exit(1);
		return 0;
	}

	TVectorD Y(NUM_POINTS_2);
	FillVectorY2(Y, y2);

	TVectorD A = F.Invert() * Y;

	//PrintVectorA2(A);

	string formula = "pol" + to_string(NUM_POINTS_2 - 1);

	TF1 *Polinom = new TF1("Polinom", formula.c_str(), LEFT_LIMIT, RIGHT_LIMIT);

	for (int i = 0; i < NUM_POINTS_2; ++i)
	{
		Polinom->SetParameter(i, A(i));
	}

	return Polinom;
}

int Factorial(int k)
{
	int result = 1;
	for (int i = 1; i <= k; ++i)
	{
		result *= i;
	}

	return result;
}

//***********Конец функций печати**********//
