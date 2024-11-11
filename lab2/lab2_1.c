#include <stdio.h>
#include <math.h>

#define NUM_POINTS  7

#define LEFT_BOARD 0.0
#define RIGHT_BOARD 3.0

#define LEFT_LIMIT -0.5
#define RIGHT_LIMIT 3.5

//#define NUM_POLINOMS_POINTS  1000

// #define POLINOM_LEFT_LIMIT -10.0
// #define POLINOM_RIGHT_LIMIT 10.0

void FillX(double*);
void FillY(double*, double*);
void FillH(double*, double*);
void FillMatrixSLAU(TMatrixD&, double*);
void FillDevDiff(double*, double*, double*);
void FillC(TMatrixD&, double*, double*);
void FillD(double*, double*, double*);
void FillB(double*, double*, double*, double*);
void FillSplain(TF1**, double*, double*, double*, double*, double*);
void FillDesign(TF1**, TLegend*, double*, double*, double*, double*, double*);
void FillUncertainities(TF1**, TF1**, TF1*, double*);
void FillUncertainitiesDesign(TF1**);

void FillMatrixSLAUThirdDerive(TMatrixD&, double*);
void FillCThirdDerivative(TMatrixD&, double*, double*);
void FillDThirdDerivative(double*, double*, double*);
void FillBTherdDerivative(double*, double*, double*, double*);
void FillDesignthThirdDerive(TF1**, double*, double*, double*, double*, double*);
void FillSplainTherdDerivative(TF1**, double*, double*, double*, double*, double*);
void FillUncertainitiesThirdDerive(TF1**, TF1**, TF1*, double*);
void FillUncertainitiesDesignThirdDerive(TF1**);

void PrintX(double*);
void PrintY(double*);
void PrintH(double*);
void PrintMatrixSLAU(TMatrixD&);
void PrintDevDiff(double*);
void PrintC(double*);
void PrintD(double*);
void PrintB(double*);

void PrintMatrixSLAUThirdDerive(TMatrixD&);
void PrintCThirdDerivative(double*);
void PrintDTherdDerivative(double*);
void PrintBTherdDerivative(double*);

int MainLab2_1()
{	
	printf("\n");

	double x[NUM_POINTS];
	double y[NUM_POINTS];

	FillX(x);
	FillY(y, x);

	PrintX(x);
	PrintY(y);

	double h[NUM_POINTS - 1];
	FillH(x, h);

	TMatrixD SLAU(NUM_POINTS, NUM_POINTS);
	FillMatrixSLAU(SLAU, h);

	PrintMatrixSLAU(SLAU);

	double DevDiff[NUM_POINTS];
	FillDevDiff(y, DevDiff, h);

	double c[NUM_POINTS];
	FillC(SLAU, DevDiff, c);

	PrintC(c);
	
	double d[NUM_POINTS - 1];
	FillD(c, d, h);

	PrintD(d);

	double b[NUM_POINTS - 1];
	FillB(b, y, h, c);
	PrintB(b);

	TCanvas *c1 = new TCanvas("c1", "c1");
	c1->SetTitle("Splain interpolation");
	c1->Divide(1, 2);
   	c1->cd(1);
   	gPad->SetGrid();

	TF1 *OriginalFunc = new TF1("OriginalFunc", "sin(x)*exp(-x)", LEFT_LIMIT, RIGHT_LIMIT);
	TGraph *OriginalFuncPoints = new TGraph(NUM_POINTS, x, y);

	OriginalFunc->SetTitle("Result of interpolation of function y(x) = sin(x)#upoint e^{-x}");
	OriginalFunc->GetXaxis()->SetTitle("X");
	OriginalFunc->GetYaxis()->SetTitle("Y");
	OriginalFunc->SetLineWidth(8);
	OriginalFunc->SetLineColor(NUM_POINTS);
	OriginalFunc->Draw();

	OriginalFuncPoints->SetMarkerColor(NUM_POINTS+2);
	OriginalFuncPoints->SetMarkerSize(3);
	OriginalFuncPoints->SetMarkerStyle(45);
	OriginalFuncPoints->Draw("P same");

	TF1 *Splains[NUM_POINTS - 1];
	FillSplain(Splains, y, b, c, d, x);

	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		cout << "s[" << i << "] = " << y[i] << " + " << b[i] << "*(x - " << x[i] << ") + " << c[i] << "* (x - " << x[i] << ")^2 + " << d[i] << "*(x - " << x[i] << ")^3" << endl;
	}

	TLegend *legend = new TLegend(0.30, 0.1, 0.9, 0.60);
    legend->SetHeader("Splain interpolation","C"); 
    legend->AddEntry(OriginalFunc,Form("OriginalFunc y(x) = sin(x) #upoint e^{-x}, x #in [%.3f, %.3f]", x[0], x[NUM_POINTS - 1]), "lp");
	FillDesign(Splains, legend, x, y, b, c, d);
	legend->Draw();

	//***********Погрешности**********//
	c1->cd(2);
	gPad->SetGrid();

	//Ещё костыль

	TF1 *func = new TF1("func", "0", LEFT_LIMIT, RIGHT_LIMIT);
	func->SetLineColor(1);
	func->SetLineWidth(0);
	func->GetYaxis()->SetRangeUser(0, 0.025);
	func->Draw();

	TF1 *Uncertainities[NUM_POINTS - 1];
	FillUncertainities(Uncertainities, Splains, OriginalFunc, x);
	FillUncertainitiesDesign(Uncertainities);

	//***********Сплайн с третьей производной**********//

	TMatrixD SLAUThirdDerive(NUM_POINTS, NUM_POINTS);
	FillMatrixSLAUThirdDerive(SLAUThirdDerive, h);

	PrintMatrixSLAUThirdDerive(SLAUThirdDerive);

	PrintDevDiff(DevDiff); 						

	double cThirdDerivative[NUM_POINTS];
	FillCThirdDerivative(SLAUThirdDerive, DevDiff, cThirdDerivative);

	PrintCThirdDerivative(cThirdDerivative);

	double dThirdDerivative[NUM_POINTS - 1];
	FillDThirdDerivative(cThirdDerivative, dThirdDerivative, h);

	PrintDTherdDerivative(dThirdDerivative);

	double bTherdDerivative[NUM_POINTS - 1];
	FillBTherdDerivative(bTherdDerivative, y, h, cThirdDerivative);
	PrintBTherdDerivative(bTherdDerivative);

	c1->cd(1);
	TF1 *SplainsThirdDerivative[NUM_POINTS - 1];
	FillSplainTherdDerivative(SplainsThirdDerivative, y, bTherdDerivative, cThirdDerivative, dThirdDerivative, x);
	FillDesignthThirdDerive(SplainsThirdDerivative, x, y, bTherdDerivative, cThirdDerivative, dThirdDerivative);


	c1->cd(2);	
	TF1 *UncertainitiesThirdDerive[NUM_POINTS - 1];
	FillUncertainitiesThirdDerive(UncertainitiesThirdDerive, SplainsThirdDerivative, OriginalFunc, x);
	FillUncertainitiesDesignThirdDerive(UncertainitiesThirdDerive);

	c1->SaveAs("CubicSplain.pdf");
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

void FillH(double* x, double* h)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		h[i] = x[i+1] - x[i];
	}
}

void FillMatrixSLAU(TMatrixD& SLAU, double* h)
{

	for (int i = 0; i < NUM_POINTS; ++i)
	{
		if (i == 0)
		{
			for (int j = 0; j < NUM_POINTS; ++j)
			{
				if (j == 0){
					SLAU(i, j) = 1;
					continue;
				}

				SLAU(i, j) = 0;
			}
			continue;
		}

		if (i == NUM_POINTS - 1)
		{
			for (int j = 0; j < NUM_POINTS; ++j)
			{
				if (j == NUM_POINTS - 1){
					SLAU(i, j) = 1;
					continue;
				}

				SLAU(i, j) = 0;
			}
			continue;
		}

		for (int j = 0; j < NUM_POINTS;)
		{
			if (j == i - 1 && j < NUM_POINTS - 2)
			{
				SLAU(i, j) = h[j];
				SLAU(i, j+1) = 2*(h[j] + h[j+1]);
				SLAU(i, j+2) = h[j+1];
				j+= 3;
			}
			else{
				SLAU(i, j) = 0;
				j++;
			}
		}
	}

}

void FillDevDiff(double* y, double* DevDiff, double* h)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		if (i == 0 || i == NUM_POINTS - 1){
			DevDiff[i] = 0;
			continue;
		}

		DevDiff[i] = 3*( (y[i+1] - y[i])/h[i] - (y[i] - y[i-1])/h[i-1] );
	}
}

void FillC(TMatrixD& SLAU, double* DevDiff, double* c)
{
	TVectorD DevDiffVec(NUM_POINTS);

	//Извиняюсь за костыли:)))

	for (int i = 0; i < NUM_POINTS; ++i)
	{
		DevDiffVec(i) = DevDiff[i];
	}

	TVectorD C = SLAU.Invert() * DevDiffVec;

	for (int i = 0; i < NUM_POINTS; ++i)
	{
		c[i] = C(i);
	}
}

void FillD(double* c, double* d, double* h)
{
	for (int i = 0; i < NUM_POINTS -1; ++i)
	{
		d[i] = (c[i+1] - c[i])/(3*h[i]);
	}
}

void FillB(double* b, double* y, double* h, double* c)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		b[i] = (y[i+1] - y[i])/h[i] - h[i]*(c[i+1] + 2*c[i])/3;
	}
}

void FillSplain(TF1** Splains, double* y, double* b, double* c, double* d, double* x)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		Splains[i] = new TF1(Form("Splain%d", i), "[0] + [1]*(x - [4]) + [2]*(x - [4])*(x - [4]) + [3]*(x - [4])*(x - [4])*(x - [4])", x[i], x[i+1]+0.01);
		Splains[i]->SetParameter(0, y[i]);
		Splains[i]->SetParameter(1, b[i]);
		Splains[i]->SetParameter(2, c[i]);
		Splains[i]->SetParameter(3, d[i]);
		Splains[i]->SetParameter(4, x[i]);
	}
}

void FillDesign(TF1** Splains, TLegend* legend, double* x, double* y, double* b, double* c, double* d)
{

	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{	
		Splains[i]->SetLineWidth(5);
		Splains[i]->SetLineColor(kBlue);
		Splains[i]->SetLineStyle(1);
		legend->AddEntry(Splains[i], Form("Splain[%d] = %.3f + %.3f*(x - %.3f) + %.3f*(x - %.3f)^{2} + %.3f*(x - %.3f)^{3}, x #in [%.3f, %.3f] ", i, y[i], b[i], x[i], c[i], x[i], d[i], x[i], x[i], x[i+1]), "lp");
		Splains[i]->Draw("same");
	}
}

void FillUncertainities(TF1** Uncertainities, TF1** Splains, TF1* OriginalFunc, double* x)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		Uncertainities[i] = new TF1(Form("Uncertainities%d", i), Form("abs(Splain%d - OriginalFunc)", i), x[i], x[i+1]);
	}
}

void FillUncertainitiesDesign(TF1** Uncertainities)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		Uncertainities[i]->SetLineWidth(2);
		Uncertainities[i]->SetLineColor(kBlue);
		Uncertainities[i]->SetLineStyle(1);
		Uncertainities[i]->Draw("same");
	}
}






void FillMatrixSLAUThirdDerive(TMatrixD& SLAUThirdDerive, double* h)
{

	for (int i = 0; i < NUM_POINTS; ++i)
	{
		if (i == 0)
		{
			for (int j = 0; j < NUM_POINTS;)
			{
				if (j == 0){
					SLAUThirdDerive(i, j) = 1/h[0];
					SLAUThirdDerive(i, j+1) = -1*(1/h[0] + 1/h[1]);
					SLAUThirdDerive(i, j+2) = 1/h[1];
					j+=3;
				}

				SLAUThirdDerive(i, j) = 0;
				j++;
			}
			continue;
		}

		if (i == NUM_POINTS - 1)
		{
			for (int j = 0; j < NUM_POINTS;)
			{
				if (j == NUM_POINTS - 3){
					SLAUThirdDerive(i, NUM_POINTS - 3) = 1/h[NUM_POINTS-3];
					SLAUThirdDerive(i, NUM_POINTS-2) = -1*(1/h[NUM_POINTS-3] + 1/h[NUM_POINTS-2]);
					SLAUThirdDerive(i, NUM_POINTS-1) = 1/h[NUM_POINTS-2];
					break;
				}

				SLAUThirdDerive(i, j) = 0;
				j++;
			}
			continue;
		}

		for (int j = 0; j < NUM_POINTS;)
		{
			if (j == i - 1 && j < NUM_POINTS - 2)
			{
				SLAUThirdDerive(i, j) = h[j];
				SLAUThirdDerive(i, j+1) = 2*(h[j] + h[j+1]);
				SLAUThirdDerive(i, j+2) = h[j+1];
				j+= 3;
			}
			else{
				SLAUThirdDerive(i, j) = 0;
				j++;
			}
		}
	}
}

void FillCThirdDerivative(TMatrixD& SLAUThirdDerive, double* DevDiff, double* cThirdDerivative)
{
	TVectorD DevDiffVec(NUM_POINTS);

	//Извиняюсь за костыли:)))

	for (int i = 0; i < NUM_POINTS; ++i)
	{
		DevDiffVec(i) = DevDiff[i];
	}

	TVectorD CThirdDerivative = SLAUThirdDerive.Invert() * DevDiffVec;

	for (int i = 0; i < NUM_POINTS; ++i)
	{
		cThirdDerivative[i] = CThirdDerivative(i);
	}
}

void FillDThirdDerivative(double* cThirdDerivative, double* dThirdDerivative, double* h)
{
	for (int i = 0; i < NUM_POINTS -1; i++)
	{
		dThirdDerivative[i] = (cThirdDerivative[i+1] - cThirdDerivative[i])/(3*h[i]);
	}
}

void FillBTherdDerivative(double* bTherdDerivative, double* y, double* h, double* cThirdDerivative)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		bTherdDerivative[i] = (y[i+1] - y[i])/h[i] - h[i]*(cThirdDerivative[i+1] + 2*cThirdDerivative[i])/3;
	}
}

void FillSplainTherdDerivative(TF1** SplainsThirdDerivative, double* y, double* bTherdDerivative, double* cThirdDerivative, double* dThirdDerivative, double* x)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		SplainsThirdDerivative[i] = new TF1(Form("SplainTherdDerivative%d", i), "[0] + [1]*(x - [4]) + [2]*(x - [4])*(x - [4]) + [3]*(x - [4])*(x - [4])*(x - [4])", x[i], x[i+1]+0.01);
		SplainsThirdDerivative[i]->SetParameter(0, y[i]);
		SplainsThirdDerivative[i]->SetParameter(1, bTherdDerivative[i]);
		SplainsThirdDerivative[i]->SetParameter(2, cThirdDerivative[i]);
		SplainsThirdDerivative[i]->SetParameter(3, dThirdDerivative[i]);
		SplainsThirdDerivative[i]->SetParameter(4, x[i]);
	}
}

void FillDesignthThirdDerive(TF1** SplainsThirdDerivative, double* x, double* y, double* bTherdDerivative, double* cThirdDerivative, double* dThirdDerivative)
{

	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{	
		SplainsThirdDerivative[i]->SetLineWidth(2);
		SplainsThirdDerivative[i]->SetLineColor(kOrange);
		SplainsThirdDerivative[i]->SetLineStyle(9);
		SplainsThirdDerivative[i]->Draw("same");
	}
}

void FillUncertainitiesThirdDerive(TF1** UncertainitiesThirdDerive, TF1** SplainsThirdDerivative, TF1* OriginalFunc, double* x)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		UncertainitiesThirdDerive[i] = new TF1(Form("UncertainitiesThirdDerive%d", i), Form("abs(SplainTherdDerivative%d - OriginalFunc)", i), x[i], x[i+1]);
	}
}

void FillUncertainitiesDesignThirdDerive(TF1** UncertainitiesThirdDerive)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		UncertainitiesThirdDerive[i]->SetLineWidth(2);
		UncertainitiesThirdDerive[i]->SetLineColor(kOrange);
		UncertainitiesThirdDerive[i]->SetLineStyle(9);
		UncertainitiesThirdDerive[i]->Draw("same");
	}
}
//***********Конец функций заполнения**********//




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

void PrintH(double* h)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		printf("h[%d] = %f\n", i, h[i]);
	}
	printf("\n");
}

void PrintMatrixSLAU(TMatrixD &SLAU)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		for (int j = 0; j < NUM_POINTS; ++j)
		{
			printf("SLAU(%d, %2d) = %.2f \t", i, j, SLAU(i, j));
		}
		printf("\n");
	}
	printf("\n");
}

void PrintDevDiff(double* DevDiff)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		printf("DevDiff[%d] = %f\n", i, DevDiff[i]);
	}
	printf("\n");
}

void PrintC(double* c)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		printf("c[%d] = %f\n", i, c[i]);
	}
	printf("\n");
}

void PrintD(double* d)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		printf("d[%d] = %f\n", i, d[i]);
	}
	printf("\n");
}

void PrintB(double* b)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		printf("b[%d] = %f\n", i, b[i]);
	}
	printf("\n");
}


void PrintMatrixSLAUThirdDerive(TMatrixD &SLAUThirdDerive)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		for (int j = 0; j < NUM_POINTS; ++j)
		{
			printf("SLAU3(%d, %2d) = %.2f \t", i, j, SLAUThirdDerive(i, j));
		}
		printf("\n");
	}
	printf("\n");
}

void PrintCThirdDerivative(double* cThirdDerivative)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		printf("cThirdDerivative[%d] = %f\n", i, cThirdDerivative[i]);
	}
	printf("\n");
}

void PrintDTherdDerivative(double* dThirdDerivative)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		printf("dThirdDerivative[%d] = %f\n", i, dThirdDerivative[i]);
	}
	printf("\n");
}

void PrintBTherdDerivative(double* bTherdDerivative)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		printf("bTherdDerivative[%d] = %f\n", i, bTherdDerivative[i]);
	}
	printf("\n");
}
//***********Конец функций печати**********//