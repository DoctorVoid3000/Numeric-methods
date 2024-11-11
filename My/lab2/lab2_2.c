#include <stdio.h>
#include <math.h>

#define NUM_POINTS  7

#define LEFT_BOARD 0.0
#define RIGHT_BOARD 3.0

#define LEFT_LIMIT -0.5
#define RIGHT_LIMIT 3.5

void FillX(TVectorD&);
void FillY(TVectorD&, TVectorD&);
void FillSplain(TF1**, TVectorD&, TVectorD&);
void FillDesign(TF1**);
void FillUncertainities(TF1**, TF1*, TVectorD&);
void FillUncertainitiesDesign(TF1**);
void FillKostyl(double*, double*, TVectorD&, TVectorD&);
void FillX3(TVectorD&);
void FillY3(TVectorD&, TVectorD&);
void FillMatrixB3(TMatrixD&, TVectorD&);
void FillA3(TVectorD&, TMatrixD&, TVectorD&);
void FillB3Splain(TF1**, TVectorD&, TVectorD&);
void FillDesignB3(TF1**);
void FillUncertainitiesB3(TF1**, TF1*, TVectorD&);
void FillUncertainitiesDesignB3(TF1**);

void PrintX(TVectorD&);
void PrintY(TVectorD&);
void PrintX3(TVectorD&);
void PrintY3(TVectorD&);
void PrintMatrixB3(TMatrixD&);
void PrintA3(TVectorD&);


int MainBSplain(){

	//***********B1Splains***********

	printf("\n");

	TVectorD x(NUM_POINTS);
	TVectorD y(NUM_POINTS);

	FillX(x);
	FillY(y, x);

	// PrintX(x);
	// PrintY(y);

	TCanvas *c1 = new TCanvas("c1", "c1");
	c1->SetTitle("B-Splain interpolation");
	c1->Divide(1, 2);
   	c1->cd(1);
   	gPad->SetGrid();

   	//Опять костыль :(( И всё из-за того что TGraph не принимает TVector...

   	double KostylX[NUM_POINTS];
   	double KostylY[NUM_POINTS];

   	FillKostyl(KostylX, KostylY, x, y);

	//конец костыля

	TF1 *OriginalFunc = new TF1("OriginalFunc", "sin(x)*exp(-x)", LEFT_LIMIT, RIGHT_LIMIT);
	TGraph *OriginalFuncPoints = new TGraph(NUM_POINTS, KostylX, KostylY);

	OriginalFunc->SetTitle("Result of B-splain interpolation of function y(x) = sin(x)#upoint e^{-x}");
	OriginalFunc->GetXaxis()->SetTitle("X");
	OriginalFunc->GetYaxis()->SetTitle("Y");
	OriginalFunc->SetLineWidth(8);
	OriginalFunc->SetLineColor(NUM_POINTS);
	OriginalFunc->Draw();

	OriginalFuncPoints->SetMarkerColor(NUM_POINTS+2);
	OriginalFuncPoints->SetMarkerSize(3);
	OriginalFuncPoints->SetMarkerStyle(45);
	OriginalFuncPoints->Draw("P same");

    TF1 *B1Splains[NUM_POINTS-1];

	FillSplain(B1Splains, x, y);

	FillDesign(B1Splains);

	//***********B1Splains***********

	


	//***********B3Splains***********
	TVectorD X3(NUM_POINTS+2);
	FillX3(X3);

	PrintX3(X3);

	TVectorD Y3(NUM_POINTS+2);
	FillY3(Y3, x);

	PrintY3(Y3);

	TMatrixD B3(NUM_POINTS+2, NUM_POINTS+2);
	FillMatrixB3(B3, X3);

	PrintMatrixB3(B3);

	TVectorD A3(NUM_POINTS+2);
	FillA3(A3, B3, Y3);

	PrintA3(A3);

	TF1 *B3Splains[NUM_POINTS-1];
	FillB3Splain(B3Splains, X3, A3);

	FillDesignB3(B3Splains);
	//***********B3Splains***********

	TLegend *legend = new TLegend(0.30, 0.1, 0.9, 0.40);
    legend->SetHeader("B-Splain interpolation","C"); 
    legend->AddEntry(OriginalFunc,Form("OriginalFunc y(x) = sin(x) #upoint e^{-x}, x #in [%.3f, %.3f]", x[0], x[NUM_POINTS - 1]), "lp");
	legend->AddEntry(B1Splains[0], "B1-splain");
	legend->AddEntry(B3Splains[0], "B3-splain");
    legend->Draw();


	//***********Погрешности**********//
	c1->cd(2);
	gPad->SetGrid();

	//Ещё костыль

	TF1 *func = new TF1("func", "0", LEFT_LIMIT, RIGHT_LIMIT);
	func->SetLineColor(1);
	func->SetLineWidth(0);
	func->GetYaxis()->SetRangeUser(0, 0.065);
	func->Draw();

	TF1 *B1Uncertainities[NUM_POINTS - 1];
	FillUncertainities(B1Uncertainities, OriginalFunc, x);
	FillUncertainitiesDesign(B1Uncertainities);

	TF1 *B3Uncertainities[NUM_POINTS - 1];
	FillUncertainitiesB3(B3Uncertainities, OriginalFunc, X3);
	FillUncertainitiesDesignB3(B3Uncertainities);

	return 0;
}


//***********Функции заполнения**********//
void FillX(TVectorD& x)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		x(i) = (RIGHT_BOARD - LEFT_BOARD)/(NUM_POINTS - 1) * i;
	}
}

void FillY(TVectorD& y, TVectorD& x)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		y(i) = sin(x(i))*exp(-x(i));
	}
}

void FillSplain(TF1** B1Splains, TVectorD& x, TVectorD& y)
{
	for (int i = 0; i < NUM_POINTS-1; ++i)
	{
		B1Splains[i] = new TF1(Form("B1Splain%d", i), "[0]*([3] - x)/([3] - [2]) + [1]*(x-[2])/([3]-[2])", x[i], x[i+1]);
		B1Splains[i]->SetParameter(0, y(i));
		B1Splains[i]->SetParameter(1, y(i+1));
		B1Splains[i]->SetParameter(2, x(i));
		B1Splains[i]->SetParameter(3, x(i+1));

	}
}

void FillDesign(TF1** B1Splains)
{

	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{	
		B1Splains[i]->SetLineWidth(5);
		B1Splains[i]->SetLineColor(kBlue);
		B1Splains[i]->SetLineStyle(1);
		B1Splains[i]->Draw("same");
	}
}

void FillUncertainities(TF1** B1Uncertainities, TF1* OriginalFunc, TVectorD& x)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		B1Uncertainities[i] = new TF1(Form("B1Uncertainities%d", i), Form("abs(B1Splain%d - OriginalFunc)", i), x[i], x[i+1]);
	}
}

void FillUncertainitiesDesign(TF1** B1Uncertainities)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		B1Uncertainities[i]->SetLineWidth(2);
		B1Uncertainities[i]->SetLineColor(kBlue);
		B1Uncertainities[i]->SetLineStyle(1);
		B1Uncertainities[i]->Draw("same");
	}
}





void FillKostyl(double* KostylX, double* KostylY, TVectorD& x, TVectorD& y)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		KostylX[i] = x(i);
		KostylY[i] = y(i);
	}
}

void FillX3(TVectorD& X3)
{
	for (int i = 0; i < NUM_POINTS+2; ++i)
	{
		X3(i) = (RIGHT_BOARD - LEFT_BOARD)/(NUM_POINTS - 1) * i - (RIGHT_BOARD - LEFT_BOARD)/(NUM_POINTS - 1);
	}	
}

void FillY3(TVectorD& Y3, TVectorD& x)
{
	for (int i = 0; i < NUM_POINTS+2; ++i)
	{
		if (i == 0 || i == NUM_POINTS+1)
		{
			Y3(i) = 0;
			continue;
		}
		Y3(i) = sin(x(i-1))*exp(-x(i-1));
	}	
}

void FillMatrixB3(TMatrixD& B3, TVectorD& X3)
{
	for (int i = 0; i < NUM_POINTS+2; ++i)
	{	
		if (i == 0)
		{
			B3(i, 0) = (2 - (X3(1) - X3(0)) / (X3(1)-X3(0)));
			B3(i, 1) = -2;
			B3(i, 2) = (2 + (X3(1) - X3(2)) / (X3(1)-X3(0)));
			for (int j = 3; j < NUM_POINTS+2; ++j)
			{
				B3(i, j) = 0;
			}
		}
		if (i==NUM_POINTS+1)
		{
			B3(i, NUM_POINTS - 1) = (2 - (X3(NUM_POINTS) - X3(NUM_POINTS-1)) / (X3(1)-X3(0)));
			B3(i, NUM_POINTS) = -2;
			B3(i, NUM_POINTS + 1) = (2 + (X3(NUM_POINTS) - X3(NUM_POINTS+1)) / (X3(1)-X3(0)));
			for (int j = 0; j < NUM_POINTS-1; ++j)
			{
				B3(i, j) = 0;
			}
		}
		if (i!=0 && i!=NUM_POINTS+1)
		{
			for (int j = 0; j < NUM_POINTS+2;)
			{	
				if (j == i-1 && j < NUM_POINTS)
				{
					B3(i, j) = 1/6.0;
					B3(i, j+1) = 2/3.0;
					B3(i, j+2) = 1/6.0;
					j+=3;
				}
				else{
				B3(i, j) = 0;
				j++;}
			}
		}
	}
}

void FillA3(TVectorD& A3, TMatrixD& B3, TVectorD& Y3)
{
	A3 = B3.Invert()*Y3;	
}

void FillB3Splain(TF1** B3Splains, TVectorD& X3, TVectorD& A3)
{
	for (int i = 0; i < NUM_POINTS-1; ++i)
	{
		B3Splains[i] = new TF1(Form("B3Splain%d", i), ("[0] * (1.0/6.0) * ( 2 - (x - [4])/[8] ) * ( 2 - (x - [4])/[8] ) * ( 2 - (x - [4])/[8] ) + "+ string("") +
													   "[1] * ( (2.0/3.0) - (x - [5])/[8]*(x - [5])/[8] + 0.5 * (x - [5])/[8]*(x - [5])/[8]*(x - [5])/[8] ) + "+
												       "[2] * ( (2.0/3.0) - (x - [6])/[8]*(x - [6])/[8] - 0.5 * (x - [6])/[8]*(x - [6])/[8]*(x - [6])/[8] ) + "+
												       "[3] * (1.0/6.0) * ( 2 + (x - [7])/[8] ) * ( 2 + (x - [7])/[8] ) * ( 2 + (x - [7])/[8] )").c_str(), X3(i+1), X3(i+2));
		B3Splains[i]->SetParameter(0, A3(i));
		B3Splains[i]->SetParameter(1, A3(i+1));
		B3Splains[i]->SetParameter(2, A3(i+2));
		B3Splains[i]->SetParameter(3, A3(i+3));
		B3Splains[i]->SetParameter(4, X3(i));
		B3Splains[i]->SetParameter(5, X3(i+1));
		B3Splains[i]->SetParameter(6, X3(i+2));
		B3Splains[i]->SetParameter(7, X3(i+3));
		B3Splains[i]->SetParameter(8, X3(i+1) - X3(i));
	}
}


void FillDesignB3(TF1** B3Splains)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{	
		B3Splains[i]->SetLineWidth(3);
		B3Splains[i]->SetLineColor(kMagenta);
		B3Splains[i]->SetLineStyle(1);
		
		B3Splains[i]->Draw("same");
	}
}

void FillUncertainitiesB3(TF1** B3Uncertainities, TF1* OriginalFunc, TVectorD& X3)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		B3Uncertainities[i] = new TF1(Form("B3Uncertainities%d", i), Form("abs(B3Splain%d - OriginalFunc)", i), X3(i+1), X3(i+2));
	}
}

void FillUncertainitiesDesignB3(TF1** B3Uncertainities)
{
	for (int i = 0; i < NUM_POINTS - 1; ++i)
	{
		B3Uncertainities[i]->SetLineWidth(2);
		B3Uncertainities[i]->SetLineColor(kMagenta);
		B3Uncertainities[i]->SetLineStyle(1);
		B3Uncertainities[i]->Draw("same");
	}
}
//***********Конец функций заполнения**********//


//***********Функции печати**********//

void PrintX(TVectorD& x)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		printf("x(%d) = %f\n", i, x(i));
	}
	printf("\n");
}

void PrintY(TVectorD& y)
{
	for (int i = 0; i < NUM_POINTS; ++i)
	{
		printf("y(%d) = %f\n", i, y(i));
	}
	printf("\n");
}

void PrintX3(TVectorD& X3)
{
	for (int i = 0; i < NUM_POINTS+2; ++i)
	{
		printf("X3(%d) = %f\n", i, X3(i));
	}
	printf("\n");
}

void PrintY3(TVectorD& Y3)
{
	for (int i = 0; i < NUM_POINTS+2; ++i)
	{
		printf("Y3(%d) = %f\n", i, Y3(i));
	}
	printf("\n");
}

void PrintMatrixB3(TMatrixD &B3)
{
	for (int i = 0; i < NUM_POINTS+2; ++i)
	{
		for (int j = 0; j < NUM_POINTS+2; ++j)
		{
			printf("B3(%d, %2d) = %.2f   ", i, j, B3(i, j));
		}
		printf("\n");
	}
	printf("\n");
}

void PrintA3(TVectorD& A3)
{
	for (int i = 0; i < NUM_POINTS+2; ++i)
	{
		printf("A3(%d) = %f\n", i, A3(i));
	}
	printf("\n");
}
//***********Конец функций печати**********//