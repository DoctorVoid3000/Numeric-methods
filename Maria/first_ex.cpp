// #include <TMatrixD.h>
// #include <TVectorD.h>
// #include <iostream>
// #include <Tmath.h>
// #include "TROOT.h"
// #include "TLegend.h"

Double_t func1(Double_t a, Double_t b, int n, int i, int j)
{
    float p = pow(a + i * (b - a) / n, j);
    return p;
}

Double_t func2(Double_t a, Double_t b, int n, int i, int s)
{
    float l = cos(3.14 / 2 / (n + 1)) * 2;
    float x;
    if (s == 0) {
        x = a + i * (b - a) / n;
        cout << "x " << x << endl;
    }
    else {
        x = (-1)* cos((3.14 / 2 + 3.14 * i) / (n + 1)) * (b - a) / l +(b + a)/2;
    }
    float p = sin(x) * exp(-x);
    return p;
}

Double_t func3(Double_t a, Double_t b, int n, Double_t x, TVectorD B, int s)
{
    float l = cos(3.14 / 2 / (n + 1)) * 2;
    Double_t fc = 0;
    Double_t mn = 1;
    float xi;
    float xj;
    for (int i = 0; i < n+1; i++) {
        if (s == 0) {
            xi = a + i * (b - a) / n;
        }
        else {
            xi = (-1)*cos((3.14 / 2 + 3.14 * i) / (n + 1)) * (b - a) / l +(b+a)/2;
        }
        for (int j = 0; j < n+1; j++) {
            if (s == 0) {
                xj = a + j * (b - a) / n;
            }
            else {
                xj = (-1) * cos((3.14 / 2 + 3.14 * j) / (n + 1)) * (b - a) / l + (b + a) / 2;
            }
            if (i != j) {
                mn = mn * (x - xj)/(xi-xj);
            }  
        }
        fc += B(i)*mn;
        mn = 1;
    }
    return fc;
}

Double_t func4(Double_t a, Double_t b, int n, Double_t x, TVectorD B)
{
    Double_t q[n];
    Double_t p = B(0);
    for (int i = 0; i < n; i++) {
        q[i] = (B(i+1) - B(i)) / ((b - a) / n);
    }
    float k = q[0];
    float r;
    for (int i = 1; i < n; i++) {
        for (int j = i; j < n; j++) {
            r = q[j];
            q[j] = (q[j] - k) / ((i+1)*(b - a) / n);
            k = r;
        }
        k = q[i];
    }
    for (int i = 1; i < n+1; i++) {
        for (int j = 0; j < i; j++) {
            q[i - 1] *= (x - a - (b - a) / n * j);
        }
        p += q[i - 1]; 
    }
    return p;
}

Double_t er(Double_t a, Double_t b, int n, float x) {
    float err = 1;
    for (int i = 0; i < n; i++) {
        err = err * abs((x - (a + i * (b - a) / (n - 1)))) / (i + 1);
    }
    err = err * pow(4, (n / 4));
    return err;
}

void first_ex() {

   TCanvas* c1 = new TCanvas();
   TCanvas* c2 = new TCanvas();
   TCanvas* c3 = new TCanvas();

   Double_t a = 0;
   Double_t b = 3;
   int n = 6;
   int n2 = n - 1;
   int z0 = 50;

   //first part of exercise
   Double_t answx0[z0];
   Double_t answy1[n], answy11[z0];
   Double_t answy2[n], answy21[z0];
   Double_t answy3[n], answy31[z0];
   Double_t answy[n], answy0[z0];
   Double_t answx[n];
   
   TMatrixD A(n, n);
   TVectorD B(n);
   for (int i = 0; i < n; i++) {
       for (int j = 0; j < n; j++) {
           A(i, j) = func1(a, b, n2, i, j);
       }
       B(i) = func2(a, b, n2, i, 0);
   }
   TVectorD w = A.Invert() * B; 
   

   for (int i = 0; i < n; i++) {
       answx[i] = a + i * (b - a) / n2;
       for (int j = 0; j < n; j++) {
           answy1[i]= answy1[i] + w(j) * pow(answx[i], j);
       }
   }
   for (int i = 0; i < n; i++) {
       answy2[i] = func3(a, b, n2, answx[i], B, 0);
       answy3[i] = func4(a, b, n2, answx[i], B);
       answy[i] = func2(a, b, n2, i, 0);
   }

   for (int i = 0; i < z0; i++) {
       answx0[i] = a + i * (b - a) / (z0 - 1);
       answy21[i] = func3(a, b, n2, answx0[i], B, 0);
       answy31[i] = func4(a, b, n2, answx0[i], B);
       for (int j = 0; j < n; j++) {
           answy11[i] = answy11[i] + w(j) * pow(answx0[i], j);
       }
       cout << "answy0" << endl;
       answy0[i] = func2(a, b, z0-1, i, 0);
   }

   c1->cd();
   TGraph* gr1 = new TGraph(n, answx, answy1);
   TGraph* gr2 = new TGraph(n, answx, answy2);
   TGraph* gr3 = new TGraph(n, answx, answy3);
   TGraph* gr4 = new TGraph(n, answx, answy);
   TGraph* gr1_1 = new TGraph(z0, answx0, answy11);
   TGraph* gr2_1 = new TGraph(z0, answx0, answy21);
   TGraph* gr3_1 = new TGraph(z0, answx0, answy31);
   TGraph* gr4_1 = new TGraph(z0, answx0, answy0);
   TMultiGraph* mg = new TMultiGraph();
   mg->Add(gr1);
   mg->Add(gr2);
   mg->Add(gr3);
   mg->Add(gr4);
   mg->Add(gr1_1);
   mg->Add(gr2_1);
   mg->Add(gr3_1);
   mg->Add(gr4_1);
   gr1->SetMarkerStyle(24);
   gr1->SetMarkerSize(2);
   gr1->SetLineWidth(0);
   gr2->SetMarkerStyle(24);
   gr2->SetMarkerSize(4);
   gr2->SetMarkerColor(kGreen);
   gr2->SetLineColor(kGreen);
   gr2->SetLineWidth(0);
   gr3->SetMarkerStyle(24);
   gr3->SetMarkerSize(6);
   gr3->SetMarkerColor(kBlue);
   gr3->SetLineColor(kBlue);
   gr3->SetLineWidth(0);
   gr4->SetMarkerStyle(24);
   gr4->SetMarkerSize(8);
   gr4->SetMarkerColor(kRed);
   gr4->SetLineColor(kRed);
   gr4->SetLineWidth(0);
   gr2_1->SetLineColor(kGreen);
   gr3_1->SetLineColor(kBlue);
   gr4_1->SetLineColor(kRed);
   mg->Draw("ACP");
   mg->SetTitle("Function interpolation");

   auto legend = new TLegend(0.6, 0.7, 0.8, 0.9);
   legend->SetHeader("Legend");
   legend->AddEntry(gr1_1, "Polynomial", "l");
   legend->AddEntry(gr2_1, "Lagrange polynomial", "l");
   legend->AddEntry(gr3_1, "Newton polynomial","l");
   legend->AddEntry(gr4_1, "Function", "l");
   legend->Draw();

   //second a part
   c2->Divide(1, 2);
   c2->cd(1);
   int z = 200;
   int n1 = 10;
   int n3 = 6;
   Double_t y1[z], y2[z], x[z];
   TVectorD B1(n1), B2(n3);
   Double_t x2_1[n1], x2_2[n3], y2_1[n1], y2_2[n3];

   for (int i = 0; i < n1; i++) {
       B1(i) = func2(a, b, n1-1, i, 0);
       x2_1[i] = a + i * (b - a) / (n1 - 1);
       y2_1[i] = func2(a, b, n1 - 1, i, 0);
   }
   for (int i = 0; i < n3; i++) {
       B2(i) = func2(a, b, n3-1, i, 0);
       x2_2[i] = a + i * (b - a) / (n3 - 1);
       y2_2[i] = func2(a, b, n3 - 1, i, 0);
   }

   for (int i = 0; i < z; i++) {
       x[i] =  a + i * (b - a) / (z-1);
       y1[i] = func3(a, b, n1 - 1, x[i], B1, 0);
       y2[i] = func3(a, b, n3 - 1, x[i], B2, 0);
   }

   TGraph* gr5 = new TGraph(z, x, y1);
   TGraph* gr6 = new TGraph(z,x, y2);
   TGraph* gr5_1 = new TGraph(n1, x2_1, y2_1);
   TGraph* gr6_1 = new TGraph(n3, x2_2, y2_2);
   TMultiGraph* mg2 = new TMultiGraph();
   mg2->Add(gr5);
   mg2->Add(gr6);
   mg2->Add(gr5_1);
   mg2->Add(gr6_1);
   gr5_1->SetLineWidth(0);
   gr5_1->SetMarkerStyle(24);
   gr5_1->SetMarkerSize(2);
   gr6_1->SetLineWidth(0);
   gr6_1->SetMarkerStyle(24);
   gr6_1->SetMarkerSize(2);
   gr6_1->SetMarkerColor(kGreen);
   gr6->SetLineColor(kGreen);
   mg2->Draw("ACP");
   mg2->SetTitle("Function interpolation");

   auto legend2 = new TLegend(0.6, 0.7, 0.8, 0.9);
   legend2->SetHeader("Legend");
   legend2->AddEntry(gr5, "Lagrange polynomial, n = 6");
   legend2->AddEntry(gr6, "Lagrange polynomial, n = 7");
   legend2->Draw();

   //second b part
   c2->cd(2);
   Double_t ery1[z], ery2[z], ery1R[z], ery2R[z];
   for (int i = 0; i < z; i++) {
       ery1[i] = abs(y1[i] - sin(x[i]) * exp(-x[i]));
       ery2[i] = abs(y2[i] - sin(x[i]) * exp(-x[i]));
       ery1R[i] = er(a, b, n1, x[i]);
       ery2R[i] = er(a, b, n3, x[i]);
   }
   TGraph* gr7 = new TGraph(z, x, ery1);
   TGraph* gr8 = new TGraph(z, x, ery2);
   TGraph* gr9 = new TGraph(z, x, ery1R);
   TGraph* gr10 = new TGraph(z, x, ery2R);
   TMultiGraph* mg3 = new TMultiGraph();
   mg3->Add(gr7);
   gr7->SetLineStyle(2);
   mg3->Add(gr8);
   gr8->SetLineStyle(2);
   gr8->SetLineColor(kGreen);
   mg3->Add(gr9);
   mg3->Add(gr10);
   gr10->SetLineColor(kGreen);
   mg3->Draw("ACP");
   mg3->SetTitle("Function errors");

   auto legend3 = new TLegend(0.6, 0.7, 0.8, 0.9);
   legend3->SetHeader("Legend");
   legend3->AddEntry(gr7, "Real function of errors for Lp (n = 6)", "lp");
   legend3->AddEntry(gr8, "Real function of errors for Lp (n = 7)");
   legend3->AddEntry(gr9, "Estimated error (n = 6)");
   legend3->AddEntry(gr10, "Estimated error (n = 7)");
   legend3->Draw();

   //third part 

   c3->Divide(1, 2);
   c3->cd(1);
   Double_t ty1[z], ty2[z];
   TVectorD B3_1(n3), B3_2(n3);
   Double_t x3_1[n3], y3_1[n3], x3_2[n3], y3_2[n3];
   float l = cos(3.14 / 2 / (n + 1)) * 2;

   for (int i = 0; i < n3; i++) {
       x3_1[i] = a + i * (b - a) / (n3 - 1);
       y3_1[i] = func2(a, b, n3 - 1, i, 0);
       x3_2[i] = (-1) * cos((3.14 / 2 + 3.14 * i) / (n3)) * (b - a) / l + (b + a) / 2;
       y3_2[i] = func2(a, b, n3 - 1, i, 1);
   }

   for (int i = 0; i < n3; i++) {
       B3_1(i) = func2(a, b, n3 - 1, i, 0);
       B3_2(i) = func2(a, b, n3 - 1, i, 1);
       cout << "B3_2 " << B3_2(i) << endl;
   }
   for (int i = 0; i < z; i++) {
       x[i] = a + i * (b - a) / (z - 1);
       ty1[i] = func3(a, b, n3 - 1, x[i], B3_1, 0);
       ty2[i] = func3(a, b, n3 - 1, x[i], B3_2, 1);
   }
   TGraph* gr11 = new TGraph(z, x, ty1);
   TGraph* gr12 = new TGraph(z, x, ty2);
   TGraph* gr11_1 = new TGraph(n3, x3_1, y3_1);
   TGraph* gr12_1 = new TGraph(n3, x3_2, y3_2);
   TMultiGraph* mg4 = new TMultiGraph();
   mg4->Add(gr11);
   mg4->Add(gr11_1);
   gr11_1->SetLineWidth(0);
   gr11_1->SetMarkerStyle(24);
   gr11_1->SetMarkerSize(2);
   gr11_1->SetMarkerColor(kGreen);
   gr11->SetLineColor(kGreen);
   mg4->Add(gr12);
   mg4->Add(gr12_1);
   gr12_1->SetLineWidth(0);
   gr12_1->SetMarkerStyle(24);
   gr12_1->SetMarkerSize(2);
   mg4->Draw("ACP");
   auto legend4 = new TLegend(0.6, 0.7, 0.8, 0.9);
   legend4->SetHeader("Legend");
   legend4->AddEntry(gr11, "Uniform grid");
   legend4->AddEntry(gr12, "Chebyshev grid");
   legend4->Draw();

   c3->cd(2);
   Double_t erty1[z], erty2[z];
   for (int i = 0; i < z; i++) {
       erty1[i] = abs(ty1[i] - sin(x[i]) * exp(-x[i]));
       erty2[i] = abs(ty2[i] - sin(x[i]) * exp(-x[i]));
   }
   TGraph* gr13 = new TGraph(z, x, erty1);
   TGraph* gr14 = new TGraph(z, x, erty2);
   TMultiGraph* mg5 = new TMultiGraph();
   mg5->Add(gr13);
   gr13->SetLineColor(kGreen);
   mg5->Add(gr14);
   mg5->Draw("ACP");
   auto legend5 = new TLegend(0.6, 0.7, 0.8, 0.9);
   legend5->SetHeader("Legend");
   legend5->AddEntry(gr13, "Error uniform grid");
   legend5->AddEntry(gr14, "Error Chebyshev grid");
   legend5->Draw();


}