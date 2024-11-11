struct FuncOptions{                                                                 //Создание структуры с набором параметров, задаваемых пользователем
    long double xmin, xmax;
    double dl, dr;
    int N, typePol;
    string func;
    bool OnFunc;
    vector<string> names;
};

// template <typename T> void TestFunc(T *x, T *y, struct FuncOptions &opt, int same = 1);
////////////////////////////////Part 1//////////////////////////////////////////////
template <typename T> void Part1(struct FuncOptions &opt, T type);                                      //Пункт 1.3
template <typename T> T** GetPointsUniform(T flag, struct FuncOptions &opt);                            //Получение координат точек функции для равномерной сетки
template <typename T> double* SomeToDobleMass(T* mass, int N);                                          //Преобразование массива значений в массив double

template <typename T> TF1** Get4Functions(T* x, T* y, struct FuncOptions &opt, bool key = false, int dcount = 1);          //Создание полиномов тремя способами
void SettingUpFunc(TF1** funcMass, TLegend* leg, int k, struct FuncOptions &opt);                       //Настройка визуальных параметров функций

template <typename T> TF1* ProstoPolin(T* x, T* y, int N);                                              //Получение полинома методом СЛАУ

template <typename T> TF1* LagranPolin(T* x, T* y, int N, int count = 0);                               //Получение полинома методом Лагранжа
template <typename T> T* PolLag(T* x, int j, int N);                                                    //Получение малого полинома в методе Лагранжа

template <typename T> TF1* NewtonPolin(T* x, T* y, int N, int count = 0);                               //Получение полинома методом Ньютона
template <typename T> T* PolNew(T* x, int j);                                                           //Получение малого полинома в методе Ньютона
template <typename T> T* FNew(T* x, T* y, int N);                                                       //Получение разделённых разностей

template <typename T> 
void SummMassive(T* arr1, int N1, T* arr2, int N2, long double k1, long double k2);                     //Сложение массивов

Double_t LagranFunc(Double_t *x, Double_t *par);                                                        //Функция, дающая значение полинома Лагранжа в точке x
Double_t NewtonFunc(Double_t *x, Double_t *par);                                                        //Функция, дающая значение полинома Ньютона в точке x

////////////////////////////////Part 2//////////////////////////////////////////////
template <typename T> void Part2(struct FuncOptions &opt, T type);                                      //Пункт 1.2
TF1* GetDeriveNforExpSinOrCos(double ksin, double kcos, int N);                                         //Возвращает функцию производной

template <typename T> TF1* RealDiff(T* x, T* y, struct FuncOptions &opt, bool key = false, int dcount = 1); //Возвращает реальзную разницу между функцией и интерполяцией
Double_t LagranDiv(Double_t *x, Double_t *par);                                                         //Возвращает реальзную разницу между функцией и интерполяцией Лагранжа
Double_t NewtonDiv(Double_t *x, Double_t *par);                                                         //Возвращает реальзную разницу между функцией и интерполяцией Ньютона

template <typename T> TF1* PolinDiv(T* x, T* y, struct FuncOptions &opt, bool key = false, int dcount = 1); //Возвращает оценочную разницу между функцией и интерполяцией
Double_t FuncDiv(Double_t *x, Double_t *par);                                                           //Возвращает оценочную разницу между функцией и интерполяцией

long double EstimMax(struct FuncOptions opt);                                                           //Возвращает максимальный модуль производной от функции 

////////////////////////////////Part 3//////////////////////////////////////////////
template <typename T> void Part3(struct FuncOptions &opt, T type);                                      //Пункт 1.3
template <typename T> T** GetPointsChebyshev(T flag, struct FuncOptions &opt);                          //Получение координат точек функции для Чебышевской сетки

void Chisl1(){
    FuncOptions opt;                                                                    //Дано
/////////For part 1/////
    opt.xmin = 0;
    opt.xmax = 3; //12
    opt.dl = 0.5; opt.dr = 10; //0.5 10
    opt.N = 7; //25
    opt.func = "exp(-x)*sin(x)";
    opt.OnFunc = true;
    opt.names = {"TEOR: "+opt.func, "SLAE", "Lagran", "Newton"};
    typedef long double myType;
    myType type;
/////////For part 2/////
    opt.typePol = 4;
////////////////////////

    Part1(opt, type);
    Part2(opt, type);
    Part3(opt, type);

/////////////////////////////////TEST///////////////////////////////////////////////
    // TCanvas* cTest = new TCanvas();                                                     //Отрисовка графиков и функций
    // cTest->SetTicks();
    // cTest->SetGrid();
    // gr1->Draw("A*same");
    // funcMass1[3]->Draw("same");

    // cTest->cd();
    // fDivReal1->Draw("same");
    // TestFunc(&x[0], &y[0], opt);
    // TestFunc(&x[0], &y[0], opt, 6);
    // Double_t test = 2.3;
    // Double_t par[] = {1, 5, 4, 0, 1, 2, 3};

    // cout<<FuncDiv(&test, par)<<endl;
////////////////////////////////////////////////////////////////////////////////////
}

////////////////////////////////Functions///////////////////////////////////////////
// template <typename T> void TestFunc(T *x, T *y, struct FuncOptions &opt, int same){
//     cout<<same<<endl;
// }

template <typename T> void Part1(struct FuncOptions &opt, T type){
    auto **points = GetPointsUniform(type, opt);                                        //Генерация векторов c координатами точек
    auto *x = points[0], *y = points[1];                                                //Разделенеи координат на два массива
    auto *x_gr = SomeToDobleMass(x, opt.N), *y_gr = SomeToDobleMass(y, opt.N);          //Создание массивов для графика

    TF1 **funcMass1 = Get4Functions(&x[0], &y[0], opt, true);                           //Создание полиномов тремя способами

    TGraph* gr1 = new TGraph(opt.N, &x_gr[0], &y_gr[0]);                                //Объявление графика с заданными точками
    gr1->SetTitle(("Interpolation of "+to_string(opt.N)+" points in the range ["+to_string(int(opt.xmin))+":"+to_string(int(opt.xmax))+"]").c_str());
    gr1->GetXaxis()->SetLimits(opt.xmin-opt.dl, opt.xmax+opt.dr);

    TLegend* leg1 = new TLegend(0.65, 0.6, 0.85, 0.8);                                  //Создание легенды
    leg1->SetBorderSize(0);

    SettingUpFunc(funcMass1, leg1, 4, opt);                                             //Настройка визуальных параметров функций и прикрепление их к легенде

    TCanvas* c1 = new TCanvas();                                                        //Отрисовка графиков и функций
    c1->SetTicks();
    c1->SetGrid();
    gr1->Draw("A*same");
    for(size_t i = 0; i<4; i++){ if(funcMass1[i]){ funcMass1[i]->Draw("same"); } }
    leg1->Draw("same");
}

template <typename T> T** GetPointsUniform(T flag, struct FuncOptions &opt){
    TF1 *f = new TF1("f", opt.func.c_str());                                            //Объявление необходимых объектов и переменных
    T *x = new T[opt.N], *y = new T[opt.N];
    T **res = new T*[2];

    long double d = (opt.xmax-opt.xmin)*1./(opt.N-1);                                   //Вычисление промежуточных значений

    for(size_t i = 0; i<opt.N; i++){                                                    //Заполнение векторов координатами точек
        x[i] = i*d+opt.xmin;
        y[i] = f->Eval(x[i]);
    }
    res[0] = x;
    res[1] = y;
    return res;
}

template <typename T> double* SomeToDobleMass(T* mass, int N){
    double *res = new double[N];
    for(size_t i = 0; i<N; i++){res[i] = mass[i];}
    return res;
}

template <typename T> TF1** Get4Functions(T* x, T* y, struct FuncOptions &opt, bool key, int dcount){
    TF1 **funcMass = new TF1*[opt.N];                                                                           //Объявление необходимых объектов и переменных
    static int count = 0;   
    if(key){count = 0;}
    else{count += dcount;}

    funcMass[0] = new TF1(("f0"+to_string(count)).c_str(), opt.func.c_str(), opt.xmin-opt.dl-1, opt.xmax+opt.dr+1);
    funcMass[1] = ProstoPolin(x, y, opt.N);                                                                     //Получение полинома методом СЛАУ
    if(opt.OnFunc){                                                                                             //Решение через функции Лагранжа и Ньютона
        Double_t *par = new Double_t[2*opt.N+1];
        par[0] = opt.N-1;
        for(size_t i = 0; i<opt.N; i++){                                                                        //Объединение N и всех x и y в один массив параметров
            par[i+1] = x[i];
            par[opt.N+1+i] = y[i];
        }

        funcMass[2] = new TF1(("LagranPol"+to_string(count)).c_str(), LagranFunc, 0., 1., 2*opt.N+1);           //Получение полинома методом Лагранжа
        funcMass[2]->SetParameters(par);
        T *f = FNew(x, y, opt.N);                                                                               //Получение разделённых разностей
        for(size_t i = 0; i<opt.N; i++){par[i+opt.N+1] = f[i];}
        funcMass[3] = new TF1(("NewtonPol"+to_string(count)).c_str(), NewtonFunc, 0., 1., 2*opt.N+1);           //Получение полинома методом Ньютона
        funcMass[3]->SetParameters(par);
    }
    else{                                                                                                       //Получение через коэффициенты полиномов Лагранжа и Ньютона
        funcMass[2] = LagranPolin(x, y, opt.N, count);
        funcMass[3] = NewtonPolin(x, y, opt.N, count);
    }

    return funcMass;
}

void SettingUpFunc(TF1** funcMass, TLegend* leg, int k, struct FuncOptions &opt){
    for(size_t i = 0; i<k; i++){
        if(funcMass[i]){
            funcMass[i]->SetLineColor(i+1);                                                                     //Настройка цвета линий функций
            funcMass[i]->SetLineStyle(11-i);                                                                    //Настройка стиля линий функций
            funcMass[i]->SetRange(opt.xmin-opt.dl-1, opt.xmax+opt.dr+1);                                        //Настройка границ отрисовки функций
            if(leg){leg->AddEntry(funcMass[i], opt.names[i].c_str(), "l");}                                     //Добавление функций в легенду
        }
    }
}

template <typename T> TF1* ProstoPolin(T* x, T* y, int N){
    TMatrixD X(N, N);                                                                                           //Объявление необходимых объектов и переменных
    TVectorD a(N);
    TVectorD f(N);

    for(size_t i = 0; i<N; i++){                                                                                //Заполнение матрицы и вектора
        f(i) = y[i];
        X(i, 0) = 1;
        for(size_t j = 1; j<N; j++){
            X(i, j) = X(i, j-1)*x[i];
        }
    }

    if(X.Determinant()){                                                                                        //Проверка на некорректный или нулевой определитель
        a = X.Invert()*f;                                                                                       //Нахождение коэффициентов полинома
    }
    else{
        cout<<"Обратной матрицы не существует. Полином не найден методом СЛАУ"<<endl;
        return 0;
    }

    TF1 *f_pol = new TF1("ProstoPol", ("pol"+to_string(N-1)+"(0)").c_str());                                    //Создание полинома
    for(size_t i = 0; i<N; i++){
        f_pol->SetParameter(i, a(i));
    }
    return f_pol;
}

template <typename T> TF1* LagranPolin(T* x, T* y, int N, int count){
    T *v1 = new T[N], *v2;                                                                                      //Объявление необходимых объектов и переменных
    fill_n(v1, N, 0);                                                                                           //Обнуление массива, который будет хранить параметры полинома

    for(size_t i = 0; i<N; i++){
        v2 = PolLag(x, i, N-1);                                                                                 //Получение параметров малого полинома Лагранжа
        if(v2){SummMassive(v1, N, v2, N, 1.0, y[i]);}                                                           //Сложение параметров малых полиномов
        else{
            cout<<"Полином не найден методом Лагранжа"<<endl;
            return 0;
        }
    }

    TF1 *f_pol = new TF1(("LagranPol"+to_string(count)).c_str(), ("pol"+to_string(N-1)+"(0)").c_str());         //Создание полинома
    for(size_t i = 0; i<N; i++){
        f_pol->SetParameter(i, v1[i]);
        //cout<<v1[i]<<endl;
    }
    //cout<<endl;
    //a.Print();

    return f_pol;
}

template <typename T> TF1* NewtonPolin(T* x, T* y, int N, int count){
    T *v1 = new T[N], *v2;                                                                                      //Объявление необходимых объектов и переменных
    fill_n(v1, N, 0);                                                                                           //Обнуление массива, который будет хранить параметры полинома

    T *f = FNew(x, y, N);                                                                                       //Получение разделённых разностей

    if(f){
        for(size_t i = 1; i<N; i++){
            v2 = PolNew(x, i);                                                                                  //Получение параметров малого полинома Лагранжа
            SummMassive(v1, N, v2, i+1, 1.0, f[i]);                                                             //Сложение параметров малых полиномов
        }
    }
    else{
        cout<<"Полином не найден методом Ньютона"<<endl;
        return 0;
    }

    TF1 *f_pol = new TF1(("NewtonPol"+to_string(count)).c_str(), ("pol"+to_string(N-1)+"(0)").c_str());         //Создание полинома
    for(size_t i = 0; i<N; i++){
        f_pol->SetParameter(i, v1[i]);
    }

    return f_pol;
}

template <typename T>
void SummMassive(T* arr1, int N1, T* arr2, int N2, long double k1, long double k2){
    if(N1>=N2){
        for(size_t i = 0; i<N2; i++){
            arr1[i]=k1*arr1[i]+k2*arr2[i];
        }
    }
}

template <typename T> T* PolLag(T* x, int j, int N){
    T* v = new T[N+1];                                                                                          //Объявление необходимых объектов и переменных
    T dx;
    fill_n(v, N+1, 0);                                                                                          //Обнуление массива, который будет хранить параметры полинома
    v[0] = 1;
    int N0 = 0;
    for(size_t i = 0; i<=N; i++){
        if(i!=j){
            dx = x[j]-x[i];
            if(dx){
                for(size_t k = N0; k>=0; k--){
                    v[k+1]+=v[k]/dx;
                    v[k]*=-x[i]/dx;
                }
                N0++;
            }
            else{return 0;}
        }
    }
    return v;
}

template <typename T> T* PolNew(T* x, int j){
    T* v = new T[j+1];                                                                                          //Объявление необходимых объектов и переменных
    fill_n(v, j+1, 0);                                                                                          //Обнуление массива, который будет хранить параметры полинома
    v[0] = 1;
    for(size_t i = 0; i<j; i++){
        for(size_t k = i; k>=0; k--){
            v[k+1]+=v[k];
            v[k]*=-x[i];
        }
    }
    return v;
}

template <typename T> T* FNew(T* x, T* y, int N){
    T *f = new T[N]; copy(y, y+N, f);                                                                           //Объявление необходимых объектов и переменных
    T dx;
    for (size_t N0 = 0; N0<N; N0++){
        for(size_t i = N-1; i>N0; i--){
            dx=x[i]-x[i-1-N0];
            if(dx){f[i]=(f[i]-f[i-1])/dx;}
            else{return 0;}
        }
    }
    return f;
}

Double_t LagranFunc(Double_t *x, Double_t *par){
    Double_t res = 0.0;                                                                                         //Объявление необходимых объектов и переменных
    Double_t fl;
    Int_t N = int(par[0]);
    for(size_t j = 0; j<N+1; j++){
        fl = par[j+N+2];
        for(size_t i = 0; i<N+1; i++){
            if(i!=j){fl*=(x[0]-par[i+1])/(par[j+1]-par[i+1]);}
        }
        res+=fl;
    }
    return res;
}

Double_t NewtonFunc(Double_t *x, Double_t *par){
    Double_t res = 0.0;                                                                                         //Объявление необходимых объектов и переменных
    Int_t N = int(par[0]);
    // cout<<x[0]<<endl;
    //Double_t *f = FNew(par+1, par+N+2, N+1);
    Double_t fl;
    for(size_t j = 0; j<N+1; j++){
        //fl = f[j];
        fl = par[j+N+2];
        for(size_t i = 0; i<j; i++){
            fl*=x[0]-par[i+1];
        }
        res+=fl;
    }
    return res;
}

template <typename T> void Part2(struct FuncOptions &opt, T type){
    double lg = 0.05, rg = 0.05;

////////////////////////////////7 Point/////////////////////////////////////////////
    auto points1 = GetPointsUniform(type, opt);                                                                 //Генерация новых точек по равномерной сетке
    auto x1 = points1[0], y1 = points1[1];                                                                      //Разделение координат на два массива
    auto *x_gr1 = SomeToDobleMass(x1, opt.N), *y_gr1 = SomeToDobleMass(y1, opt.N);

    TF1 **funcMass1 = Get4Functions(&x1[0], &y1[0], opt, true);                                                 //Перезагрузка полиномов

    TGraph* gr1 = new TGraph(opt.N, &x_gr1[0], &y_gr1[0]);                                                      //Объявление графика с новыми точками
    gr1->SetTitle(("Interpolation of "+to_string(opt.N-1)+" and "+to_string(opt.N)+" points in the range ["+to_string(int(opt.xmin))+":"+to_string(int(opt.xmax))+"]").c_str());
    gr1->GetXaxis()->SetLimits(opt.xmin-lg, opt.xmax+rg);

    TF1 *fDivReal1 = RealDiff(&x1[0], &y1[0], opt, true);                                                       //Фактическая разница интерполяции
    //new TF1(("DivReal"+to_string(opt.N)).c_str(), ("abs(f0-NewtonPol"+to_string(opt.N)+")").c_str());
    fDivReal1->SetRange(opt.xmin-0.05, opt.xmax+0.05);
    fDivReal1->SetLineColor(2);
    fDivReal1->SetLineStyle(11);

    TF1 *fDivMax1 = PolinDiv(&x1[0], &y1[0], opt, true);                                                        //Оценочная разница интерполяции
    fDivMax1->SetRange(opt.xmin-lg, opt.xmax+rg);
    fDivMax1->SetLineColor(2);
    fDivMax1->SetLineStyle(10);

////////////////////////////////8 Point/////////////////////////////////////////////
    opt.N++;                                                                                                    //Увеличение количества точек
    auto points2 = GetPointsUniform(type, opt);                                                                 //Генерация новых точек
    auto x2 = points2[0], y2 = points2[1];                                                                      //Разделение координат на два массива
    auto *x_gr2 = SomeToDobleMass(x2, opt.N), *y_gr2 = SomeToDobleMass(y2, opt.N);

    TF1 **funcMass2 = Get4Functions(&x2[0], &y2[0], opt);                                                       //Перезагрузка полиномов

    TGraph* gr2 = new TGraph(opt.N, &x_gr2[0], &y_gr2[0]);                                                      //Объявление графика с новыми точками
    gr2->SetTitle(("Interpolation of "+to_string(opt.N-1)+" and "+to_string(opt.N)+" points in the range ["+to_string(int(opt.xmin))+":"+to_string(int(opt.xmax))+"]").c_str());
    gr2->GetXaxis()->SetLimits(opt.xmin-lg, opt.xmax+rg);

    TF1 *fDivReal2 = RealDiff(&x2[0], &y2[0], opt);                                                             //Фактическая разница интерполяции
    fDivReal2->SetRange(opt.xmin-lg, opt.xmax+rg);
    fDivReal2->SetLineColor(4);
    fDivReal2->SetLineStyle(11);

    TF1 *fDivMax2 = PolinDiv(&x2[0], &y2[0], opt);                                                              //Оценочная разница интерполяции
    fDivMax2->SetRange(opt.xmin-lg, opt.xmax+rg);
    fDivMax2->SetLineColor(4);
    fDivMax2->SetLineStyle(10);
    opt.N--;

    TLegend* leg1 = new TLegend(0.3, 0.5, 0.7, 0.8);                                                            //Создание легенды
    leg1->SetBorderSize(0);
    leg1->AddEntry(fDivReal1, ("Real diff for "+to_string(opt.N)).c_str(), "l");
    leg1->AddEntry(fDivMax1, ("Estimated diff for "+to_string(opt.N)).c_str(), "l");
    leg1->AddEntry(fDivReal2, ("Real diff for "+to_string(opt.N+1)).c_str(), "l");
    leg1->AddEntry(fDivMax2, ("Estimated diff for "+to_string(opt.N+1)).c_str(), "l");

    TCanvas* c = new TCanvas();                                                                                 //Отрисовка графиков и функций
    c->Divide(1, 2, 0.001, 0.001);
    c->cd(1);
    gPad->SetTicks();
    gPad->SetGrid();
    fDivMax1->SetTitle(("R(x) of "+to_string(opt.N)+" and "+to_string(opt.N+1)
                                  +" points in the range ["+to_string(int(opt.xmin))
                                  +":"+to_string(int(opt.xmax))+"]").c_str());
    fDivMax1->Draw();
    fDivReal1->Draw("same");
    fDivMax2->Draw("same");
    fDivReal2->Draw("same");
    leg1->Draw("same");

    TLegend* leg2 = new TLegend(0.65, 0.5, 0.85, 0.8);                                                          //Создание легенды
    leg2->SetBorderSize(0);

    funcMass1[3]->SetLineColor(2);
    funcMass1[3]->SetLineStyle(10);
    funcMass1[3]->SetRange(opt.xmin-opt.dl-1, opt.xmax+opt.dr+1);
    SettingUpFunc(funcMass2, 0, 4, opt);

    TF1 *fClone = (TF1*)funcMass1[3]->Clone();

    leg2->AddEntry(fClone, ("Newton for "+to_string(opt.N)+" point").c_str(), "l");
    leg2->AddEntry(funcMass2[3], ("Newton for "+to_string(opt.N+1)+" point").c_str(), "l");

    gr1->SetMarkerColor(2);
    gr2->SetMarkerColor(4);

    c->cd(2);                                                                                                   //Отрисовка графиков и функций
    gPad->SetTicks();
    gPad->SetGrid();
    gr2->Draw("A*");
    gr1->Draw("*same");
    fClone->Draw("same");
    funcMass2[3]->Draw("same");
    //funcMass1[0]->Draw("same");
    leg2->Draw("same");
}

TF1* GetDeriveNforExpSinOrCos(double ksin, double kcos, int N){
    for(size_t i = 0; i<N; i++){
        ksin = -(ksin+kcos);
        kcos = -(ksin+2*kcos);
    }
    //cout<<"ksin = "<<ksin<<endl;
    //cout<<"kcos = "<<kcos<<endl;

    TF1 *f = new TF1("DeriveN", "abs(exp(-x)*([0]*sin(x)+[1]*cos(x)))");                                        //Создание функции N-ой производной
    f->SetParameter(0, ksin);   f->SetParameter(1, kcos);

    return f;
}

template <typename T> TF1* RealDiff(T* x, T* y, struct FuncOptions &opt, bool key, int dcount){
    static int count = 0;
    cout<<count;
    if(key){count = 0;}
    else{count += dcount;}
    cout<<" "<<count<<" "<<key<<endl;

    TF1 *func;                                                                                                  //Объявление необходимых объектов и переменных
    int N = opt.func.size();
    string Name = opt.names[opt.typePol-1];

    if(opt.OnFunc){                                                                                             //Решение через функции Лагранжа и Ньютона
        Double_t *par = new Double_t[N+1+2*opt.N+1];
        par[0] = N; par[N+1] = opt.N-1;

        for(size_t i = 0; i<N; i++){
            par[i+1] = opt.func[i];
        }
        for(size_t i = 1; i<2*opt.N+1; i++){
            if(i>opt.N){par[i+N+1] = y[i-opt.N-1];}
            else{par[i+N+1] = x[i-1];}
        }

        // for(size_t i = 0; i<N+1+opt.N+1; i++){
        //     if(i>0 && i<N+1){cout<<i<<" "<<char(int(par[i]))<<endl;}
        //     else{cout<<i<<" "<<par[i]<<endl;}
        // }
        if(Name=="Newton"){
            T *f = FNew(x, y, opt.N);                                                                           //Получение разделённых разностей
            for(size_t i = 0; i<opt.N; i++){par[i+N+1+opt.N+1] = f[i];}                                         //Объединение N и всех x и y в один массив параметров
            func = new TF1(("DivReal"+to_string(count)).c_str(), NewtonDiv, 0., 1., N+1+2*opt.N+1);
        }
        else if(Name=="Lagran"){
            for(size_t i = 0; i<opt.N; i++){par[i+N+1+opt.N+1] = y[i];}                                         //Объединение N и всех x и y в один массив параметров
            func = new TF1(("DivReal"+to_string(count)).c_str(), LagranDiv, 0., 1., N+1+2*opt.N+1);
        }
        else{return 0;}

        func->SetParameters(par);
    }
    else{                                                                                                       //Получение через коэффициенты полиномов Лагранжа и Ньютона
        func = new TF1(("DivReal"+to_string(count)).c_str(), ("abs(f0"+to_string(count)+"-"+Name+"Pol"+to_string(count)+")").c_str());
    }

    return func;
}

Double_t LagranDiv(Double_t *x, Double_t *par){
    Double_t res;                                                                                               //Объявление необходимых объектов и переменных

    Int_t N = int(par[0]);
    string str = "";

    for(size_t i = 0; i<N; i++){
        str+=char(par[i+1]);
    }

    TF1 *func = new TF1("", str.c_str());
    res = fabs(func->Eval(x[0])-LagranFunc(x, par+N+1));
    return res;
}

Double_t NewtonDiv(Double_t *x, Double_t *par){
    Double_t res;                                                                                               //Объявление необходимых объектов и переменных
    Int_t N = int(par[0]);
    string str = "";

    for(size_t i = 0; i<N; i++){
        str+=char(par[i+1]);
    }

    TF1 *func = new TF1("", str.c_str());
    res = fabs(func->Eval(x[0])-NewtonFunc(x, par+N+1));
    return res;
}

template <typename T> TF1* PolinDiv(T* x, T* y, struct FuncOptions &opt, bool key, int dcount){
    static int count = 0;
    if(key){count = 0;}
    else{count += dcount;}

    TF1 *f_pol;

    long double f = EstimMax(opt);                                                                              //Получение максимального модуля N-ой производной
    if(!f){
        cout<<"Полином не найден методом Ньютона"<<endl;
        return 0;
    }
    //cout<<f<<endl;

    if(opt.OnFunc){
        f_pol = new TF1(("PolinDiv"+to_string(count)).c_str(), FuncDiv, 0., 1., 3+opt.N);
        f_pol->SetParameter(0, 1);  f_pol->SetParameter(1, f); f_pol->SetParameter(2, opt.N);
        //cout<<f<<endl;
        for(size_t i = 0; i<opt.N; i++){
            f_pol->SetParameter(i+3, x[i]);
            //cout<<i<<" "<<x[i]<<endl;
        }
    }
    else{
        T *v2;                                                                                                  //Объявление необходимых объектов и переменных
        v2 = PolNew(x, opt.N);                                                                                  //Получение параметров малого полинома Лагранжа

        f_pol = new TF1(("PolinDiv"+to_string(count)).c_str(), ("[0]*abs(pol"+to_string(opt.N)+"(1))").c_str());                      //Создание полинома
        f_pol->SetParameter(0, f);
        for(size_t i = 0; i<opt.N+1; i++){
            f_pol->SetParameter(i+1, v2[i]);
            //cout<<v2[i]<<endl;
        }
    }

    return f_pol;
}

Double_t FuncDiv(Double_t *x, Double_t *par){
    Double_t res = 1.;                                                                                          //Объявление необходимых объектов и переменных

    Int_t Nf = int(par[0]);
    for(size_t i = 1; i<Nf+1; i++){res *= par[i];}

    Int_t Nx = int(par[Nf+1]);
    for(size_t i = Nf+2; i<Nf+Nx+2; i++){res *= x[0] - par[i];}

    return fabs(res);
}

long double EstimMax(struct FuncOptions opt){
    long double f;                                                                                              //Объявление необходимых объектов и переменных
    if(opt.func=="exp(-x)*sin(x)"){                                                                             //Получение максимального значения N-ой производной
        f = GetDeriveNforExpSinOrCos(1, 0, opt.N)->GetMaximum(opt.xmin, opt.xmax);
    }
    else if(opt.func=="exp(-x)*cos(x)"){
        f = GetDeriveNforExpSinOrCos(0, 1, opt.N)->GetMaximum(opt.xmin, opt.xmax);
    }
    else{return 0;}

    for(size_t i = 1; i<opt.N; i++) {                                                                           //Деление на N!
        //cout<<"f["<<i-1<<"] = "<<f<<endl;
        f /= i+1;
        //cout<<i+1<<endl;
        //cout<<"f_res = "<<f<<endl;
    }
    return abs(f);
}

template <typename T> void Part3(struct FuncOptions &opt, T type){
    double lg = 0.05, rg = 0.05;

////////////////////////////////Uniform/////////////////////////////////////////////
    auto points1 = GetPointsUniform(type, opt);                                                                 //Генерация новых точек
    auto x1 = points1[0], y1 = points1[1];                                                                      //Разделение координат на два массива
    auto *x_gr1 = SomeToDobleMass(x1, opt.N), *y_gr1 = SomeToDobleMass(y1, opt.N);

    TF1 **funcMass1 = Get4Functions(&x1[0], &y1[0], opt, true);                                                 //Перезагрузка полиномов

    TGraph* gr1 = new TGraph(opt.N, &x_gr1[0], &y_gr1[0]);                                                      //Объявление графика с новыми точками
    gr1->SetTitle(("Interpolation of "+to_string(opt.N)+" points in the range ["+to_string(int(opt.xmin))+":"+to_string(int(opt.xmax))+"]").c_str());
    gr1->GetXaxis()->SetLimits(opt.xmin-lg, opt.xmax+rg);

    TF1 *fDivReal1 = RealDiff(&x1[0], &y1[0], opt, true);                                                       //Фактическая разница интерполяции
    fDivReal1->SetRange(opt.xmin-lg, opt.xmax+rg);
    fDivReal1->SetLineColor(2);
    fDivReal1->SetLineStyle(11);

    TF1 *fDivMax1 = PolinDiv(&x1[0], &y1[0], opt, true);                                                        //Оценочная разница интерполяции
    fDivMax1->SetRange(opt.xmin-lg, opt.xmax+rg);
    fDivMax1->SetLineColor(2);
    fDivMax1->SetLineStyle(10);

////////////////////////////////Chebyshev///////////////////////////////////////////
    auto points2 = GetPointsChebyshev(type, opt);                                                               //Генерация новых точек
    auto x2 = points2[0], y2 = points2[1];                                                                      //Разделение координат на два массива
    auto *x_gr2 = SomeToDobleMass(x2, opt.N), *y_gr2 = SomeToDobleMass(y2, opt.N);

    TF1 **funcMass2 = Get4Functions(&x2[0], &y2[0], opt);                                                       //Перезагрузка полиномов

    TGraph* gr2 = new TGraph(opt.N, &x_gr2[0], &y_gr2[0]);                                                      //Объявление графика с новыми точками
    gr2->SetTitle(("Interpolation of "+to_string(opt.N)+" points in the range ["+to_string(int(opt.xmin))+":"+to_string(int(opt.xmax))+"]").c_str());
    gr2->GetXaxis()->SetLimits(opt.xmin-lg, opt.xmax+rg);

    TF1 *fDivReal2 = RealDiff(&x2[0], &y2[0], opt);                                                             //Фактическая разница интерполяции
    fDivReal2->SetRange(opt.xmin-lg, opt.xmax+rg);
    fDivReal2->SetLineColor(4);
    fDivReal2->SetLineStyle(11);

    TF1 *fDivMax2 = PolinDiv(&x2[0], &y2[0], opt);                                                              //Оценочная разница интерполяции
    fDivMax2->SetRange(opt.xmin-lg, opt.xmax+rg);
    fDivMax2->SetLineColor(4);
    fDivMax2->SetLineStyle(10);

    SettingUpFunc(funcMass2, 0, 4, opt);                                                                        //Настройка параметров отрисовки
    gr1->SetMarkerColor(2);
    gr2->SetMarkerColor(4);
    fDivMax1->SetTitle(("R(x) for Chebyshev and Uniform points in the range ["
                            +to_string(int(opt.xmin))+":"
                            +to_string(int(opt.xmax))+"]").c_str());
    funcMass1[3]->SetLineColor(2);
    funcMass1[3]->SetLineStyle(10);
    funcMass1[3]->SetRange(opt.xmin-opt.dl-1, opt.xmax+opt.dr+1);

    TLegend* leg1 = new TLegend(0.3, 0.5, 0.7, 0.8);                                                            //Создание легенды
    leg1->SetBorderSize(0);
    leg1->AddEntry(fDivReal1, ("Real diff for "+to_string(opt.N)+" point of Uniform").c_str(), "l");
    leg1->AddEntry(fDivMax1, ("Estimated diff for "+to_string(opt.N)+" point of Uniform").c_str(), "l");
    leg1->AddEntry(fDivReal2, ("Real diff for "+to_string(opt.N)+" point of Chebyshev").c_str(), "l");
    leg1->AddEntry(fDivMax2, ("Estimated diff for "+to_string(opt.N)+" point of Chebyshev").c_str(), "l");

    TLegend* leg2 = new TLegend(0.65, 0.5, 0.85, 0.8);                                                          //Создание легенды
    leg2->SetBorderSize(0);
    leg2->AddEntry(funcMass1[3], ("Newton for "+to_string(opt.N)+" point of Uniform").c_str(), "l");
    leg2->AddEntry(funcMass2[3], ("Newton for "+to_string(opt.N)+" point of Chebyshev").c_str(), "l");

    TCanvas* c = new TCanvas();                                                                                 //Отрисовка графиков и функций
    c->Divide(1, 2, 0.001, 0.001);

    c->cd(1);
    gPad->SetTicks();
    gPad->SetGrid();
    fDivMax1->Draw();
    fDivReal1->Draw("same");
    fDivMax2->Draw("same");
    fDivReal2->Draw("same");
    leg1->Draw("same");

    c->cd(2);                                                                                                   //Отрисовка графиков и функций
    gPad->SetTicks();
    gPad->SetGrid();
    gr1->Draw("A*");
    gr2->Draw("*same");
    funcMass1[3]->Draw("same");
    funcMass2[3]->Draw("same");
    leg2->Draw("same");
}

template <typename T> T** GetPointsChebyshev(T flag, struct FuncOptions &opt){
    TF1 *f = new TF1("f", opt.func.c_str());                                                                    //Объявление необходимых объектов и переменных
    T *x = new T[opt.N], *y = new T[opt.N];
    T **res = new T*[2];

    for(size_t i = 0; i<opt.N; i++){                                                                            //Заполнение векторов координатами точек
        x[i] = (opt.xmax+opt.xmin)/2. + (opt.xmax-opt.xmin)/2.*cos((2*i+1)*TMath::Pi()/opt.N/2);
        //cout<<x[i]<<endl;
        y[i] = f->Eval(x[i]);
    }
    res[0] = x;
    res[1] = y;
    return res;
}