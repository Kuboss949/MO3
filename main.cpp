#include <iostream>
#include <math.h>
#include <functional>

#define ERR 1.3e-6
#define FUNERR 1.3e-6


using namespace std;


void picard(function<double(double)>fx, double initialGuess, int MaxIteracji){
    double xn=initialGuess, xn1;
    for(int i=0; i<MaxIteracji; i++){
        xn1=fx(xn)+xn;
        if(i+1==MaxIteracji){
            cout << "Zakończono z powodu wyczerpania liczby iteracji";
            return;
        }
        cout << "Iteracja nr " << i+1 << " xn+1= " << xn1 << ". Estymator bledu: " << xn1-xn << ". Residuum rownania: "<< fx(xn1) <<  endl;
        if(abs(xn1-xn)<=ERR && abs(fx(xn))<=FUNERR){
            cout << "Przerwane ze wzgledu na kryterium dokladnosci wyznaczenia xn oraz na kryterium wiarygodnosci xn jako przyblizenia pierwiastka" << endl;
            break;
        }
        xn=xn1;
    }
    cout << "Pierwiastek: " << xn1 << endl;
}

void bisekcja(function<double(double)>fx, double a, double b, int MaxIteracji){
    double x;
    if(fx(a)<0 && fx(b)<0 || fx(a)>0 && fx(b)>0){
        cout << "Brak pierwiastków w podanym przedziale"<< endl;
    }
    if(a>b){
        swap(a,b);
    }
    for(int i=0; i<MaxIteracji; i++){
        x=(a+b)/2.0;
        if(i+1==MaxIteracji){
            cout << "Zakonczono z powodu wyczerpania liczby iteracji"<< endl;
            return;
        }
        cout << "Iteracja nr " << i+1 << "  Przedzial: ["<< a << "," << b << "]" << " xn= " << x  << ". Estymator bledu: " << (b-a)/2 << ". Residuum rownania: "<< fx(x) << endl;
        if(abs((b-a)/2)<=ERR && abs(fx(x))<=FUNERR){
            cout << "Przerwane ze wzgledu na kryterium dokladnosci wyznaczenia xn oraz na kryterium wiarygodnosci xn jako przyblizenia pierwiastka" << endl;
            break;
        }

        if(x<0 && b>0 || x>0 && b<0){
            a=x;
        }else{
            b=x;
        }
    }
    cout << "Pierwiastek: " << x << endl;
}

void newton(function<double(double)>fx, function<double(double)>derivative, double initialGuess, int MaxIteracji){
    double xn=initialGuess, xn1;
    for(int i=0; i<MaxIteracji; i++){
        if(derivative(xn)==0){
            cout << "Dzielenie przez 0! Algorytm konczy dzialanie!";
            return;
        }
        xn1=xn - fx(xn)/derivative(xn);
        if(i+1==MaxIteracji){
            cout << "Zakonczono z powodu wyczerpania liczby iteracji"<< endl;
            return;
        }
        cout << "Iteracja nr " << i+1 << " xn+1= " << xn1 << ". Estymator bledu: " << xn1-xn << ". Residuum rownania: "<< fx(xn1) <<  endl;

        if(abs(xn1-xn)<=ERR && abs(fx(xn))<=FUNERR){
            cout << "Przerwane ze wzgledu na kryterium dokladnosci wyznaczenia xn oraz na kryterium wiarygodnosci xn jako przyblizenia pierwiastka" << endl;
            break;
        }

        xn=xn1;
    }
    cout << "Pierwiastek: " << xn1 << endl;
}

void sieczne(function<double(double)>fx, double initialGuess1, double initialGuess2, int MaxIteracji){
    double xn=initialGuess1, xn1 = initialGuess2, xn2;
    for(int i=0; i<MaxIteracji; i++){
        if(((fx(xn1)-fx(xn))/(xn1-xn))==0){
            cout << "Dzielenie przez 0! Algorytm konczy dzialanie!";
            return;
        }

        xn2=xn1 - fx(xn1)/((fx(xn1)-fx(xn))/(xn1-xn));
        if(i+1==MaxIteracji){
            cout << "Zakonczono z powodu wyczerpania liczby iteracji"<< endl;
            return;
        }
        cout << "Iteracja nr " << i+1 << " xn+1= " << xn2 << ". Estymator bledu: " << xn2-xn1 << ". Residuum rownania: "<< fx(xn2) <<  endl;
        if(abs(xn2-xn1)<=ERR && abs(fx(xn2))<=FUNERR){
            cout << "Przerwane ze wzgledu na kryterium dokladnosci wyznaczenia xn oraz na kryterium wiarygodnosci xn jako przyblizenia pierwiastka" << endl;
            break;
        }
        xn=xn1;
        xn1=xn2;
    }
    cout << "Pierwiastek: " << xn2 << endl;
}



int main() {

    cout << "Picard:" << endl;

    picard([](double x)->double{return sin(x/4.0)*sin(x/4.0)-x;}, 323.0, 1000);

    cout << "--------------------------------------------------" << endl;

    picard([](double x)->double{return tan(2.0*x)-x-1.0;}, 0.0, 6);

    cout << "--------------------------------------------------" << endl << endl << endl;




    cout << "Bisekcja:"<< endl;
    bisekcja([](double x)->double{return sin(x/4.0)*sin(x/4.0)-x;}, -4.0, 2.0, 1000);
    cout << "--------------------------------------------------" << endl;
    bisekcja([](double x)->double{return tan(2.0*x)-x-1.0;}, -1.0, 2.0, 1000);
    cout << "--------------------------------------------------" << endl << endl << endl;



    cout << "Newton:" << endl;
    newton([](double x)->double{return sin(x/4.0)*sin(x/4.0)-x;},[](double x)->double{return 0.25*sin(x/2.0)-1.0;}, 323.0, 1000);
    cout << "--------------------------------------------------" << endl;
    newton([](double x)->double{return tan(2.0*x)-x-1.0;},[](double x)->double{return -1.0+2.0/(cos(2.0*x)*cos(2.0*x));}, 2.0, 1000);
    cout << "--------------------------------------------------" << endl << endl << endl;

    cout << "Sieczne:" << endl;
    sieczne([](double x)->double{return sin(x/4.0)*sin(x/4.0)-x;}, 23.0, 12.0, 1000);
    cout << "--------------------------------------------------" << endl;
    sieczne([](double x)->double{return tan(2.0*x)-x-1.0;}, 0.0, 1.0,  1000);
    cout << "--------------------------------------------------" << endl << endl << endl;

    return 0;
}
