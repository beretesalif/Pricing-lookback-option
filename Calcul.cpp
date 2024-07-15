#include "pch.h"
#include<iostream>
#include <vector>
#include <cmath>
#include<algorithm>
#include <cstdlib>  // Pour rand()
#include <ctime>    // Pour srand()
#include <cstring>
#include "Calcul.h"
#include "PathSimulator.h"
#include "LookBack.h"



int __stdcall MonteCarloSimulation(double S_t, double r, double T, double sigma, int NbrePas, int NmbreSim, double eps, double proba, const char* optionType,
    double* ptrPrice, double* ptrPriceApprox, double* ptrDelta, double* ptrGamma, double* ptrVega, double* ptrTheta, double* ptrRho, double* ptrbornInf, double* ptrbornSup) {

    // Utilisez la fonction rand() de <cstdlib> avec srand() pour initialiser la séquence de nombres aléatoires.
    // La fonction rand() n'est pas idéale pour les simulations Monte Carlo avancées, mais nous la garderons pour cet exemple.
    srand(static_cast<unsigned>(time(0)));

    PathSimulator path(S_t, r, T, sigma, NbrePas);

    LookBack look(optionType, eps);

    double sumLookbackPayoff = 0.0;
    double deltaSum = 0.0;
    double gammaSum = 0.0;
    double vegaSum = 0.0;
    double RhoSum = 0.0;
    double ThetaSum = 0.0;
    double sumVaR = 0.0;

    double MinMoyen = 0.0, MaxMoyen = 0.0, S_finalMoyen = 0.0;

    for (int i = 0; i < NmbreSim; ++i) {
        path.simulatePath();
        sumLookbackPayoff += look.Payoff(path);
        sumVaR += look.Payoff(path)*look.Payoff(path);

        deltaSum += look.Delta(path);
        gammaSum += look.Gamma(path);
        vegaSum += look.Vega(path);
        RhoSum += look.Rho(path);
        ThetaSum += look.Theta(path);

        MinMoyen += path.getMin();
        MaxMoyen += path.getMax();
    }
    //path.simulatePath();

    // Calcul et retourne le prix actualisé de l'option en fonction du type spécifié
    *ptrPrice = std::exp(-r * T) * (sumLookbackPayoff / NmbreSim);
    *ptrDelta = std::exp(-r * T) * (deltaSum / NmbreSim);
    *ptrGamma = std::exp(-r * T) * (gammaSum / NmbreSim);
    *ptrVega = std::exp(-r * T) * (vegaSum / NmbreSim);
    *ptrRho = RhoSum / NmbreSim;
    *ptrTheta = -ThetaSum / NmbreSim; 


    double a = MinMoyen/NmbreSim, b = MaxMoyen /NmbreSim ;
    double var = (NmbreSim / (NmbreSim - 1)) * (
        (1 / NmbreSim) * sumVaR * std::exp(-r * 2 * T) -
        (1 / NmbreSim) * sumLookbackPayoff * std::exp(-r * 2 * T) * (1 / NmbreSim) * sumLookbackPayoff);

    double pvalue = NormSInv((1+ proba)/2);

    *ptrbornInf = std::exp(-r * T) * ((sumLookbackPayoff / NmbreSim) - (pvalue * sigma) / sqrt(NmbreSim));
    *ptrbornSup = std::exp(-r * T) * ((sumLookbackPayoff / NmbreSim) + (pvalue * sigma) / sqrt(NmbreSim));

    *ptrPriceApprox = look.PayoffApproxim(path);

    return 0; // Return an integer, not a double
}