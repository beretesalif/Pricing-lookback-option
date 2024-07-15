#pragma once
#include <cmath>
#include <ctime>
#include <random>
#include<iostream>
#include <cstring>



#ifdef CALCUL_EXPORTS
#define CALCUL_API __declspec(dllexport)
#else
#define CALCUL_API __declspec(dllimport)
#endif

extern "C" CALCUL_API int __stdcall MonteCarloSimulation(double S_t, double r, double T, double sigma, int NbrePas, int NmbreSim, double eps, double proba, const char* optionType,
	double* ptrPrice, double* ptrPriceApprox, double* ptrDelta, double* ptrGamma, double* ptrVega, double* ptrTheta, double* ptrRho,double* ptrbornInf, double* ptrbornSup);

double NormSInv(double p) {
    double a1 = -39.6968302866538, a2 = 220.946098424521, a3 = -275.928510446969;
    double a4 = 138.357751867269, a5 = -30.6647980661472, a6 = 2.50662827745924;
    double b1 = -54.4760987982241, b2 = 161.585836858041, b3 = -155.698979859887;
    double b4 = 66.8013118877197, b5 = -13.2806815528857, c1 = -7.78489400243029E-03;
    double c2 = -0.322396458041136, c3 = -2.40075827716184, c4 = -2.54973253934373;
    double c5 = 4.37466414146497, c6 = 2.93816398269878, d1 = 7.78469570904146E-03;
    double d2 = 0.32246712907004, d3 = 2.445134137143, d4 = 3.75440866190742;
    double p_low = 0.02425, p_high = 1 - p_low;
    double q, r;
    double retVal;

    if ((p < 0) || (p > 1)) {
        std::cerr << "NormSInv: Argument out of range." << std::endl;
        retVal = 0;
    }
    else if (p < p_low) {
        q = std::sqrt(-2 * std::log(p));
        retVal = (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
    }
    else if (p <= p_high) {
        q = p - 0.5;
        r = q * q;
        retVal = (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q / (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1);
    }
    else {
        q = std::sqrt(-2 * std::log(1 - p));
        retVal = -(((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6) / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1);
    }

    return retVal;
}