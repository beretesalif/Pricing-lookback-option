#pragma once
// LookBack.h
#ifndef LOOKBACK_H
#define LOOKBACK_H

#define PI 3.14159265358979
#define MAX(A,B) ( (A) > (B) ? (A):(B) )
#define MIN(A,B) ( (A) < (B) ? (A):(B) )
#define SQR(X) ((X)*(X))

#include "PathSimulator.h"

class LookBack {
public:
    // Constructeur paramétrique
    LookBack(const char* optionType, double epsilon);

    // Méthodes de calcul
    double Payoff(const PathSimulator& pathSimulator) const;
    double Delta(PathSimulator pathSimulator) const;
    double Gamma(PathSimulator pathSimulator) const;
    double Vega(PathSimulator pathSimulator) const;
    double Theta(PathSimulator pathSimulator) const;
    double Rho(PathSimulator pathSimulator) const;

    // Nouvelle méthode avec des paramètres supplémentaires
    double PayoffApproxim(PathSimulator pathSimulator) ;
    // Attributs
private: 
    const char* optionType_;
    double epsilon_;
    double N(double x) const;
    double a_1(const double S,    // Spot price
        const double H,    // Min/max of asset price over period
        const double r,    // Risk free rate
        const double v,    // Volatility of underlying asset
        const double T);

    double a_2(const double S,    // Spot price
        const double H,    // Min/max of asset price over period
        const double r,    // Risk free rate
        const double v,    // Volatility of underlying asset
        const double T);
    double a_3(const double S,    // Spot price
        const double H,    // Min/max of asset price over period
        const double r,    // Risk free rate
        const double v,    // Volatility of underlying asset
        const double T);

};

#endif // LOOKBACK_H

