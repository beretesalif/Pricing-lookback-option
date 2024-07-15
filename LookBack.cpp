#include "pch.h"
#include "LookBack.h"
#include <cmath>
#include <string>

// Implémentation du constructeur paramétrique
LookBack::LookBack(const char* optionType, double epsilon)
    : optionType_(optionType), epsilon_(epsilon) {}

// Implémentation des méthodes de calcul
double LookBack::Payoff(const PathSimulator& pathSimulator) const {
    double finalValue = pathSimulator.getFinalValue();
    double minPath = pathSimulator.getMin();
    double maxPath = pathSimulator.getMax();

    if (std::string(optionType_) == "Call") {
        return MAX(finalValue - minPath, 0.0);
    }
    else if (std::string(optionType_) == "Put") {
        return MAX(maxPath - finalValue, 0.0);
    }
    else {
        // Gestion d'une option non reconnue
        return 0.0;
    }
}

double LookBack::Delta(PathSimulator pathSimulator) const {
    double originalS = pathSimulator.GetS();

    // Calcul du payoff avec S + epsilon
    pathSimulator.setS(originalS * (1 + epsilon_));
    double payoffPlus = Payoff(pathSimulator);

     // Réinitialisez la valeur de S à sa valeur initiale
    pathSimulator.setS(originalS);

    // Calcul du payoff avec S - epsilon
    pathSimulator.setS(originalS * (1 - epsilon_));
    double payoffMinus = Payoff(pathSimulator);

    return (payoffPlus - payoffMinus) / (2 * epsilon_ * originalS);
}

double LookBack::Gamma(PathSimulator pathSimulator) const {
    double originalS = pathSimulator.GetS();
    double originalPayoff = Payoff(pathSimulator);  // Calcul du payoff avec le vecteur initial

    // Calcul du payoff avec S + epsilon
    pathSimulator.setS(originalS * (1 + epsilon_));
    double payoffPlus = Payoff(pathSimulator);

    // Réinitialisez la valeur de S à sa valeur initiale
    pathSimulator.setS(originalS);

    // Calcul du payoff avec S - epsilon
    pathSimulator.setS(originalS * (1 - epsilon_));
    double payoffMinus = Payoff(pathSimulator);

    // Utilisez le payoff initial stocké pour éviter le recalcul redondant
    return (payoffPlus - 2 * originalPayoff + payoffMinus) / ((epsilon_ * originalS) * (epsilon_ * originalS));
}

double LookBack::Vega(PathSimulator pathSimulator) const {
    double originalSigma = pathSimulator.GetSigma();

    // Calcul du payoff avec Sigma + epsilon
    pathSimulator.setsigma(originalSigma * (1 + epsilon_));
    double payoffPlus = Payoff(pathSimulator);

    // Réinitialisez la valeur de Sigma à sa valeur initiale
    pathSimulator.setsigma(originalSigma);

    // Calcul du payoff avec Sigma - epsilon
    pathSimulator.setsigma(originalSigma * (1 - epsilon_));
    double payoffMinus = Payoff(pathSimulator);

    // Utilisez le payoff initial stocké pour éviter le recalcul redondant
    return (payoffPlus - payoffMinus) / (2 * epsilon_ * originalSigma*100);
}

double LookBack::Theta(PathSimulator pathSimulator) const {
    double originalT = pathSimulator.GetT();

    // Calcul du payoff avec T + epsilon
    pathSimulator.setT(originalT * (1 + epsilon_));
    double payoffPlus = Payoff(pathSimulator);

    // Réinitialisez la valeur de T à sa valeur initiale
    pathSimulator.setT(originalT);

    // Calcul du payoff avec T - epsilon
    pathSimulator.setT(originalT * (1 - epsilon_));
    double payoffMinus = Payoff(pathSimulator);

    // Utilisez le payoff initial stocké pour éviter le recalcul redondant
    return (exp(-pathSimulator.GetR() * pathSimulator.GetT() * (1 + epsilon_)) * payoffPlus - exp(-pathSimulator.GetR() * pathSimulator.GetT() * (1 - epsilon_)) * payoffMinus) / (2 * epsilon_ * originalT*365);
}

double LookBack::Rho(PathSimulator pathSimulator) const {
    double originalR = pathSimulator.GetR();

    // Calcul du payoff avec T + epsilon
    pathSimulator.setr(originalR * (1 + epsilon_));
    double payoffPlus = Payoff(pathSimulator);

    // Réinitialisez la valeur de R à sa valeur initiale
    pathSimulator.setr(originalR);

    // Calcul du payoff avec T - epsilon
    pathSimulator.setr(originalR * (1 - epsilon_));
    double payoffMinus = Payoff(pathSimulator);

    // Utilisez le payoff initial stocké pour éviter le recalcul redondant
    return (exp(-pathSimulator.GetR() * pathSimulator.GetT() * (1 + epsilon_)) * payoffPlus - exp(-pathSimulator.GetR() * pathSimulator.GetT() * (1 - epsilon_)) * payoffMinus) / (2 * epsilon_ * originalR*100);
}


double LookBack::N(double x) const {
    double k = 1.0 / (1.0 + 0.2316419 * x);
    double k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));

    if (x >= 0.0) {
        return (1.0 - (1.0 / (pow(2 * PI, 0.5))) * exp(-0.5 * x * x) * k_sum);
    }
    else {
        return 1.0 - N(-x);
    }
}

double LookBack::a_1(const double S,    // Spot price
    const double H,    // Min/max of asset price over period
    const double r,    // Risk free rate
    const double v,    // Volatility of underlying asset
    const double T) {  // Time to expiry
    double num = log(S / H) + (r + 0.5 * v * v) * T;
    double denom = v * sqrt(T);
    return num / denom;
}

double LookBack::a_2(const double S,    // Spot price
    const double H,    // Min/max of asset price over period
    const double r,    // Risk free rate
    const double v,    // Volatility of underlying asset
    const double T) {
    return a_1(S, H, r, v, T) - v * sqrt(T);
}

double LookBack::a_3(const double S,    // Spot price
    const double H,    // Min/max of asset price over period
    const double r,    // Risk free rate
    const double v,    // Volatility of underlying asset
    const double T) {
    return a_1(S, H, r, v, T) - (2.0 * r * sqrt(T) / v);
}


double LookBack::PayoffApproxim(PathSimulator pathSimulator) {
    double r = pathSimulator.GetR();
    double T = pathSimulator.GetT();
    double sigma = pathSimulator.GetSigma();
    double S = pathSimulator.GetS();


    
    double a1 = a_1(S, S, r, sigma, T);
    double a2 = a_2(S, S, r, sigma, T);
    double a3 = a_3(S, S, r, sigma, T);


    double term1 = S * N(a1);
    double term2 = S * exp(-r * T) * N(a2);
    double mult = S * sigma * sigma / (2.0 * r);
    double term3 = N(-a1) - exp(-r * T) * pow((S / S), ((2 * r) / (sigma * sigma))) * N(-a3);
    double Callfloat = term1 - term2 - mult * term3;

    double b1 = a_1(S, S, r, sigma, T);
    double b2 = a_2(S, S, r, sigma, T);
    double b3 = a_3(S, S, r, sigma, T);

    double term11 = -S * N(-b1);
    double term21 = S * exp(-r * T) * N(-b2);
    double mult1 = S * sigma * sigma / (2.0 * r);
    double term31 = N(b1) - exp(-r * T) * pow((S / S), ((2 * r) / (sigma * sigma))) * N(b3);

    double Putfloat = term11 + term21 + mult1 * term31;

    // Choose to return Callfloat or Putfloat based on the option type
    if (std::string(optionType_) == "Call") {
        return Callfloat;
    }
    else if (std::string(optionType_) == "Put") {
        return Putfloat;
    }
    else {
        // Gestion d'une option non reconnue
        return 0.0;
    }
}
