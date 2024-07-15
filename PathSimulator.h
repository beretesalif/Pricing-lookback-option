#ifndef PATH_SIMULATOR_H
#define PATH_SIMULATOR_H

#include <vector>

class PathSimulator {
public:
    // Constructeur par d�faut
    PathSimulator();

    // Constructeur param�trique
    PathSimulator(double S, double r, double T, double sigma, int discretization_steps);

    // Constructeur par copie
    PathSimulator(const PathSimulator& other);
    double getMin() const;
    double getMax() const;
    double getFinalValue() const;

    // Destructeur
    ~PathSimulator();

    // Fonction pour simuler le chemin
    void simulatePath();
    void recalculatePathWithPreviousValues();

    // Setters pour les param�tres
    void setT(double T);
    void setr(double r);
    void setsigma(double sigma);
    void setS(double S);

    // Getter pour S_vect
    const std::vector<double>& getSvect() const;
    const int getNomberOfDiscratizations()const;
    const double GetS()const;

    const double GetSigma() const;
    const double GetR() const;
    const double GetT() const;

private:
    // Fonction pour g�n�rer des nombres al�atoires gaussiens en utilisant la m�thode de Box-Muller

    // Param�tres
    double r_;
    double T_;
    double sigma_;
    double S_;

    // Nombre de discr�tisations
    size_t numberOfDiscretizations_;

    // Vecteur de chemin
    std::vector<double> S_vect_;
    std::vector<double> randomValues_;

};

#endif // PATH_SIMULATOR_H
