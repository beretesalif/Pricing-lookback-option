#include "pch.h"
#include "PathSimulator.h"
#include <iostream>
#include <algorithm>
#include <cmath>
#include <tuple>
#include <random>
#include <sstream>


PathSimulator::PathSimulator() : S_(100), r_(0.05), T_(1.0), sigma_(0.2), numberOfDiscretizations_(100), S_vect_(100,10), randomValues_(10, 0) {
    // Initialise S_vect_ avec 100 éléments
}

PathSimulator::PathSimulator(double S, double r, double T, double sigma, int N)
    : S_(S), r_(r), T_(T), sigma_(sigma), numberOfDiscretizations_(N), S_vect_(N, S), randomValues_(N, S) {
    // Initialisation du vecteur de chemin ou autres opérations nécessaires
}

PathSimulator::PathSimulator(const PathSimulator& other)
    : S_(other.S_), r_(other.r_), T_(other.T_), sigma_(other.sigma_), S_vect_(other.S_vect_), numberOfDiscretizations_(other.numberOfDiscretizations_), randomValues_(other.randomValues_) {
    // Copie profonde du vecteur de chemin ou autres opérations nécessaires
}

PathSimulator::~PathSimulator() {
    // Libération des ressources ou autres opérations nécessaires
}

void PathSimulator::simulatePath() {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);

    randomValues_.clear();


    for (size_t i = 1; i < S_vect_.size(); i++) {
        double gauss_bm = distribution(gen);
        randomValues_.push_back(gauss_bm);  // Stocke chaque nombre généré

        double drift = std::exp((r_ - 0.5 * sigma_ * sigma_) * (T_ / numberOfDiscretizations_));
        double vol = std::sqrt(sigma_ * sigma_ * (T_ / numberOfDiscretizations_));

        S_vect_[i] = S_vect_[i - 1] * drift * std::exp(vol * gauss_bm);
    }
}


double PathSimulator::getMin() const {
    if (S_vect_.empty()) {
        return 0.0; // Gérer le cas où le vecteur est vide
    }
    return *std::min_element(S_vect_.begin(), S_vect_.end());
}

double PathSimulator::getMax() const {
    if (S_vect_.empty()) {
        return 0.0; // Gérer le cas où le vecteur est vide
    }
    return *std::max_element(S_vect_.begin(), S_vect_.end());
}

double PathSimulator::getFinalValue() const {
    if (S_vect_.empty()) {
        return 0.0; // Gérer le cas où le vecteur est vide
    }
    return S_vect_.back();
}


void PathSimulator::recalculatePathWithPreviousValues() {
    S_vect_.clear();
    S_vect_.reserve(getNomberOfDiscratizations());

    S_vect_.push_back(S_);  // Ajoutez la valeur initiale

    double epsilon = (T_) / getNomberOfDiscratizations();
    double drift = std::exp(epsilon * (r_ - 0.5 * sigma_ * sigma_));
    double vol = std::sqrt(epsilon * sigma_ * sigma_);

    // Vérifiez que randomValues_ a suffisamment d'éléments
    if (randomValues_.size() < getNomberOfDiscratizations() - 1) {
        // Gérer l'erreur, peut-être lever une exception ou prendre une autre mesure appropriée
        std::cerr << "Erreur : Pas assez de valeurs aléatoires pour la discrétisation." << std::endl;
        return;  // Quittez la fonction en cas d'erreur
    }

    for (int i = 1; i < getNomberOfDiscratizations(); i++) {
        double gauss_bm = randomValues_[i - 1];
        S_vect_.push_back(S_vect_.back() * drift * std::exp(vol * gauss_bm));
    }
}


void PathSimulator::setT(double T) {
    T_ = T;
    recalculatePathWithPreviousValues();
}

void PathSimulator::setr(double r) {
    r_ = r;
    recalculatePathWithPreviousValues();
}

void PathSimulator::setsigma(double sigma) {
    sigma_ = sigma;
    recalculatePathWithPreviousValues();
}

void PathSimulator::setS(double S) {
    S_ = S;
    recalculatePathWithPreviousValues();
}

const double PathSimulator::GetS() const {
    return S_;
}

const double PathSimulator::GetSigma() const {
    return sigma_;
}

const double PathSimulator::GetR() const {
    return r_;
}
const double PathSimulator::GetT() const {
    return T_;
}



const std::vector<double>& PathSimulator::getSvect() const {
    return S_vect_;
}

const int PathSimulator::getNomberOfDiscratizations() const {
    return numberOfDiscretizations_;
}

