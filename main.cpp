#include <iostream>
#include <vector>
#include <cmath>

class PathSimulator {
public:
    // Constructor
    PathSimulator() {}

    // Destructor
    ~PathSimulator() {}

    // Member function for simulating the path
    void simulatePath() {
        double epsilon = (T_ - t_) / static_cast<double>(S_vect_.size());
        double drift = std::exp(epsilon * (r_ - 0.5 * sigma_ * sigma_));
        double vol = std::sqrt(epsilon * sigma_ * sigma_);

        for (int i = 1; i < S_vect_.size(); i++) {
            double gauss_bm = gaussianBoxMuller();
            S_vect_[i] = S_vect_[i - 1] * drift * std::exp(vol * gauss_bm);
        }
    }

    // Setters for parameters
    void setT(double T) {
        T_ = T;
    }

    void sett(double t) {
        t_ = t;
    }

    void setr(double r) {
        r_ = r;
    }

    void setsigma(double sigma) {
        sigma_ = sigma;
    }

    // Getter for S_vect
    const std::vector<double>& getSvect() const {
        return S_vect_;
    }

private:
    // Member function for generating Gaussian random numbers using Box-Muller method
    double gaussianBoxMuller() {
        double x = 0.0;
        double y = 0.0;
        double euclidsq = 0.0;

        do {
            x = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
            y = 2.0 * rand() / static_cast<double>(RAND_MAX) - 1;
            euclidsq = x * x + y * y;
        } while (euclidsq >= 1.0);

        return x * std::sqrt(-2 * std::log(euclidsq) / euclidsq);
    }

    // Parameters
    double t_ = 0.0;
    double r_ = 0.05;
    double T_ = 1.0;
    double sigma_ = 0.2;

    // Path vector
    std::vector<double> S_vect_;
};

int main() {
    // Example of using the PathSimulator class
    PathSimulator pathSimulator;

    // Initialize parameters and simulate the path
    pathSimulator.simulatePath();

    // Display the initial results
    std::cout << "Initial Path: ";
    for (const auto& value : pathSimulator.getSvect()) {
        std::cout << value << " ";
    }

    // Modify parameters and simulate the path again
    pathSimulator.setT(2.0);
    pathSimulator.setr(0.1);
    pathSimulator.simulatePath();

    // Display the modified results
    std::cout << "\nModified Path: ";
    for (const auto& value : pathSimulator.getSvect()) {
        std::cout << value << " ";
    }

    return 0;
}
