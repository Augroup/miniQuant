#ifndef LENGTH_DENSITY_ESTIMATOR_H
#define LENGTH_DENSITY_ESTIMATOR_H

#include <vector>
#include <map>
#include <utility>
#include <cstddef>
#include <EM.hh>

class LengthDensityEstimator {
public:
    // Constructor
    LengthDensityEstimator();

    // Calculate the Kaplan-Meier estimator
    void calculateKaplanMeier(const std::vector<Read>& data, const std::vector<bool>& censored);

    // Print the Kaplan-Meier death function
    void printKMDeathFunction() const;

    // Get death probability for a specific read length
    double getDeathProbability(size_t readLength) const;

private:
    std::vector<Read> data;
    std::vector<bool> censored;
    std::map<size_t, double> kmEstimateMap; // Store KM estimates
};

#endif