#include <lengthDensityEstimator.hh>
#include <EM.hh>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <util.hh>

// LengthDensityEstimator constructor
LengthDensityEstimator::LengthDensityEstimator(){
}

// Calculate the Kaplan-Meier death estimate
void LengthDensityEstimator::calculateKaplanMeier(const std::vector<Read>& data, const std::vector<bool>& censored) {
    // program_log_message(Info, "Estimating long reads length distribution by Kaplan-Meier...");
    std::vector<Read> sortedData = data;
    std::sort(sortedData.begin(), sortedData.end(), [](const Read& a, const Read& b) {
        return a.length < b.length;
    });

    double n = static_cast<double>(sortedData.size());
    double cumNumReads = 0.0;
    double deathProb = 1.0;
    for (size_t i = 0; i < sortedData.size(); ++i) {
        // kmEstimateMap[sortedData[i].length] = i/n;
        // std::cout << sortedData[i].length << std::endl;
        // std::cout << censored[i] << std::endl;
        if (!censored[i]) {
            // Count events
            size_t eventCount = 1;
            while (i + 1 < sortedData.size() && sortedData[i].length == sortedData[i + 1].length && !censored[i + 1]) {
                ++eventCount;
                ++i;
            }
            deathProb *= ((n - eventCount) / n);
            kmEstimateMap[sortedData[i].length] = 1 - deathProb;
            n -= eventCount;
        } else {
            n -= 1.0;
        }
    }
    // program_log_message(Info, "Done.");
}

// Print the Kaplan-Meier death function
void LengthDensityEstimator::printKMDeathFunction() const {
    std::cout << "Read Length\tDeath Probability\n";
    for (const auto& entry : kmEstimateMap) {
        std::cout << entry.first << "\t\t" << std::fixed << std::setprecision(2) << entry.second << "\n";
    }
}

// Get death probability for a specific read length
double LengthDensityEstimator::getDeathProbability(size_t readLength) const {
    auto it = kmEstimateMap.find(readLength);
    if (it != kmEstimateMap.end()) {
        return it->second;
    } else {
        // If the exact length is not found, find the closest higher length
        double lastProb = 1.0;
        for (const auto& entry : kmEstimateMap) {
            if (entry.first > readLength) {
                lastProb = entry.second;
                break;
            }
            lastProb = entry.second;
        }
        return lastProb;
    }
}
