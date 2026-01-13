#ifndef COND_PROB_ESTIMATOR_H
#define COND_PROB_ESTIMATOR_H

#include <vector>
#include <unordered_map>
#include <utility>
#include <cstddef>
#include <EM.hh>
#include <lengthDensityEstimator.hh>

class CondProbEstimator {
public:
    void estimateCondProb(LengthDensityEstimator &estimator,std::vector<Transcript>& transcripts, std::vector<Read>& reads, size_t start, size_t end);
    void parallelEstimateCondProb(LengthDensityEstimator &estimator,std::vector<Transcript>& transcripts, std::vector<Read>& reads, ProgramOptions &opt);

};

#endif