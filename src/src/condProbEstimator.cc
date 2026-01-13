
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <thread>
#include <cmath>
#include <condProbEstimator.hh>
#include <lengthDensityEstimator.hh>
#include <EM.hh>
#include <common.hh>
#include <util.hh>

size_t calLongReadEffectiveLength(size_t transcript_length,size_t read_length){
    size_t eff_length = transcript_length - read_length + 1;
    if (eff_length <= 0) {
        eff_length = 1;
    }
    return eff_length;
}
// void calAllShortReadEffectiveLength(std::vector<Transcript>& transcripts, double fraglen, double fragsd) {

//     for (size_t j = 0; j < transcripts.size(); ++j) {
//         auto & transcript = transcripts[j];
//         transcript.short_reads_effective_length = std::max(transcript.length -  fraglen + 1.0, 1.0);
//     }
// }
void calAllShortReadEffectiveLength(std::vector<Transcript>& transcripts, double fraglen, double fragsd) {
    std::vector<size_t> all_transcript_length(transcripts.size());
    std::transform(transcripts.begin(), transcripts.end(), all_transcript_length.begin(), [](Transcript transcript) {
        return transcript.length;
    });
    std::vector<double> precomputed_numerators;  // Precompute these based on longest transcript
    std::vector<double> precomputed_denominators;
    const double sqrt_2pi = std::sqrt(2 * M_PI);
    const double normalization_factor = 1 / (fragsd * sqrt_2pi);

    // Find the maximum transcript length to precompute terms
    size_t max_transcript_length = *std::max_element(all_transcript_length.begin(), all_transcript_length.end());
    
    // Precompute for all possible transcript lengths encountered
    precomputed_numerators.reserve(max_transcript_length);
    precomputed_denominators.reserve(max_transcript_length);
    double cumsum_numerators = 0;
    double cumsum_denominators = 0;

    for (size_t length = 1; length <= max_transcript_length; ++length) {
        double exp_term = std::exp(-0.5 * std::pow((length - fraglen) / fragsd, 2));
        double weighted_exp = exp_term * normalization_factor;
        cumsum_numerators += length * weighted_exp;
        cumsum_denominators += weighted_exp;
        precomputed_numerators.push_back(cumsum_numerators);
        precomputed_denominators.push_back(cumsum_denominators);
    }

    for (size_t j = 0; j < transcripts.size(); ++j) {
        auto & transcript = transcripts[j];
        if (transcript.length < fraglen - 5 * fragsd) {
            transcript.short_reads_effective_length = 1.0;
            continue;
        }

        double sum_numerator = precomputed_numerators[transcript.length-1];
        double sum_denominator =  precomputed_denominators[transcript.length-1];

        transcript.short_reads_effective_length = std::max(transcript.length - sum_numerator / sum_denominator + 1.0, 1.0);
    }
}
void CondProbEstimator::estimateCondProb(LengthDensityEstimator &estimator,std::vector<Transcript>& transcripts, std::vector<Read>& reads, size_t start, size_t end){
    for (size_t i = start; i < end; ++i) {
        auto &read = reads[i];
        size_t read_length = read.length;
        for (auto & hit : read.hits){
            size_t transcript_length = transcripts[hit.transcript_id].length;
            size_t transcript_effective_length = 1;
            double read_cond_prob = 0;
            if (read.readType == LongReadType){
                // If read is long reads, cal LR cond prob
                transcript_effective_length = calLongReadEffectiveLength(transcript_length,read_length);
                
                double length_cond_prob = estimator.getDeathProbability(read.length+1) - estimator.getDeathProbability(read.length);
                double cum_length_prob = estimator.getDeathProbability(transcript_length+1) - estimator.getDeathProbability(1);
                if (cum_length_prob > 0){
                    length_cond_prob /= cum_length_prob;
                } 
                double start_site_cond_prob= 1.0/(double)transcript_effective_length;
                read_cond_prob = length_cond_prob * start_site_cond_prob;
            } else {
                // If read is short reads, cal SR cond prob
                transcript_effective_length = transcripts[hit.transcript_id].short_reads_effective_length;
                read_cond_prob = 1.0/(double)transcript_effective_length;
            }
            // Get cond prob and read weight
            hit.cond_prob = read_cond_prob;
        }
    }
}

void CondProbEstimator::parallelEstimateCondProb(LengthDensityEstimator &estimator,std::vector<Transcript>& transcripts, std::vector<Read>& reads,ProgramOptions &opt) {
   if (opt.use_short_reads){
        // If we need to use short reads 
        calAllShortReadEffectiveLength(std::ref(transcripts),opt.fraglen,opt.fragsd);
    }
    size_t num_reads = reads.size();
    std::vector<std::thread> threads;
    size_t readsPerThread = num_reads / opt.threads;
    size_t additionalReads = num_reads % opt.threads;
    // std::cout << "Cal cond prob now:" << std::endl;
    // estimateCondProb(estimator,transcripts,reads,hits,0,num_reads);
    size_t start = 0;
    for (int i = 0; i < opt.threads; ++i) {
        size_t end = start + readsPerThread + (i < additionalReads ? 1 : 0);
        if (start < num_reads) {
                threads.emplace_back(
                &CondProbEstimator::estimateCondProb, // Member function pointer
                this,                                  // Instance of the class
                std::ref(estimator),                   // Pass by reference using std::ref
                std::ref(transcripts),
                std::ref(reads),
                start,
                end
            );
        }
        start = end;
    }

    for (auto &t : threads) {
        t.join();
    }
}
