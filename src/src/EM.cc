#include <iostream>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <vector>
#include <EM.hh>

#include <alignmentReader.hh>
#include <condProbEstimator.hh>
#include <util.hh>
#include <fileWriter.hh>
std::mutex mtx;
double min_relative_abund = 1e-10;
int maxIterations = 1000;
double tolerance = 1e-6;

void initializeAbundance(std::vector<Transcript>& transcripts) {
    // initialize abundance and effective count by uniform values
    double initialAbundance = 1.0 / transcripts.size();
    double sumEffectiveCount = 0.0;
    for (auto& transcript : transcripts) {
        transcript.relative_abund = initialAbundance;
        transcript.short_reads_effective_count = transcript.relative_abund * transcript.short_reads_effective_length;
        sumEffectiveCount += transcript.short_reads_effective_count;
    }
    if (sumEffectiveCount > 0.0) {
        for (auto& transcript : transcripts) {
            transcript.short_reads_effective_count /= sumEffectiveCount;
        }
    }
}

void calTPM(std::vector<Transcript>& transcripts) {
    double sumFractions = 0.0;
    double uniformTPM = 1e6/transcripts.size();
    for (auto& transcript : transcripts) {
        if (transcript.relative_abund >= 0){
            sumFractions += transcript.relative_abund;
        } else {
            // std::cout << transcript.toJSONL() << std::endl;
        }
        transcript.abundance = 0;
    }
    if (sumFractions > 0) {
        for (auto& transcript : transcripts) {
            transcript.relative_abund /= sumFractions;
            transcript.abundance = transcript.relative_abund * 1e6;
        }
     }
}

void EMmanager(std::vector<Transcript>& transcripts, std::vector<Read>& reads, std::unordered_map<size_t, double>& readFraction,ProgramOptions& opt) {
    std::unordered_map<size_t, double> previousAbundance;
    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        double previousRelativeAbundEffLengthProduct= 0;
        for (auto& transcript : transcripts) {
            readFraction[transcript.id] = 0.0;
            previousRelativeAbundEffLengthProduct += transcript.short_reads_effective_length * transcript.relative_abund;
        }
        parallelExpectationStep(transcripts, reads, readFraction, opt.threads);
        // M step
        double sumRelativeAbund = 0.0;
        double sumEffectiveCount = 0.0;
        for (auto& transcript : transcripts) {
            previousAbundance[transcript.id] = transcript.relative_abund;
            if (previousAbundance[transcript.id] > 0){
                double denominator = opt.long_reads_weight_sum + (opt.short_reads_weight_sum * transcript.short_reads_effective_length/previousRelativeAbundEffLengthProduct);
                if (denominator <= 0)
                    denominator = 1;
                transcript.relative_abund = readFraction[transcript.id]/denominator;
                if (transcript.relative_abund < 0){
                    assert((transcript.relative_abund > 0)||(transcript.relative_abund==0));
                }
                sumRelativeAbund += transcript.relative_abund;   
                sumEffectiveCount += transcript.relative_abund * transcript.short_reads_effective_length;
            }
        }
        double diff = 0.0;
        if (sumRelativeAbund > 0.0) {
            for (auto& transcript : transcripts) {
                transcript.relative_abund /= sumRelativeAbund;
                transcript.short_reads_effective_count = transcript.relative_abund * transcript.short_reads_effective_length / sumEffectiveCount;   
                double transcript_diff = std::abs(transcript.relative_abund - previousAbundance[transcript.id]);
                // if the transcript is very small and has a very low diff, set to zero
                if (transcript_diff < tolerance && transcript.relative_abund < min_relative_abund){
                    transcript.relative_abund = 0;
                    transcript.short_reads_effective_count = 0;
                };
                diff = std::max(diff, transcript_diff);
            }
            if (diff < tolerance) {
                // if max diff < tol, converge.
                program_log_message(Info, ("Converged after "+std::to_string(iteration + 1)+ " iterations.").c_str());
                break;
            }
        }  else {
            // if all transcripts have zero relative abundance, converge.
            program_log_message(Info, ("Converged after "+std::to_string(iteration + 1)+ " iterations.").c_str());
            break;
        }
        if (iteration % 100 == 99) {
            program_log_message(Info, ("Iteration "+std::to_string(iteration)+ ": max relative abund diff: " + std::to_string(diff)).c_str());
            if (opt.debug){
                FileWriter debug_writer;
                debug_writer.writeTranscript(opt.outputFolder+"/"+std::to_string(iteration)+"_transcripts.jsonl",transcripts);
            }
        } 
    }
}

void expectationStep(const std::vector<Transcript>& transcripts, std::vector<Read>& reads, std::unordered_map<size_t, double>& readFraction, size_t start, size_t end) {
    std::unordered_map<size_t, double> localReadFraction;

    for (size_t i = start; i < end; ++i) {
        const auto& read = reads[i];
        double totalCompatibility = 0.0;
        for (const auto& hit : read.hits) {
            const auto& transcript = transcripts[hit.transcript_id];
            if (read.readType == LongReadType){
                totalCompatibility += hit.cond_prob * transcript.relative_abund;
            } else {
                 totalCompatibility += hit.cond_prob * transcript.short_reads_effective_count;
            }
            
        }

        if (totalCompatibility > 0.0) {
            for (const auto& hit : read.hits) {
                const auto& transcript = transcripts[hit.transcript_id];
                if (localReadFraction.count(hit.transcript_id) == 0) {
                    localReadFraction[hit.transcript_id] = 0.0;
                }
                if (read.readType == LongReadType){
                    localReadFraction[hit.transcript_id] += (hit.cond_prob * transcript.relative_abund) / totalCompatibility;
                } else {
                    localReadFraction[hit.transcript_id] += (hit.cond_prob * transcript.short_reads_effective_count) / totalCompatibility;
                }           
            }
        }
    }

    std::lock_guard<std::mutex> lock(mtx);
    for (const auto& entry : localReadFraction) {
        readFraction[entry.first] += entry.second;
    }
}

void parallelExpectationStep(std::vector<Transcript>& transcripts, std::vector<Read>& reads, std::unordered_map<size_t, double>& readFraction, int numThreads) {
    std::vector<std::thread> threads;
    size_t num_reads = reads.size();
    size_t readsPerThread = num_reads / numThreads;
    size_t additionalReads = num_reads % numThreads;

    size_t start = 0;
    for (int i = 0; i < numThreads; ++i) {
        size_t end = start + readsPerThread + (i < additionalReads ? 1 : 0);
        if (start < num_reads) {
            threads.emplace_back(expectationStep, std::cref(transcripts), std::ref(reads), std::ref(readFraction), start, end);
        }
        start = end;
    }

    for (auto& t : threads) {
        t.join();
    }
}

void runEM(std::vector<Read> &reads,std::vector<Transcript>& transcripts, ProgramOptions& opt) {
    std::unordered_map<size_t, double> readFraction; // Declare readFraction locally
    FileWriter writer;
    LengthDensityEstimator length_estimator;
    AlignmentReader alignment_reader;
    if (opt.use_long_reads){
        std::vector<Read> long_reads;
        alignment_reader.parallelProcessAlignments(transcripts, long_reads, opt,LongReadType);
        std::vector<bool> censored(long_reads.size(),false);
        length_estimator.calculateKaplanMeier(long_reads, censored);
        reads.insert(end(reads), begin(long_reads), end(long_reads));
    }
    if (opt.use_short_reads){
        std::vector<Read> short_reads;
        alignment_reader.parallelProcessAlignments(transcripts, short_reads, opt,ShortReadType);
        reads.insert(end(reads), begin(short_reads), end(short_reads));

    }
    program_log_message(Info, "Estimating cond prob...");
    CondProbEstimator cond_prob_estimator;
    cond_prob_estimator.parallelEstimateCondProb(length_estimator,transcripts, std::ref(reads), opt);

    program_log_message(Info, "Done.");
    // writer.readWeight(opt.inputWeightFile,reads,opt);
    program_log_message(Info, "Quantifying...");
    initializeAbundance(transcripts);
    if (opt.debug){
        writer.writeTranscript(opt.outputFolder+"/initial_transcripts.jsonl",transcripts);
    }
    EMmanager(transcripts, reads, readFraction, opt);

    // Output results in TPM
    calTPM(transcripts);
    program_log_message(Info, "Done.");
    program_log_message(Info, "Writing results to file...");
    writer.writeAbundance(opt.outputFolder+"/abundance.tsv",transcripts,opt);
    if (opt.debug){
        writer.writeRead(opt.outputFolder+"/reads.jsonl",reads);
        writer.writeTranscript(opt.outputFolder+"/final_transcripts.jsonl",transcripts);
    }
    program_log_message(Info, "Done.");
}