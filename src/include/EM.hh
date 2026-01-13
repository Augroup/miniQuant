#ifndef EM_H // Include guard to prevent multiple inclusions
#define EM_H

#include <iostream>
#include <vector>
#include <unordered_map>
#include <thread>
#include <mutex>
#include <common.hh> // Ensure this header is properly set up for your project
#include <sstream>

// Assuming Transcript is defined like this
struct Transcript {
    // internal ID for transcript
    size_t id;
    // transcript ID shown in annotation
    std::string name;
    // transcript length
    size_t length;
    // effective length for short reads, only dependes on the fragment length and transcript length
    double short_reads_effective_length;
    // relative abundance, with all transcript sum to 1
    double relative_abund;
    // relative abundance multiply by effective length, and normalized for all transcripts.
    double short_reads_effective_count;
    // double true_relative_abund;
    // double true_short_reads_effective_count;
    // relative abundance in TPM
    double abundance;
    Transcript()
        : id(-1), name(""), length(1), short_reads_effective_length(1.0),relative_abund(0.0), short_reads_effective_count(0.0),abundance(0) {}
    Transcript(size_t id, std::string name, size_t length)
        : id(id), name(name), length(length), short_reads_effective_length(1.0), relative_abund(0.0), short_reads_effective_count(0.0),abundance(0) {}
    std::string toString(){
        std::ostringstream outs;
        outs << "Transcript(ID:" << id << ", Name:" << name << ", Eff_length:"<< short_reads_effective_length << ");";
        return outs.str(); 
    }
    std::string toJSONL(){
        std::ostringstream outs;
        outs << "{\"transcript_id\":" << id << ",\"transcript_name\":\"" << name<<"\",\"length\":"<<length << ",\"short_reads_effective_length\":"<<short_reads_effective_length << ",\"relative_abund\":"<<relative_abund <<",\"abundance\":"<<abundance << ",\"short_reads_effective_count\":"<<short_reads_effective_count << "}\n";
        return outs.str();
    }
};

// Assuming Hit is defined like this
struct Hit {
    size_t transcript_id;
    size_t mapped_kmer;
    size_t mapped_kmer_R1;
    size_t mapped_kmer_R2;
    double cond_prob;
    Hit()
        : transcript_id(-1), cond_prob(0) {}
    Hit(size_t transcript_id, double cond_prob)
        : transcript_id(transcript_id), cond_prob(cond_prob) {}
    std::string toString(){
        std::ostringstream outs;
        outs << "Hit(Transcript ID:" << transcript_id << ", Cond prob:" << cond_prob << ");";
        return outs.str();
        
    }
    std::string toJSONL(ReadType readType){
        std::ostringstream outs;
        if (readType == LongReadType){
            outs << "\"" << transcript_id << "\"" << ":{\"cond prob\":" << cond_prob <<",\"mapped_kmer\":"<<mapped_kmer<<"}";
        } else if (readType == ShortReadType){
            outs << "\"" << transcript_id << "\"" << ":{\"cond prob\":" << cond_prob <<",\"mapped_kmer_R1\":"<<mapped_kmer_R1<<",\"mapped_kmer_R2\":"<<mapped_kmer_R2<<"}";
        }
        return outs.str();
    }
    void set_mapped_kmer(size_t _mapped_kmer){
        mapped_kmer = _mapped_kmer;
    }
    void set_mapped_kmer(size_t _mapped_kmer_R1,size_t _mapped_kmer_R2){
        mapped_kmer_R1 = _mapped_kmer_R1;
        mapped_kmer_R2 = _mapped_kmer_R2;
    }
};


// Assuming Read is defined like this
struct Read {
    size_t id;
    size_t length;
    size_t all_mapped_kmer;
    size_t all_kmer;
    // debug
    size_t all_mapped_kmer_R1;
    size_t all_kmer_R1;
    size_t all_mapped_kmer_R2;
    size_t all_kmer_R2;
    
    std::vector<Hit> hits;
    ReadType readType;
    Read()
        : length(-1), hits({}),readType(LongReadType) {}
    Read(size_t length, std::vector<Hit> hits,ReadType readType)
        : length(length), hits(hits), readType(readType){}
    std::string toString(){
        std::ostringstream outs;
        outs << "Read(Length " << length;
        return outs.str();
    }
    std::string toJSONL(){
        std::ostringstream outs;
        if (readType == LongReadType){
            outs << "{\"read_id\":"<< id <<",\"length\":" << length << ",\"all_mapped_kmer\":" <<all_mapped_kmer<<",\"all_kmer\":"<<all_kmer << ",\"hits\":{";
            bool first_hit = true;
            for (auto & hit : hits){
                if (first_hit){
                    outs << hit.toJSONL(readType);
                    first_hit = false;
                } else{
                    outs << "," << hit.toJSONL(readType);
                } 
            }
            outs << "}}\n";
        } else if (readType == ShortReadType){
            outs << "{\"read_id\":"<< id <<",\"length\":" << length << ",\"all_mapped_kmer_R1\":" <<all_mapped_kmer_R1<<",\"all_kmer_R1\":"<<all_kmer_R1 <<",\"all_mapped_kmer_R2\":" <<all_mapped_kmer_R2<<",\"all_kmer_R2\":"<<all_kmer_R2 << ",\"hits\":{";
            bool first_hit = true;
            for (auto & hit : hits){
                if (first_hit){
                    outs << hit.toJSONL(readType);
                    first_hit = false;
                } else{
                    outs << "," << hit.toJSONL(readType);
                } 
            }
            outs << "}}\n";
        } 
        return outs.str();
    }
    void set_id(size_t _id){
        id = _id;
    }
    void set_mapped_kmer(size_t _all_mapped_kmer,size_t _all_kmer){
        all_mapped_kmer = _all_mapped_kmer;
        all_kmer = _all_kmer;
    }
    void set_mapped_kmer(size_t _all_mapped_kmer_R1,size_t _all_kmer_R1,size_t _all_mapped_kmer_R2,size_t _all_kmer_R2){
        all_mapped_kmer_R1 = _all_mapped_kmer_R1;
        all_kmer_R1 = _all_kmer_R1;
        all_mapped_kmer_R2 = _all_mapped_kmer_R2;
        all_kmer_R2 = _all_kmer_R2;
    }
    // void set_weight(double _weight){
    //     weight = _weight;
    // }
};

void runEM(std::vector<Read>& reads,std::vector<Transcript>& transcripts, ProgramOptions& opt);
void initializeAbundance(std::vector<Transcript>& transcripts);
void EMmanager(std::vector<Transcript>& transcripts, std::vector<Read>& reads,  std::unordered_map<size_t, double>& readFraction,ProgramOptions& opt);
void expectationStep(const std::vector<Transcript>& transcripts, std::vector<Read>& reads,  std::unordered_map<size_t, double>& readFraction, size_t start, size_t end);
void parallelExpectationStep(std::vector<Transcript>& transcripts, std::vector<Read>& reads,  std::unordered_map<size_t, double>& readFraction, int numThreads);
void calTPM(std::vector<Transcript>& transcripts);

#endif // EM_H