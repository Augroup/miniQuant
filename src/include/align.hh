#ifndef ALIGN_H
#define ALIGN_H


#include <cstddef>
#include <string>
#include <vector>
#include <memory>
#include <cstdint>
#include <common.hh>
#include <EM.hh>
#include <ggcat.hh>
enum GraphStrand {
    // forward: same strand as reference sequences
    // reverse: opposite strand as reference sequence
    forward = 1,
    reverse = 2
};
struct GraphParam{
    size_t kmer_size;
    GraphStrand graph_strand;
    bool same_strand_only;
    GraphParam(size_t kmer_size,GraphStrand graph_strand,bool same_strand_only) 
        : kmer_size(kmer_size), graph_strand(graph_strand), same_strand_only(same_strand_only){}
    bool operator==(const GraphParam& other) const {
        return kmer_size == other.kmer_size && graph_strand == other.graph_strand && same_strand_only == other.same_strand_only;
    }
    std::string toString(){
        std::ostringstream outs;
        outs << "k" << kmer_size<< "_" << "strand" << graph_strand << "_"<<"same_strand_only" << same_strand_only <<".graph.fa.lz4";
        return outs.str();
    }
};
template <>
struct std::hash<GraphParam>{
    std::size_t operator()(const GraphParam& graph) const {
        // Use the hash function of std::string and XOR the results
        return ((std::hash<size_t>()(graph.kmer_size)
             ^ (std::hash<size_t>()(graph.graph_strand) << 1)) >> 1)
             ^ (std::hash<bool>()(graph.same_strand_only) << 1);
    }
};
void runAlign(std::vector<Transcript>& transcripts,ProgramOptions &opt);
void create_graph(std::string &graph_file_name,size_t kmer_size,GraphStrand strand,bool same_strand_only,ggcat::GGCATInstance *instance,ProgramOptions &opt);
   

#endif