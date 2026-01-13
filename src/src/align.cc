
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <string>
#include <mutex>
#include <cstring>
#include <align.hh>
#include <EM.hh>
#include <util.hh>
#include <sequenceReader.hh>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/utility/views/chunk.hpp>
using namespace ggcat;
// whether GGCAT output the information
bool debug_GGCAT = false;
// parse reference sequences

class MemorySequencesReader : public ggcat::StreamReader
{
public:
    MemorySequencesReader() {}
    static uint64_t estimated_base_count(void *block)
    {
        std::istringstream iss((std::string)(char *)block);
        std::string color_id;
        std::string seq;
        std::getline(iss,color_id);
        std::getline(iss,seq);
        return strlen(seq.data());
    }

    // This is a toy example where the block corresponds to a single sequence.
    // In practice a block should be a chunk of sequences, e.g. a pointer to a file descriptor
    // that can be read to retrieve the sequences.
    void read_block(
        void *block,
        bool copy_ident_data,
        size_t partial_read_copyback,
        void (*callback)(DnaSequence sequence, SequenceInfo info)) override
    {
        std::istringstream iss((std::string)(char *)block);
        std::string color_id;
        std::string seq;
        std::getline(iss,color_id);
        std::getline(iss,seq);
        callback(DnaSequence{
                     // To build a graph the ident_data can be empty
                     .ident_data = Slice<char>(nullptr, 0),
                     .seq = Slice<char>((char *)seq.data(), strlen(seq.data())),
                 },
                 SequenceInfo{
                     // This is the index of the color of the current sequence
                     .color = static_cast<unsigned int>(std::stoi(color_id.data())),
                 });
    }
};
void log_message(MessageLevel level, const char *message)
{
    switch (level)
    {
    case MessageLevel_Info:
    {
        if (debug_GGCAT){
            printf("[INFO][GGCAT] %s\n", message);
        }
       
    }
    break;
    case MessageLevel_Warning:
    {
        printf("[WARNING][GGCAT] %s\n", message);
    }
    break;
    case MessageLevel_Error:
    {
        printf("[ERROR][GGCAT] %s\n", message);
    }
    break;
    case MessageLevel_UnrecoverableError:
    {
        printf("[FATAL ERROR][GGCAT] %s\n", message);
        abort();
    }
    }
}
std::unordered_map<GraphParam,std::string> graph_file_path_map;
std::vector<std::string> reference_sequences;
std::vector<std::string> sequenceIDs;
char ** sequences_char = nullptr;
char ** reverse_complement_sequences_char = nullptr;
void create_graph(std::string &graph_file_name,size_t kmer_size,GraphStrand strand,bool same_strand_only,GGCATInstance *instance,ProgramOptions &opt){
    GraphParam new_graph_param(kmer_size,strand,same_strand_only);
    if (graph_file_path_map.count(new_graph_param) == 0){
        // initialize the sequence pointer if not
        if (strand == forward & !sequences_char){
            getSequencesCharPtr(reference_sequences,&sequences_char);
        }
        if (strand == reverse & !reverse_complement_sequences_char){
            std::vector<std::string> reverse_complement_reference_sequences;
            reverse_complement_reference_sequences.reserve(reference_sequences.size());
            reverse_complementary_for_all_sequences(reference_sequences,reverse_complement_reference_sequences);
            getSequencesCharPtr(reverse_complement_reference_sequences,&reverse_complement_sequences_char);
        };
        char** sequences = nullptr;
        if (strand == reverse){
            sequences = reverse_complement_sequences_char;
        } else if (strand == forward){
            sequences = sequences_char;
        }
        std::string graph_file_path = opt.tempFolder+"/"+new_graph_param.toString();
        // Graph not exist, create it
        size_t min_multiplicity = 1;
        if (kmer_size > 64) {
            min_multiplicity = 2;
        }
        instance->build_graph_from_streams<MemorySequencesReader>(
            Slice<void *>((void **)sequences,sequenceIDs.size()),
            graph_file_path,
            kmer_size,
            opt.threads,
            same_strand_only,
            min_multiplicity,
            opt.extra_elaboration_step,
            true,
            Slice<std::string>(sequenceIDs.data(),sequenceIDs.size()),opt.minimizer_size);
        graph_file_path_map[new_graph_param] = graph_file_path;
    };
    graph_file_name = graph_file_path_map[new_graph_param];

}
void runAlign(std::vector<Transcript>& transcripts, ProgramOptions& opt)
{
    debug_GGCAT = opt.debug;
    size_t ggcat_threads = opt.threads;
    GGCATConfig config;

    config.use_temp_dir = true;
    config.temp_dir = opt.tempFolder,
    config.memory = opt.mem * 0.5,
    config.prefer_memory = true,
    config.total_threads_count = ggcat_threads,
    config.intermediate_compression_level = 0,
    config.use_stats_file = false;
    config.stats_file = opt.outputFolder+"/stats.log";
    config.messages_callback = log_message;
    
    GGCATInstance *instance = GGCATInstance::create(config);

    std::string short_reads_R1_graph_file = "";
    std::string short_reads_R2_graph_file = "";
    std::string long_reads_graph_file = "";

    program_log_message(Info,"Building index on reference sequences...");
    readFASTA(opt.transcriptsFile,transcripts, reference_sequences);
    sequenceIDs = std::vector<std::string>(transcripts.size());
    std::transform(transcripts.begin(), transcripts.end(), sequenceIDs.begin(), [](Transcript transcript) {
        return transcript.name;
    });
    
    if (opt.use_long_reads){
        if (opt.long_reads_library_prep == dRNA_ONT){
            create_graph(long_reads_graph_file,opt.long_reads_kmer_size,forward,true,instance,opt);
        } else{
            // cDNA-ONT and cDNA-PacBio are unstranded, simply use same graph
            create_graph(long_reads_graph_file,opt.long_reads_kmer_size,forward,false,instance,opt);
        }
    }
    if (opt.use_short_reads){
        if (opt.short_reads_strandness == unstranded){
        // if unstranded, simply use the same graph for both ends
            create_graph(short_reads_R1_graph_file,opt.short_reads_kmer_size,forward,false,instance,opt);
            short_reads_R2_graph_file = short_reads_R1_graph_file;
        } else if (opt.short_reads_strandness == fr_stranded){
            // if first pair is forward
            create_graph(short_reads_R1_graph_file,opt.short_reads_kmer_size,forward,true,instance,opt);
            create_graph(short_reads_R2_graph_file,opt.short_reads_kmer_size,reverse,true,instance,opt);
        } else if (opt.short_reads_strandness == rf_stranded){
            // if first pair is reverse
            create_graph(short_reads_R1_graph_file,opt.short_reads_kmer_size,reverse,true,instance,opt);
            create_graph(short_reads_R2_graph_file,opt.short_reads_kmer_size,forward,true,instance,opt);
        };

    }

    program_log_message(Info,"Done.");
    if (opt.use_short_reads){
        program_log_message(Info,"Aligning short reads...");
        // TODO: add error handling if aligning short reads failed
        emptyAndCreateFolder(opt.tempFolder+"/short_reads_alignments_R1/");
        emptyAndCreateFolder(opt.tempFolder+"/short_reads_alignments_R2/");
        instance->query_graph(
            short_reads_R1_graph_file,
            opt.shortreadsFilePairOne,
            opt.tempFolder+"/short_reads_alignments_R1/bucket_.jsonl.lz4",
            opt.short_reads_kmer_size,
            ggcat_threads,
            opt.forward_only,
            true,
            ColoredQueryOutputFormat_JsonLinesWithNumbers,opt.minimizer_size);
        instance->query_graph(
            short_reads_R2_graph_file,
            opt.shortreadsFilePairTwo,
            opt.tempFolder+"/short_reads_alignments_R2/bucket_.jsonl.lz4",
            opt.short_reads_kmer_size,
            ggcat_threads,
            opt.forward_only,
            true,
            ColoredQueryOutputFormat_JsonLinesWithNumbers,opt.minimizer_size);
        program_log_message(Info,"Done.");
    }

    if (opt.use_long_reads){
        program_log_message(Info,"Aligning long reads...");
        emptyAndCreateFolder(opt.tempFolder+"/long_reads_alignments/");
        size_t single_chunk_num_long_reads = 10000000;
         // TODO: add error handling if aligning long reads failed
        seqan3::sequence_file_input fin{opt.longreadsFile,seqan3::fields<seqan3::field::id, seqan3::field::seq>{}}; // auto-gz supported
        opt.num_long_reads_chunks = 0;
        for (auto && records : fin | seqan3::views::chunk(single_chunk_num_long_reads)){
            seqan3::sequence_file_output fout{opt.tempFolder+"/temp_long_reads.fasta",
            seqan3::fields<seqan3::field::id, seqan3::field::seq>{}};
            fout = records;
            fout.get_stream().flush();
            instance->query_graph(
                long_reads_graph_file,
                opt.tempFolder+"/temp_long_reads.fasta",
                opt.tempFolder+"/long_reads_alignments/bucket_" + std::to_string(opt.num_long_reads_chunks) +"_.jsonl.lz4",
                opt.long_reads_kmer_size,
                ggcat_threads,
                opt.forward_only,
                true,
                ColoredQueryOutputFormat_JsonLinesWithNumbers,opt.minimizer_size);
            opt.num_long_reads_chunks += 1;
        }
        
        program_log_message(Info,"Done.");
    }

    program_log_message(Info,"Aligning reads done.");
}

