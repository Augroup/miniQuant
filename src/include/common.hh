#ifndef COMMON_H
#define COMMON_H

#include <string>
#include <vector>
#include <iostream>
#include <ggcat.hh>
enum ReadType{
    ShortReadType = 0,
    LongReadType = 1
};
enum ShortReadsStrandness{
  unstranded = 0,
  fr_stranded = 1,
  rf_stranded = 2,
};
enum LongReadsLibraryPrep{
  cDNA_ONT = 0,
  dRNA_ONT = 1,
  cDNA_PacBio = 2,
};
enum KvalueEntryType {
  kvalue_entry_type_effective_length = 0,
  kvalue_entry_type_binary = 1,
};
struct ProgramOptions {

  // quant options
  size_t threads;
  bool forward_only;
  bool debug;
  float mem;
  size_t short_reads_kmer_size;
  size_t long_reads_kmer_size;
  int minimizer_size;
  bool use_short_reads;
  bool use_long_reads;
  // Mean fragment length of short reads
  double fraglen;
  // Standard deviation of fragment length of short reads
  double fragsd;
  double long_reads_weight_sum;
  double short_reads_weight_sum;
  size_t num_long_reads;
  size_t num_short_reads;
  size_t num_long_reads_chunks;
  std::string transcriptsFile;
  std::string longreadsFile;
  std::string shortreadsFilePairOne;
  std::string shortreadsFilePairTwo;
  std::string tempFolder;
  std::string outputFolder;
  std::string inputWeightFile;
  std::string inputAbundanceFile;
  bool use_input_weight;
  bool use_input_abundance;
  ggcat::ExtraElaborationStep extra_elaboration_step;
  ShortReadsStrandness short_reads_strandness;
  LongReadsLibraryPrep long_reads_library_prep;

  // K-value options
  size_t minimum_split_exon_map_len;
  std::string annotationFile;
  bool not_normalize_entry;
  KvalueEntryType kvalue_entry_type;
  ProgramOptions():
    threads(1),
    short_reads_kmer_size(31),
    long_reads_kmer_size(31),
    minimizer_size(-1),
    mem(5.0),
    debug(false),
    forward_only(false),
    extra_elaboration_step(ggcat::ExtraElaborationStep_None),
    use_input_weight(false),
    use_long_reads(false),
    use_short_reads(false),
    fraglen(235.0),
    fragsd(23.0),
    num_long_reads(0),
    num_short_reads(0),
    long_reads_weight_sum(1.0),
    short_reads_weight_sum(1.0),
    short_reads_strandness(unstranded),
    long_reads_library_prep(cDNA_ONT),
    minimum_split_exon_map_len(1),
    not_normalize_entry(false),
    num_long_reads_chunks(0),
    kvalue_entry_type(kvalue_entry_type_effective_length)
    {}


};

#endif