#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cxxopts.hpp>
#include <ggcat.hh>
#include <align.hh>
#include <EM.hh>
#include <common.hh>
#include <fileWriter.hh>
#include <calKvalue.hh>
#include <util.hh>
bool ParseKvalueOptions(int argc, char** argv,ProgramOptions& opt){
    cxxopts::Options options("miniQuant kvalue", "\nCalculate K-value to identify a problematic set of gene isoforms with erroneous quantification\n");
    try{
    options.add_options("Required")
        ("a,annotation", "Gene isoform annotation file in GTF, GFF or genePred format", cxxopts::value<std::string>());
    options.add_options("Optional")
        ("o,output", "Output folder", cxxopts::value<std::string>()->default_value("./miniQuant_kvalue/"))
        ("t,threads", "Num of threads", cxxopts::value<size_t>()->default_value("1"))
        ("short_reads_mean_fragment_length", "Mean value of short reads fragment lengths", cxxopts::value<double>()->default_value("235.0"))
        ("debug", "Whether output debugging information", cxxopts::value<bool>()->default_value("false"))
        ("not_normalize_entry", "Whether NOT normalize region-isoform matrix (A matrix) before calculating K-value", cxxopts::value<bool>()->default_value("false"))
        ("kvalue_entry_type", "What kind of entry to use for region-isoform matrix (A matrix). Choices:[effective_length,binary]", cxxopts::value<std::string>()->default_value("effective_length"))
        // ("short_reads_sd_fragment_length", "Standard deviation of short reads fragment lengths", cxxopts::value<double>()->default_value("23.0"))
        ("h,help", "Print usage");
    auto result = options.parse(argc, argv);
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return false;
    }
    if (!result.count("annotation")) {
        std::ostringstream error_msg;
        error_msg << "Gene isoform annotation file (-a,--annotation) is required for calculating the K-value! Exit...";
        program_log_message(Error,error_msg.str().c_str());
        return false;
    }
    opt.threads = result["threads"].as<size_t>();
    opt.annotationFile = result["annotation"].as<std::string>();
    opt.outputFolder = result["output"].as<std::string>();
    opt.fraglen = result["short_reads_mean_fragment_length"].as<double>();
    // opt.fragsd = result["short_reads_sd_fragment_length"].as<double>();
    opt.tempFolder = result["output"].as<std::string>()+"/temp";
    opt.debug =result["debug"].as<bool>();
    opt.not_normalize_entry = result["not_normalize_entry"].as<bool>();
    std::string kvalue_entry_type = result["kvalue_entry_type"].as<std::string>();
    if (kvalue_entry_type == "effective_length"){
        opt.kvalue_entry_type = kvalue_entry_type_effective_length;
    } else if (kvalue_entry_type == "binary"){
        opt.kvalue_entry_type = kvalue_entry_type_binary;
    } else{
        program_log_message(Error,"--kvalue_entry_type is invalid! Choices:[effective_length,binary]");
        std::cout << options.help() << std::endl;
        return false;
    };
    return true;
    }
    catch (...){
        std::cout << options.help() << std::endl;
        return false;
        
    }
}
bool ParseQuantOptions(int argc, char** argv,ProgramOptions& opt){

    cxxopts::Options options("miniQuant quant -r <ref.fa> -l <lr.fq.gz> -o <./output> [-1 <sr_r1.fq.gz> -2 <sr_r2.fq.gz>]\n[--long_reads_library_prep CHOICE] [--short_reads_strandness CHOICE] [-t THREADS]", "\nQuantify gene isoform abundance by long reads and/or short reads\n");
    try{
    options.add_options("Required")
        ("r,reference", "Reference transcripts sequence file in FASTA format", cxxopts::value<std::string>())
        ("l,long_reads", "Input long reads file in plain or gzipped FASTA/FASTQ format.", cxxopts::value<std::string>()->default_value(""))
        ("1,short_reads_pair_1", "Input short reads pair 1 in plain or gzipped FASTA/FASTQ format. Leave blank if using only long reads.", cxxopts::value<std::string>()->default_value(""))
        ("2,short_reads_pair_2", "Input short reads pair 2 in plain or gzipped FASTA/FASTQ format. Leave blank if using only long reads.", cxxopts::value<std::string>()->default_value(""));
     options.add_options("Optional")
        ("o,output", "Output folder", cxxopts::value<std::string>()->default_value("./miniQuant_res/"))
        // ("temp", "Temp folder", cxxopts::value<std::string>()->default_value("/temp"))
        // ("debug", "Whether output debugging information", cxxopts::value<bool>()->default_value("false"))
        // ("k,kmer_size", "K-mer size", cxxopts::value<size_t>()->default_value("31"))
        ("short_reads_kmer_size", "Short reads K-mer size", cxxopts::value<size_t>()->default_value("39"))
        ("long_reads_kmer_size", "Long reads K-mer size", cxxopts::value<size_t>()->default_value("40"))
        // ("minimizer_size", "Minimizer size", cxxopts::value<int>()->default_value("-1"))
        // ("extra_elaboration_step", "Extra elaboration step", cxxopts::value<std::string>()->default_value("none"))
        ("long_reads_library_prep", "The library preparation for long reads. Choices:[cDNA-ONT,dRNA-ONT,cDNA-PacBio]", cxxopts::value<std::string>()->default_value("cDNA-ONT"))
        ("short_reads_strandness", "The strandness of short reads. Choices:[unstranded,fr-stranded,rf-stranded]\n\n*fr-stranded: Strand specific reads, first read forward\n*rf-stranded: Strand specific reads, first read reverse\n", cxxopts::value<std::string>()->default_value("unstranded"))
        // ("input_weight", "The input weight for debugging", cxxopts::value<std::string>()->default_value(""))
        // ("input_abundance", "The input true abundance for debugging", cxxopts::value<std::string>()->default_value(""))
        // ("long_reads_weight", "The weight for long reads", cxxopts::value<double>()->default_value("1.0"))
        // ("short_reads_weight", "The weight for short reads", cxxopts::value<double>()->default_value("1.0"))
        ("short_reads_mean_fragment_length", "Mean value of short reads fragment lengths", cxxopts::value<double>()->default_value("235.0"))
        ("short_reads_sd_fragment_length", "Standard deviation of short reads fragment lengths", cxxopts::value<double>()->default_value("23.0"))
        ("t,threads", "Num of threads", cxxopts::value<size_t>()->default_value("1"))
        ("mem", "Max RAM usage in GB allowed when aligning the reads", cxxopts::value<double>()->default_value("20.0"))
        ("h,help", "Print usage");
    options.parse_positional({"subcommand"});
    auto result = options.parse(argc, argv);
    // required arguments
    if (result.count("help")) {
        std::cout << options.help() << std::endl;
        return false;
    }
    if (!result.count("reference")) {
        std::ostringstream error_msg;
        error_msg << "Reference transcripts file (-r,--reference) is required for quantification! Exit...";
        program_log_message(Error,error_msg.str().c_str());
        return false;
    }
    // || !result.count("output") || !result.count("reference")
    if (!(result.count("long_reads") || (result.count("short_reads_pair_1") && result.count("short_reads_pair_2")))) {
        std::ostringstream error_msg;
        error_msg << "Either short reads (-1,--short_reads_pair_1 and -2,--short_reads_pair_2) or long reads (-l,--long_reads) need to be provided! Exit...";
        program_log_message(Error,error_msg.str().c_str());
        return false;
    }
    if ((result.count("short_reads_pair_1") && !result.count("short_reads_pair_2")) || (result.count("short_reads_pair_2") && !result.count("short_reads_pair_1"))) {
        std::ostringstream error_msg;
        error_msg << "Both pairs for short reads (-1,--short_reads_pair_1  -2,--short_reads_pair_2) need to be provided! Exit...";
        program_log_message(Error,error_msg.str().c_str());
        return false;
    }

    opt.transcriptsFile = result["reference"].as<std::string>();
    if (!isFile(opt.transcriptsFile) || !(isReferenceSequenceFile(opt.transcriptsFile))){
        std::ostringstream error_msg;
        error_msg << "Reference transcripts file "<< opt.transcriptsFile << " doesn't exit or in invalid format(plain .fasta,.fa only)! Exit...";
        program_log_message(Error,error_msg.str().c_str());
        return false;
    }
    opt.use_long_reads = false;
    opt.use_short_reads = false;

    if (result.count("long_reads")){
        opt.longreadsFile = result["long_reads"].as<std::string>();
        if (!isFile(opt.longreadsFile) || !(isReadSequenceFile(opt.longreadsFile))){
            std::ostringstream error_msg;
            error_msg << "Long reads sequences "<< opt.longreadsFile << " doesn't exit or in invalid format(plain or gzipped FASTA/FASTQ only)! Exit...";
            program_log_message(Error,error_msg.str().c_str());
            return false;
        }
        opt.use_long_reads = true;
    }
    if (result.count("short_reads_pair_1") && result.count("short_reads_pair_2")){
        opt.shortreadsFilePairOne = result["short_reads_pair_1"].as<std::string>();
        opt.shortreadsFilePairTwo = result["short_reads_pair_2"].as<std::string>();
        if (!isFile(opt.shortreadsFilePairOne) || !isFile(opt.shortreadsFilePairTwo) || !(isReadSequenceFile(opt.shortreadsFilePairOne) || !(isReadSequenceFile(opt.shortreadsFilePairTwo)))){
            std::ostringstream error_msg;
            error_msg << "Short reads sequences "<< opt.shortreadsFilePairOne << " and/or " << opt.shortreadsFilePairTwo << " don't exit or in invalid format(plain or gzipped FASTA/FASTQ only)! Exit...";
            program_log_message(Error,error_msg.str().c_str());
            return false;
        }
        opt.use_short_reads = true;
    }

    // opt.tempFolder = result["output"].as<std::string>()+"/"+result["temp"].as<std::string>();
    opt.tempFolder = result["output"].as<std::string>()+"/temp";
    opt.outputFolder = result["output"].as<std::string>();
    // opt.debug =result["debug"].as<bool>();
    opt.debug = false;
    
    // // set weight for short and long reads
    // opt.long_reads_weight = result["long_reads_weight"].as<double>();
    // opt.short_reads_weight = result["short_reads_weight"].as<double>();
    // // TODO: remember checking whether this holds, when weight = 0, do not use the reads
    // if (opt.long_reads_weight == 0){
    //     opt.use_long_reads = false;
    // }
    // if (opt.short_reads_weight == 0){
    //     opt.use_short_reads = false;
    // }
    opt.fraglen = result["short_reads_mean_fragment_length"].as<double>();
    opt.fragsd = result["short_reads_sd_fragment_length"].as<double>();
    opt.minimizer_size = -1;
    opt.threads = result["threads"].as<size_t>();
    opt.threads = smallestPowerOfTwoLargerOrEqual(opt.threads);
    opt.mem = result["mem"].as<double>();
    
    // opt.minimizer_size = result["minimizer_size"].as<int>();
    // Short reads strandness
    std::string short_reads_strandness = result["short_reads_strandness"].as<std::string>();
    if (short_reads_strandness == "unstranded"){
        opt.short_reads_strandness = unstranded;
    } else if (short_reads_strandness == "fr-stranded"){
        opt.short_reads_strandness = fr_stranded;
    } else if (short_reads_strandness == "rf-stranded"){
        opt.short_reads_strandness = rf_stranded;
    } else {
        program_log_message(Error,"--short_reads_strandness is invalid!");
        std::cout << options.help() << std::endl;
        return false;
    };
    // set library prep for short and long reads
    std::string long_reads_library_prep = result["long_reads_library_prep"].as<std::string>();
    if (long_reads_library_prep == "cDNA-ONT"){
        opt.long_reads_kmer_size = 40;
        opt.long_reads_library_prep = cDNA_ONT;
    } else if (long_reads_library_prep == "dRNA-ONT"){
        opt.long_reads_kmer_size = 40;
        opt.long_reads_library_prep = dRNA_ONT;
    } else if (long_reads_library_prep == "cDNA-PacBio"){
        opt.long_reads_kmer_size = 40;
        opt.long_reads_library_prep = cDNA_PacBio;
    } else {
        program_log_message(Error,"--long_reads_library_prep is invalid! Choices:[cDNA-ONT,dRNA-ONT,cDNA-PacBio]");
        std::cout << options.help() << std::endl;
        return false;
    };
    if (result.count("long_reads_kmer_size")) {
        opt.long_reads_kmer_size = result["long_reads_kmer_size"].as<size_t>();
    }
    if (result.count("short_reads_kmer_size")) {
        opt.short_reads_kmer_size = result["short_reads_kmer_size"].as<size_t>();
    }
    // opt.inputWeightFile = result["input_weight"].as<std::string>();
    // if (!isFile(opt.inputWeightFile) || !(has_suffix(opt.inputWeightFile,".tsv"))){
    //     opt.use_input_weight = false;
    // } else {
    opt.use_input_weight = false;
    // }
    // opt.inputAbundanceFile = result["input_abundance"].as<std::string>();
    // if (!isFile(opt.inputAbundanceFile) || !(has_suffix(opt.inputAbundanceFile,".tsv"))){
    opt.use_input_abundance = false;
    // } else {
    //     opt.use_input_abundance = true;
    // }

    // TODO: Using none as default. No impact for quantification
    // std::string extra_elaboration_step = result["extra_elaboration_step"].as<std::string>();
    opt.extra_elaboration_step = ggcat::ExtraElaborationStep_None;
    // if (extra_elaboration_step=="none"){
    //     opt.extra_elaboration_step = ggcat::ExtraElaborationStep_None;
    // } else if (extra_elaboration_step == "maximal"){
    //     opt.extra_elaboration_step = ggcat::ExtraElaborationStep_UnitigLinks;
    // } else if (extra_elaboration_step == "euler"){
    //     opt.extra_elaboration_step = ggcat::ExtraElaborationStep_Eulertigs;
    // } else if (extra_elaboration_step == "path"){
    //     opt.extra_elaboration_step = ggcat::ExtraElaborationStep_Pathtigs;
    // } else {
    //     std::cout << options.help() << std::endl;
    //     return false;
    // }
    return true;
    }
    catch (...){
        std::cout << options.help() << std::endl;
        return false;
        
    }

}
bool ParseOptions(int argc, char** argv,std::string &subcommand,ProgramOptions& opt){
    if ((argc > 1) && (!strcmp(argv[1], "quant") || !strcmp(argv[1],"kvalue") || !strcmp(argv[1],"train"))) {
        subcommand = argv[1];
        return true;
    } else {
        std::cout << std::endl;
        std::cout << "Program: miniQuant (Optimal and fast gene isoform abundance estimation)" << std::endl;
        std::cout << "Version: 1.4.1" << std::endl<< std::endl;
        std::cout << "Usage: miniQuant <subcommand> [arguments] .." << std::endl<< std::endl;
        std::cout << "Available subcommands:" << std::endl<< std::endl;
        std::cout << "  quant       Quantify gene isoform abundance by long reads and/or short reads" << std::endl;
        std::cout << "  kvalue      Calculate K-value to identify a problematic set of gene isoforms with erroneous quantification" << std::endl<< std::endl;
        // Perhaps list available subcommands here
        return false;
        
    }
}

int execQuant(ProgramOptions &opt){
    emptyAndCreateFolder(opt.outputFolder);
    emptyAndCreateFolder(opt.tempFolder);

    std::vector<Transcript> transcripts;
    runAlign(transcripts,opt);
    std::vector<Read> reads;
    runEM(std::ref(reads),transcripts,opt);            
    if (!opt.debug){
        removeFolder(opt.tempFolder);
    }
    return 0;
}
int execKvalue(ProgramOptions &opt){
    emptyAndCreateFolder(opt.outputFolder);
    // emptyAndCreateFolder(opt.tempFolder);
    parallelCalKvalue(opt);
    return 0;
}

int main(int argc, char** argv) {
    ProgramOptions opt;
    // try{
    std::string subcommand;
    if (ParseOptions(argc,argv,subcommand,opt)){
        if (subcommand == "quant"){
            // below for quantification
            if (ParseQuantOptions(argc-1, argv+1,opt)){
                return execQuant(opt);
            }
        } else if (subcommand == "kvalue"){
            if (ParseKvalueOptions(argc-1, argv+1,opt)){
                return execKvalue(opt);
            }
        } else if (subcommand == "train"){
            if (ParseQuantOptions(argc-1, argv+1,opt)){
                opt.debug = true;
                return execQuant(opt);
            }
        }
    }

    // } catch(...){
    //     return 1;
    // }
    return 1;
       
}