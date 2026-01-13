
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <fileWriter.hh>
#include <EM.hh>
#include <regex>
#include <cmath>
#include <util.hh>
#include <annotationParser.hh>
#include <regionGenerator.hh>

// Method to write data to a TSV file
void FileWriter::writeAbundance(const std::string& filename,std::vector<Transcript>& transcripts,ProgramOptions &opt) {

    // Open file for writing
    std::ofstream outFile(filename);

    // Check if the file was opened successfully
    if (!outFile.is_open()) {
        std::ostringstream error_msg;
        error_msg << "Output folder " << opt.outputFolder << " not exist! Exit now!";
        program_log_message(UnrecoverableError,error_msg.str().c_str());
        
        return;
    }

    if (opt.use_short_reads && opt.use_long_reads)
        outFile << "Transcript_id\tTPM\tExpected_num_long_reads\tExpected_num_short_read_pairs\tEffective_length\n";
    if (opt.use_short_reads && !opt.use_long_reads)
        outFile << "Transcript_id\tTPM\tExpected_num_short_read_pairs\tEffective_length\n";
    if (!opt.use_short_reads && opt.use_long_reads)
        outFile << "Transcript_id\tTPM\tExpected_num_long_reads\n";

    for (size_t i = 0; i < transcripts.size(); ++i) {
        double num_expected_long_reads = opt.num_long_reads * transcripts[i].relative_abund;
        double num_expected_short_reads = opt.num_short_reads * transcripts[i].short_reads_effective_count;
        if (std::isnan(num_expected_long_reads) || num_expected_long_reads < 0){
            num_expected_long_reads = 0;
        }
        if (std::isnan(num_expected_short_reads) || num_expected_short_reads < 0){
            num_expected_short_reads = 0;
        }
        outFile << transcripts[i].name << '\t' << transcripts[i].abundance;
        if (opt.use_long_reads)
            outFile << '\t' << num_expected_long_reads;
        if (opt.use_short_reads)
            outFile << '\t' << num_expected_short_reads << '\t' << transcripts[i].short_reads_effective_length;
        
        outFile << '\n';
    }
    
    // Write data to file


    // Close the file
    outFile.close();
    
    // if (outFile) {
    //     std::cout << "Data successfully written to " << filename << std::endl;
    // }
}
void FileWriter::writeKvalue(const std::string& filename,std::vector<Gene> &genes,ProgramOptions &opt) {
    // Open file for writing
    std::ofstream outFile(filename);

    // Check if the file was opened successfully
    if (!outFile.is_open()) {
        std::ostringstream error_msg;
        error_msg << "Output folder " << opt.outputFolder << " not exist! Exit now!";
        program_log_message(UnrecoverableError,error_msg.str().c_str());
        
        return;
    };
    outFile << "Gene_id\tK-value\n";
    for (auto &gene : genes){
        std::string kval_str;
        if (gene.kvalue >0){
            kval_str = std::to_string(gene.kvalue);
            // if (gene.isFullRank){
            //     kval_str = std::to_string(gene.kvalue) + "\tFullRank";
            // } else {
            //     kval_str = std::to_string(gene.kvalue) + "\tNonFullRank";
            // }
        } else {
            kval_str = "NA";
        }
        outFile << gene.id << '\t' << kval_str <<'\n';
    };

    outFile.close();

}
void FileWriter::writeDesignMatrix(const std::string& filename,std::string gene_name,std::map<std::string, Isoform> isoforms,Eigen::MatrixXd regionMatrix,ProgramOptions &opt) {
    // Open file for writing
    std::ofstream outFile(filename);
    Eigen::IOFormat CleanFmt(Eigen::FullPrecision, 0, ", ", "\n", "", "");

    // Check if the file was opened successfully
    if (!outFile.is_open()) {
        std::ostringstream error_msg;
        error_msg << "Output folder " << opt.outputFolder << " not exist! Exit now!";
        program_log_message(UnrecoverableError,error_msg.str().c_str());
        
        return;
    };
    outFile << regionMatrix.format(CleanFmt) << "\n";
    outFile << gene_name << "\n";
    for (auto &[isoform_id,_]:isoforms){
        outFile << isoform_id << ",";
    }
    outFile.close();

}
void FileWriter::writeIsoformKvalue(const std::string& filename,std::vector<Gene> &genes,ProgramOptions &opt) {
    // Open file for writing
    std::ofstream outFile(filename);

    // Check if the file was opened successfully
    if (!outFile.is_open()) {
        std::ostringstream error_msg;
        error_msg << "Output folder " << opt.outputFolder << " not exist! Exit now!";
        program_log_message(UnrecoverableError,error_msg.str().c_str());
        
        return;
    };
    outFile << "Transcript_id\tGene_id\tIsoform_length\tIsoform_num_exons\tIsoform_num_split_exons\tIsoform_num_regions\tGene_K-value\tGene_isFullRank\tGene_num_isoforms\tGene_num_split_exons\tGene_num_regions\n";
    for (auto &gene : genes){
        std::string kval_str;
        if (gene.kvalue >0){
            kval_str = std::to_string(gene.kvalue);
        } else {
            kval_str = "NA";
        }
        for (auto [isoform_id,isoform]:gene.isoforms){
            outFile << isoform_id << '\t' << gene.id << '\t' \
            << isoform.isoform_length <<'\t' <<isoform.exons.size()<<'\t' << isoform.num_split_exons << '\t' << isoform.num_regions << '\t' \
            <<kval_str << '\t' << gene.isFullRank << '\t' << gene.isoforms.size() << '\t' << gene.num_split_exons << '\t' << gene.num_regions<<'\n';
        }
        
    };

    outFile.close();

}
void FileWriter::writeRegion(const std::string& filename,std::vector<Region> &regions,ProgramOptions &opt) {
    // Open file for writing
    std::ofstream outFile(filename);

    // Check if the file was opened successfully
    if (!outFile.is_open()) {
        std::ostringstream error_msg;
        error_msg << "Output folder " << opt.outputFolder << " not exist! Exit now!";
        program_log_message(UnrecoverableError,error_msg.str().c_str());
        
        return;
    };
    for (auto &region : regions){
        outFile << region.toJSON();
    };

    outFile.close();

}
void FileWriter::writeSplitExon(const std::string& filename,std::vector<SplitExon> &splitExons,ProgramOptions &opt) {
    // Open file for writing
    std::ofstream outFile(filename);

    // Check if the file was opened successfully
    if (!outFile.is_open()) {
        std::ostringstream error_msg;
        error_msg << "Output folder " << opt.outputFolder << " not exist! Exit now!";
        program_log_message(UnrecoverableError,error_msg.str().c_str());
        
        return;
    };
    for (auto &splitExon : splitExons){
        outFile << splitExon.toJSON();
    };

    outFile.close();

}

void FileWriter::writeRead(const std::string& filename,std::vector<Read>& reads) {

    // Open file for writing
    std::ofstream outFile(filename);

    // Check if the file was opened successfully
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    for (auto &read : reads){
        outFile << read.toJSONL();
    }

    // Close the file
    outFile.close();
}
void FileWriter::writeTranscript(const std::string& filename,std::vector<Transcript>& transcripts){

    // Open file for writing
    std::ofstream outFile(filename);

    // Check if the file was opened successfully
    if (!outFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    for (auto &transcript : transcripts){
        outFile << transcript.toJSONL();
    }

    // Close the file
    outFile.close();

}
void FileWriter::readWeight(const std::string& filename,std::vector<Read>& reads,ProgramOptions &opt){
    std::map<int,double> short_read_weight;
    std::map<int,double> long_read_weight;
    opt.long_reads_weight_sum = 0;
    opt.short_reads_weight_sum = 0;
    double abundance_sum;
    // Open file for writing
    std::ifstream inFile(filename);

    // Check if the file was opened successfully
    if (!inFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    std::string line;
    while (std::getline(inFile, line)){
        std::regex re("\t");

        std::vector<std::string> fields(
            std::sregex_token_iterator(line.begin(), line.end(), re, -1),
            std::sregex_token_iterator()
        );
        try{
            double weight = std::stod(fields[2]);
            int readType_int = std::stoi(fields[1]);
            int read_id = std::stoi(fields[0]);
            if (readType_int == 0){
                // SR
                short_read_weight[read_id] = weight;
            } else {
                long_read_weight[read_id] = weight;

            }
        }
        catch (...){
        }
    }

    // Close the file
    inFile.close();
    // for (size_t i =0;i<reads.size();i++){
    //     if (reads[i].readType == ShortReadType){
    //         if (short_read_weight.count(reads[i].id) > 0)
    //             reads[i].weight = short_read_weight[reads[i].id];
    //         else
    //             reads[i].weight = 1;
    //         opt.short_reads_weight_sum += reads[i].weight;
    //     } else {
    //         if (long_read_weight.count(reads[i].id) > 0)
    //             reads[i].weight = long_read_weight[reads[i].id];
    //         else
    //             reads[i].weight = 1;
    //         opt.long_reads_weight_sum += reads[i].weight;
    //     }
    // }

}
void FileWriter::readAbundance(const std::string& filename,std::vector<Transcript>& transcripts){
    std::unordered_map<std::string,double> transcript_name_abundance;
    double abundance_sum;
    // Open file for writing
    std::ifstream inFile(filename);

    // Check if the file was opened successfully
    if (!inFile.is_open()) {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return;
    }
    std::string line;
    while (std::getline(inFile, line)){
        std::regex re("\t");

        std::vector<std::string> fields(
            std::sregex_token_iterator(line.begin(), line.end(), re, -1),
            std::sregex_token_iterator()
        );
        try{
            double TPM = std::stod(fields[1]);
            transcript_name_abundance[fields[0]] = TPM;
            abundance_sum += TPM;
        }
        catch (...){
        }
    }

    // Close the file
    inFile.close();
    if (abundance_sum > 0){
        double sumEffectiveCount = 0.0;
        for (auto & transcript : transcripts){
            if (transcript_name_abundance.count(transcript.name)>0){
                transcript.relative_abund = transcript_name_abundance[transcript.name]/abundance_sum;
                transcript.short_reads_effective_count = transcript.relative_abund * transcript.short_reads_effective_length;
                sumEffectiveCount += transcript.short_reads_effective_count;
            }
        }
        if (sumEffectiveCount > 0.0) {
            for (auto& transcript : transcripts) {
                transcript.short_reads_effective_count /= sumEffectiveCount;
            }
        }

    }
}