#ifndef FILE_WRITER_H
#define FILE_WRITER_H

#include <vector>
#include <string>
#include <EM.hh>
#include <Eigen/Dense>
#include <regionGenerator.hh>
class FileWriter {
public:

    // Method to write abundance to a TSV file
    void writeAbundance(const std::string& filename,std::vector<Transcript>& transcripts,ProgramOptions &opt);
    // Method to write kvalue to a TSV file
    void writeKvalue(const std::string& filename,std::vector<Gene> &genes,ProgramOptions &opt);
    void writeDesignMatrix(const std::string& filename,std::string gene_name,std::map<std::string, Isoform> isoforms,Eigen::MatrixXd regionMatrix,ProgramOptions &opt);
    void writeIsoformKvalue(const std::string& filename,std::vector<Gene> &genes,ProgramOptions &opt);

    // Method to write region to a TXT file
    void writeRegion(const std::string& filename,std::vector<Region> &regions,ProgramOptions &opt);
    void writeSplitExon(const std::string& filename,std::vector<SplitExon> &splitExons,ProgramOptions &opt);
    // Method to write transcript to a JSONL file
    void writeTranscript(const std::string& filename,std::vector<Transcript>& transcripts);

    // Method to write read to a JSONL file
    void writeRead(const std::string& filename,std::vector<Read>& reads);
    void readWeight(const std::string& filename,std::vector<Read>& reads,ProgramOptions &opt);
    void readAbundance(const std::string& filename,std::vector<Transcript>& transcripts);
};

#endif