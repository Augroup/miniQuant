#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <set>
#include <map>
#include <regex>
#include <algorithm>
#include <utility>
#include <annotationParser.hh>
#include <util.hh>
#include <common.hh>
// All the coordinate stored will be 0-based. So convert 1-based gtf and gff3 to 0-based.

std::vector<std::string> split(const std::string& line, char delimiter) {
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(line);
    while (std::getline(tokenStream, token, delimiter)) {
        tokens.push_back(token);
    }
    return tokens;
}

std::string extractGTFAttribute(const std::string& attributes, const std::string& key) {
    size_t start = attributes.find(key + " \"");
    if (start == std::string::npos) {
        return "";
    }

    start += key.length() + 2;  // Skip past key=" opening
    size_t end = attributes.find("\"", start);
    return attributes.substr(start, end - start);
}

std::string extractGFFAttribute(const std::string& attributes, const std::string& key) {
    size_t start = attributes.find(key + "=");
    if (start == std::string::npos) {
        return "";
    }

    start += key.length() + 1;  // Skip past key= opening
    size_t end = attributes.find(";", start);
    if (end == std::string::npos) {
        end = attributes.length();  // If not found, take until the end of string
    }
    return attributes.substr(start, end - start);
}

void AnnotationParser::parseGTForGFF(std::ifstream& file, std::map<std::string, Gene>& genes,bool isGTF) {
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty() || line[0] == '#') {
            continue;  // Skip empty lines and comments
        }

        std::vector<std::string> columns = split(line, '\t');
        if (columns.size() < 8) {
            std::ostringstream error_msg;
            error_msg << "Invalid line encountered in annotation file: " << line;
            program_log_message(Error,error_msg.str().c_str());
            continue;
        }

        std::string seqname = columns[0];
        std::string source = columns[1];
        std::string feature = columns[2];
        if (feature != "exon")
            continue;
        int start = std::stoi(columns[3]);
        // change start from 1-based to 0-based
        start -= 1;
        int end = std::stoi(columns[4]);
        if (start >= end){
            std::ostringstream error_msg;
            error_msg << "Invalid line encountered in annotation file: " << line;
            program_log_message(Error,error_msg.str().c_str());
            continue;
        }

        std::string score = columns[5];
        char strand = columns[6][0];
        std::string frame = columns[7];

        if (columns.size() > 8) {
            std::string attributes = columns[8];
            
            std::string transcript_id;
            std::string gene_id;

            if (isGTF) {
                transcript_id = extractGTFAttribute(attributes, "transcript_id");
                gene_id = extractGTFAttribute(attributes, "gene_id");
            } else {
                transcript_id = extractGFFAttribute(attributes, "transcript_id");
                gene_id = extractGFFAttribute(attributes, "gene_id");
            }
            if (!gene_id.empty()) {
                if (genes.find(gene_id) == genes.end())
                    genes[gene_id] = Gene(gene_id);
                if (!transcript_id.empty()) {
                    genes[gene_id].addIsoform(transcript_id, start, end);
                }
            }

        }
    }

    file.close();
}

// Parsing logic for genePred
void AnnotationParser::parseGenePredFile(std::ifstream& file, std::map<std::string, Gene>& genes) {
    std::string line;
    while (std::getline(file, line)) {
        if (line.empty()) continue;

        std::istringstream ss(line);
        std::regex re("\t");

        std::vector<std::string> fields(
            std::sregex_token_iterator(line.begin(), line.end(), re, -1),
            std::sregex_token_iterator()
        );
        std::string geneId = fields[0];
        std::string transcriptId = fields[1];
        std::string chrom = fields[2];
        std::string strand = fields[3];
        int txStart = std::stoi(fields[4]);
        int txEnd = std::stoi(fields[5]);
        int exonCount = std::stoi(fields[8]);
        std::string exonStarts = fields[9];
        std::string exonEnds = fields[10];
        

        std::vector<int> starts, ends;
        std::stringstream startsStream(exonStarts);
        std::stringstream endsStream(exonEnds);

        // Parse exon start and end positions
        for (int i = 0; i < exonCount; ++i) {
            std::getline(startsStream, exonStarts, ',');
            std::getline(endsStream, exonEnds, ',');
            starts.push_back(std::stoi(exonStarts));
            ends.push_back(std::stoi(exonEnds));
        }



        if (genes.find(geneId) == genes.end()) {
            genes[geneId] = Gene(geneId);
        }

        for (int i = 0; i < exonCount; ++i) {
            if (starts[i] >= ends[i]){
                std::ostringstream error_msg;
                error_msg << "Invalid line encountered in annotation file: " << line;
                program_log_message(Error,error_msg.str().c_str());
                break;
            }
            genes[geneId].addIsoform(transcriptId, starts[i], ends[i]);
        }
    }
}
void AnnotationParser::decideAndParseFile(const std::string& fileName,std::map<std::string, Gene> &genes) {
    
    if (!isFile(fileName)) {
        program_log_message(UnrecoverableError,"Gene isoform annotation file not exist!");
    }
    std::ifstream file(fileName);
    if (has_suffix(fileName,".gpd")||has_suffix(fileName,".genepred"))
        parseGenePredFile(file, genes);
    else if (has_suffix(fileName,".gtf"))
        parseGTForGFF(file,genes,true);
    else if (has_suffix(fileName,".gff")||has_suffix(fileName,".gff3"))
        parseGTForGFF(file,genes,false);
    else{
        std::ostringstream error_msg;
        error_msg << "Gene isoform annotation file in unsupported format! Exit now.";
        program_log_message(UnrecoverableError,error_msg.str().c_str());
    }

}