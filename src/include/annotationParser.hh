#ifndef ANNOTATION_PARSER_HH
#define ANNOTATION_PARSER_HH

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <set>

struct Exon {
    int start;
    int end;
    Exon()=default;
    Exon(int s, int e) : start(s), end(e) {}

    void printExon() const {
        std::cout << "Exon Start: " << start << ", End: " << end << std::endl;
    }
};

struct Isoform {
    std::string id;
    std::vector<Exon> exons;
    size_t num_split_exons;
    size_t num_regions;
    size_t isoform_length;

    Isoform()=default;
    Isoform(const std::string& id) : id(id),num_split_exons(0),num_regions(0),isoform_length(0) {}

    void addExon(int start, int end) {
        exons.emplace_back(start, end);
        // start and end are 0-based
        isoform_length += end - start;
    }

    void printIsoform() const {
        std::cout << "Isoform ID: " << id << std::endl;
        for (const auto& exon : exons) {
            exon.printExon();
        }
    }
};

struct Gene {
    std::string id;
    std::map<std::string, Isoform> isoforms;
    double kvalue;
    size_t num_regions;
    size_t num_split_exons;
    bool isFullRank;
    Gene()=default;
    Gene(const std::string& id) : id(id),kvalue(-1),isFullRank(false),num_split_exons(0),num_regions(0) {}

    void addIsoform(const std::string& isoformId, int exonStart, int exonEnd) {
        if (isoforms.find(isoformId) == isoforms.end()) {
            isoforms[isoformId] = Isoform(isoformId);
        }
        isoforms[isoformId].addExon(exonStart, exonEnd);
    }

    void printGene() const {
        std::cout << "Gene ID: " << id << std::endl;
        for (const auto& isoformPair : isoforms) {
            isoformPair.second.printIsoform();
        }
    }
};

class AnnotationParser{
public:
    // Parsing functions
    void parseGTForGFF(std::ifstream& file, std::map<std::string, Gene>& genes,bool isGTF);
    void parseGenePredFile(std::ifstream& file, std::map<std::string, Gene>& genes);
    void decideAndParseFile(const std::string& fileName,std::map<std::string, Gene>& genes);
};



#endif