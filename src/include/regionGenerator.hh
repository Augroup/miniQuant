#ifndef REGION_GENERATOR_HH
#define REGION_GENERATOR_HH

#include <iostream>
#include <vector>
#include <map>
#include <string>
#include <set>
#include <common.hh>
#include <sstream>
#include <annotationParser.hh>

struct SplitExon {
    size_t start;
    size_t end;
    size_t length;
    std::set<std::string> belong_isoform_id_set;
    std::map<std::string,int> isoform_split_exon_next_map;
    SplitExon()=default;
    SplitExon(size_t s, size_t e,std::set<std::string> set_id) : start(s), end(e),belong_isoform_id_set(set_id),length(end-start){}

    void printSplitExon() const {
        std::cout << "(" << start << "," << end <<")"<<"\n";
         for (auto & isoform : belong_isoform_id_set){
            std::cout << isoform <<";";
        }
        std::cout << std::endl;
    }
    std::string toJSON() const {
        std::ostringstream outs;
        bool firstIsoform = true;
        outs << "{\"split_exon_start\":" << start << ",\"split_exon_end\":" << end << ",\"eff_len\":" << length << ",\"isoforms\":[";
        for (auto & isoform : belong_isoform_id_set){
            if (firstIsoform)
                outs << "\"" << isoform << "\"";
            else
                outs << ",\"" << isoform << "\"";
            firstIsoform = false;
        }
        outs << "]}\n";
        return outs.str();
    }
};
struct SplitExonIndicesCompare {
    bool operator()(const std::vector<int>& lhs, const std::vector<int>& rhs) const {
        return lhs < rhs; // lexicographical comparison
    }
};
struct Region {
    size_t start;
    size_t end;
    int eff_len;
    std::vector<int> split_exon_indices;
    std::set<std::string> belong_isoform_id_set;
    Region():eff_len(-1){}
    Region(size_t s, size_t e,int eff_len) : start(s), end(e),eff_len(eff_len),split_exon_indices({}){}
    void printRegion() const {
        std::cout << "(" << start << "," << end <<","<<eff_len<<")"<<"\n";
        for (auto & isoform : belong_isoform_id_set){
            std::cout << isoform <<";";
        }
        std::cout << std::endl << std::endl;
    }
    std::string toJSON() const {
        std::ostringstream outs;
        bool firstIsoform = true;
        outs << "{\"region_start\":" << start << ",\"region_end\":" << end << ",\"eff_len\":" << eff_len << ",\"isoforms\":[";
        for (auto & isoform : belong_isoform_id_set){
            if (firstIsoform)
                outs << "\"" << isoform << "\"";
            else
                outs << ",\"" << isoform << "\"";
            firstIsoform = false;
        }
        outs << "]}\n";
        return outs.str();
    }
};
struct Point {
    int pos;
    bool isStart;
    std::string isoform_id;
    Point()=default;
    Point(int p, bool start,std::string id) : pos(p), isStart(start),isoform_id(id) {}
};
class RegionGenerator{
public:
    void traverseSplitExon(std::vector<SplitExon> &splittedExons,std::map<std::vector<int>,Region,SplitExonIndicesCompare> &regions_map,int index,int inner_length,int start_length,int end_length,int region_start,int region_end,Region prev_region);
    void splitExons(std::vector<SplitExon> &splittedExons,const std::vector<Isoform>& isoforms);
    void generateRegion(std::vector<SplitExon> &splittedExons,std::vector<Region> &regions);
    void generateRegionForGene(Gene &gene,std::vector<Region> &regions,ProgramOptions &opt);
};
#endif