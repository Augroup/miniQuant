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
#include <cmath>
#include <cassert>
#include <regionGenerator.hh>
#include <annotationParser.hh>
#include <util.hh>
#include <common.hh>
#include <fileWriter.hh>


int SHORT_READ_LENGTH = 235;
int READ_MIN_MAP_LEN = 1;


void RegionGenerator::traverseSplitExon(std::vector<SplitExon> &splittedExons,std::map<std::vector<int>,Region,SplitExonIndicesCompare> &regions_map,int index,int inner_length,int start_length,int end_length,int region_start,int region_end,Region prev_region){
    // The region cannot be covered by short read
    if (inner_length >  SHORT_READ_LENGTH){
        return;
    }
    std::set<int> visited;
    int region_eff_length = -1;
    // if only an exon
    if (inner_length == -1){
        region_eff_length = start_length - SHORT_READ_LENGTH + 1;
    } else {
        std::vector<int> region_eff_length_arr = {SHORT_READ_LENGTH - inner_length - 2*READ_MIN_MAP_LEN,start_length-READ_MIN_MAP_LEN,end_length-READ_MIN_MAP_LEN};
        region_eff_length = *std::min_element(region_eff_length_arr.begin(), region_eff_length_arr.end());
        region_eff_length += 1;
    }
    // std::cout <<region_eff_length << ";" << inner_length << ";" << start_length << ";" << end_length<< ";" << region_start << ";" << region_end << std::endl;
    // if (region_eff_length > 0) {
        // add all isoforms it maps to
    // create the region
    Region region(region_start,region_end,region_eff_length);
    region.split_exon_indices.push_back(index);
    if (inner_length == -1){
        for (auto &mapped_isoform : splittedExons[index].belong_isoform_id_set){
        // If this is the single exon region
            region.belong_isoform_id_set.insert(mapped_isoform);
            
        }
    } else {
        for (auto &mapped_isoform : splittedExons[index].belong_isoform_id_set){
            // If previous region has this isoform
            if (prev_region.belong_isoform_id_set.find(mapped_isoform) !=prev_region.belong_isoform_id_set.end()){
                region.belong_isoform_id_set.insert(mapped_isoform);
            }   
        }
    }
   
    if (region_eff_length > 0){
        std::vector<int> region_split_exon_indices = prev_region.split_exon_indices;
        region_split_exon_indices.push_back(index);
        regions_map[region_split_exon_indices] = region;
    }
    prev_region = region;

    for (auto & isoform_pair :splittedExons[index].isoform_split_exon_next_map){        
        int nextIndex = isoform_pair.second;
        if (nextIndex == -1){
            continue;
        } else {
            if (visited.count(nextIndex) > 0)
                continue;
            else
                visited.insert(nextIndex);
            int new_inner_length = 0;
            // if only an exon
            if (inner_length == -1){
                new_inner_length = 0;
            } else {
                // the last exon added to inner length
                new_inner_length = inner_length + end_length;
            }
            assert (new_inner_length>=0);
            // start should keep the same
            int new_start_length = start_length;
            // the end length is current exon
            int new_end_length = splittedExons[nextIndex].length;
            int new_region_start = region_start;
            // the new region end is the current exon end
            int new_region_end = splittedExons[nextIndex].end;
            traverseSplitExon(splittedExons,regions_map,nextIndex,new_inner_length,new_start_length,new_end_length,new_region_start,new_region_end,prev_region);
        }
    }

}
void RegionGenerator::splitExons(std::vector<SplitExon> &splittedExons,const std::vector<Isoform>& isoforms) {
    std::vector<Point> points;
    // Collect all exon start and end points
    for (const auto& isoform : isoforms) {
        for (const auto& exon : isoform.exons) {
            points.emplace_back(exon.start, true,isoform.id);
            points.emplace_back(exon.end, false,isoform.id);
        }
    }

    // Sort points: by position first, endpoints come after start points
    std::sort(points.begin(), points.end(), [](const Point& a, const Point& b) {
        if (a.pos == b.pos) return a.isStart < b.isStart;
        return a.pos < b.pos;
    });

    // int overlapCount = 0;
    // int currentStart = -1;
    size_t last_start = -1;
    std::set<std::string> current_isoform_ids_set;
    // Traverse points to build non-overlapping intervals
    for (const auto& point : points) {
        if (point.isStart) {
            if (last_start != -1) {
                if (last_start != point.pos && !current_isoform_ids_set.empty())
                    splittedExons.emplace_back(last_start, point.pos,current_isoform_ids_set);
            };
            current_isoform_ids_set.insert(point.isoform_id);
            last_start  = point.pos;
        } else {
            if (last_start != point.pos && !current_isoform_ids_set.empty())
                splittedExons.emplace_back(last_start, point.pos,current_isoform_ids_set);
            current_isoform_ids_set.erase(point.isoform_id);
            last_start = point.pos;
        }
    };
    if (last_start != -1){
        if (last_start != points.back().pos && !current_isoform_ids_set.empty())
            splittedExons.emplace_back(last_start,points.back().pos,current_isoform_ids_set);
    }

    
}
void RegionGenerator::generateRegion(std::vector<SplitExon> &splittedExons,std::vector<Region> &regions) {

    std::map<std::string,size_t> last_isoform_split_exon_index;
    for (size_t i = splittedExons.size()-1;;i--){
        for (auto &isoform :splittedExons[i].belong_isoform_id_set){
            if (last_isoform_split_exon_index.count(isoform)>0){
                splittedExons[i].isoform_split_exon_next_map[isoform] = last_isoform_split_exon_index[isoform];
            } else {
                // current is the last split exon
                splittedExons[i].isoform_split_exon_next_map[isoform] = -1;
            }
            last_isoform_split_exon_index[isoform] = i;   
        }
        if (i == 0) break;
    }
    std::map<std::vector<int>,Region,SplitExonIndicesCompare> regions_map;
    for (size_t i = 0;i<splittedExons.size();i++){
        traverseSplitExon(splittedExons,regions_map,i,-1,splittedExons[i].length,-1,splittedExons[i].start,splittedExons[i].end,Region());
        // break;
        // splittedExons[i].printSplitExon();
    }
    // std::cout << std::endl;
    for (auto & [_,region]: regions_map){
        if (region.eff_len > 0 && !region.belong_isoform_id_set.empty())
            regions.push_back(region);
    }
    // for (auto & region:regions){
    //     region.printRegion();
    // }
    // std::cout << "===================="<<std::endl;

}
void RegionGenerator::generateRegionForGene(Gene &gene,std::vector<Region> &regions,ProgramOptions &opt){
    SHORT_READ_LENGTH = std::round(opt.fraglen);
    READ_MIN_MAP_LEN = std::round(opt.minimum_split_exon_map_len);
    // for (auto & gene : genes){
    std::vector<SplitExon> splittedExons;
    std::vector<Isoform> isoforms;
    for (auto & isoform : gene.isoforms){
        isoforms.push_back(isoform.second);
    }
    splitExons(splittedExons,isoforms);
    // count number of split_exons for gene and isoform
    for (auto & splittedExon : splittedExons){
        for (auto & isoform_id : splittedExon.belong_isoform_id_set){
            gene.isoforms[isoform_id].num_split_exons += 1;
        }
    }
    gene.num_split_exons += splittedExons.size();
    if (opt.debug){
        FileWriter fileWriter;
        fileWriter.writeSplitExon(opt.outputFolder+"/split_exons_"+gene.id+".jsonl",splittedExons,opt);
    }
    
    generateRegion(splittedExons,regions);
    // count number of regions for gene and isoform
    for (auto & region : regions){
        for (auto & isoform_id : region.belong_isoform_id_set){
            gene.isoforms[isoform_id].num_regions += 1;
        }
    }
    gene.num_regions += regions.size();

    // }
}
