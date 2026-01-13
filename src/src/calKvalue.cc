
#include <Eigen/Dense>
#include <thread>
#include <calKvalue.hh>
#include <fileWriter.hh>
#include <annotationParser.hh>
#include <regionGenerator.hh>
#include <common.hh>
#include <util.hh>
Eigen::MatrixXd normalizeCols(const Eigen::MatrixXd& matrix) {
    Eigen::MatrixXd normalizedMatrix = matrix;
    for (int j = 0; j < matrix.cols(); ++j) {
        double colSum = matrix.col(j).sum();
        if (colSum != 0) {
            normalizedMatrix.col(j) /= colSum;
        } else {
            // Handle the case where the column sum is zero 
            normalizedMatrix.col(j).setConstant(0);
        }
    }
    return normalizedMatrix;
}
void createRegionMatrix(std::vector<Region> &regions,std::map<std::string, Isoform> isoforms,Eigen::MatrixXd &regionMatrix,ProgramOptions &opt){
    // create a index for the isoforms
    std::map<std::string, int> isoform_order;
    int isoform_index = 0;
    for (auto &[isoform_id,_]:isoforms){
        isoform_order[isoform_id] = isoform_index;
        isoform_index += 1;
    }
    int numRows = regions.size();
    int numCols = isoforms.size();
    Eigen::MatrixXd matrix = Eigen::MatrixXd::Constant(numRows, numCols, 0.0);
    
    for (int i = 0; i < numRows; ++i) {
        for (auto &isoform_id : regions[i].belong_isoform_id_set){
            if (isoform_order.count(isoform_id)!=0){
                if (opt.kvalue_entry_type == kvalue_entry_type_effective_length){
                    matrix(i, isoform_order[isoform_id]) = regions[i].eff_len;
                } else if (opt.kvalue_entry_type == kvalue_entry_type_binary){
                    matrix(i, isoform_order[isoform_id]) = 1;
                }
               
            }
        }
    }
    if (opt.not_normalize_entry){
        regionMatrix = matrix;
    } else{
        regionMatrix = normalizeCols(matrix);
    }
    

}
void getKvalue(Eigen::MatrixXd &regionMatrix,Gene & gene){
    Eigen::MatrixXd multiply_transpose_matrix = regionMatrix.transpose() * regionMatrix;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(multiply_transpose_matrix);
    Eigen::VectorXd singularValues = svd.singularValues();
    // double eps = 2.220446049250313e-16;
    double eps = 1.1920929e-07;
    double tolerance = eps * 
                       std::max(multiply_transpose_matrix.rows(), multiply_transpose_matrix.cols()) * singularValues.array().abs().maxCoeff();
    int rank = (singularValues.array() > tolerance).count();
    double svd_val_max = sqrt(singularValues(0));
    double svd_val_pos_min = sqrt(singularValues(rank-1));
    gene.kvalue = svd_val_max/svd_val_pos_min;
    gene.isFullRank = (rank == multiply_transpose_matrix.rows());
}
void calKvalue(std::vector<Gene> &genes,size_t start, size_t end,ProgramOptions &opt){
    FileWriter fileWriter;
    RegionGenerator regionGenerator;
    for (size_t i = start; i < end; ++i) {
        double kval;
        std::vector<Region> regions;
        regionGenerator.generateRegionForGene(genes[i],regions,opt);
        // fileWriter.writeRegion(opt.tempFolder+"/regions_"+genes[i].id+".jsonl",regions,opt);
        if (!genes[i].isoforms.empty() && !regions.empty()){
            Eigen::MatrixXd regionMatrix;
            createRegionMatrix(regions,genes[i].isoforms,regionMatrix,opt);
            if (opt.debug){
                fileWriter.writeRegion(opt.outputFolder+"/regions_"+genes[i].id+".jsonl",regions,opt);
                fileWriter.writeDesignMatrix(opt.outputFolder+"/design_matrix_"+genes[i].id+".txt",genes[i].id,genes[i].isoforms,regionMatrix,opt);
            }      
            getKvalue(regionMatrix,genes[i]);
        } else {
            genes[i].kvalue = kval;
        }
        
    }

}
void parallelCalKvalue(ProgramOptions &opt){
    program_log_message(Info,"Reading the annotation file...");
    std::map<std::string, Gene> gene_map;
    AnnotationParser annotationParser;
    annotationParser.decideAndParseFile(opt.annotationFile,gene_map);
    program_log_message(Info,"Done.");
    program_log_message(Info,"Calculating the K-value...");
    std::vector<Gene> genes;
    for (const auto& pair : gene_map) {
        genes.push_back(pair.second);
    };
    gene_map.clear();

    size_t numThreads = opt.threads;
    std::vector<std::thread> threads;
    size_t num_genes = genes.size();
    size_t genesPerThread = num_genes / numThreads;
    size_t additionalGenes = num_genes % numThreads;

    size_t start = 0;
    for (int i = 0; i < numThreads; ++i) {
        size_t end = start + genesPerThread + (i < additionalGenes ? 1 : 0);
        if (start < num_genes) {
            // calKvalue(genes,start,end,opt);
            threads.emplace_back(calKvalue, std::ref(genes), start, end,std::ref(opt));
        }
        start = end;
    }
    for (auto& t : threads) {
        t.join();
    }


    FileWriter fileWriter;
    fileWriter.writeKvalue(opt.outputFolder+"/kvalues.tsv",genes,opt);
    // if (opt.debug)
    fileWriter.writeIsoformKvalue(opt.outputFolder+"/gene_isoform_characteristics.tsv",genes,opt);
    program_log_message(Info,"Done.");
}