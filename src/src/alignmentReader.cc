#include <iostream>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>
#include <alignmentReader.hh>
#include <EM.hh>
#include <util.hh>

// Method to read and process JSON into Read and Hit directly
void AlignmentReader::processLongReadsAlignments(std::vector<Transcript>& transcripts,
                                        std::vector<Read>& reads,
                                        std::vector<size_t>& thread_bucket_id,
                                        int chunk_id,
                                        ProgramOptions& opt) {
    size_t num_long_reads = 0;
    for (const auto & bucket_id : thread_bucket_id){
        std::istringstream decompressedStream;
        decompressLZ4File(opt.tempFolder+"/long_reads_alignments/bucket_"+std::to_string(chunk_id)+'_'+std::to_string(bucket_id)+".jsonl.lz4", decompressedStream);
        std::string line;
        while (std::getline(decompressedStream, line)) {
            try {
                auto jsonObject = nlohmann::json::parse(line);
                
                size_t read_id = fromBase64(jsonObject.at("i").get<std::string>());
                size_t num_kmers = jsonObject.at("t").get<size_t>();
                size_t num_mapped_kmers = jsonObject.at("c").get<size_t>();
                std::vector<Hit> hits;

                if (jsonObject.contains("m") && jsonObject["m"].is_object()) {
                    for (const auto& [key, _] : jsonObject["m"].items()) {
                        size_t transcript_id = static_cast<size_t>(fromBase64(key));
                        if (transcript_id < transcripts.size()) {
                            Hit hit = Hit(transcript_id, 1.0);
                            if (opt.debug){
                                hit.set_mapped_kmer(jsonObject["m"][toBase64(transcript_id)]);
                            }
                            hits.push_back(hit);
                        }
                    }
                }

                // Create Read with its hits
                if (hits.size() > 0){
                    Read read = Read(num_kmers + opt.long_reads_kmer_size - 1, hits,LongReadType);
                    if (opt.debug){
                        read.set_mapped_kmer(num_mapped_kmers, num_kmers);
                        read.set_id(read_id);
                    }
                    reads.push_back(read);
                    num_long_reads += 1;
                }
            } catch (const nlohmann::json::exception& e) {
                std:: cout << line <<std::endl;
                std::cout << bucket_id<< std::endl;
                std::cerr << "JSON error: " << e.what() << std::endl;
            } catch (const std::invalid_argument& ia) {
                std::cerr << "Base64 conversion error: " << ia.what() << std::endl;
            }
        }

    }

}

// Method to read and process JSON into Read and Hit directly
void AlignmentReader::processPairedShortReadsAlignments(
    std::vector<Transcript>& transcripts,
    std::vector<Read>& reads,
    std::vector<size_t>& thread_bucket_id,
    ProgramOptions& opt) {

    size_t num_short_reads =  0;
    for (const auto & bucket_id : thread_bucket_id){
        bool eof1 = false, eof2 = false;
        std::istringstream decompressedStream_R1;
        std::istringstream decompressedStream_R2;
        decompressLZ4File(opt.tempFolder+"/short_reads_alignments_R1/bucket_"+std::to_string(bucket_id)+".jsonl.lz4", decompressedStream_R1);
        decompressLZ4File(opt.tempFolder+"/short_reads_alignments_R2/bucket_"+std::to_string(bucket_id)+".jsonl.lz4", decompressedStream_R2);
        std::string line1;
        std::string line2;
        // Read the first line from both files to initialize
        if (!std::getline(decompressedStream_R1, line1)) eof1 = true;
        if (!std::getline(decompressedStream_R2, line2)) eof2 = true;

        while (!eof1 && !eof2) {
            try {
                auto jsonObject1 = nlohmann::json::parse(line1);
                auto jsonObject2 = nlohmann::json::parse(line2);

                size_t read_id1 = fromBase64(jsonObject1.at("i").get<std::string>());
                size_t read_id2 = fromBase64(jsonObject2.at("i").get<std::string>());

                if (read_id1 == read_id2) {
                    // Process and merge both reads as they have the same read_id

                    size_t num_kmers1 = jsonObject1.at("t").get<size_t>();
                    size_t num_kmers2 = jsonObject2.at("t").get<size_t>();
                    size_t num_mapped_kmers1 = jsonObject1.at("c").get<size_t>();
                    size_t num_mapped_kmers2 = jsonObject2.at("c").get<size_t>();

                    std::vector<size_t> transcript_id_set1;
                    std::vector<size_t> transcript_id_set2;
                    
                    if (jsonObject1.contains("m") && jsonObject1["m"].is_object()) {
                        for (const auto& [key, _] : jsonObject1["m"].items()) {
                            size_t transcript_id = static_cast<size_t>(fromBase64(key));
                            if (transcript_id < transcripts.size()) {
                                transcript_id_set1.push_back(transcript_id);
                            }
                        }
                    }

                    if (jsonObject2.contains("m") && jsonObject2["m"].is_object()) {
                        for (const auto& [key, _] : jsonObject2["m"].items()) {
                            size_t transcript_id = static_cast<size_t>(fromBase64(key));
                            if (transcript_id < transcripts.size()) {
                                transcript_id_set2.push_back(transcript_id);
                            }
                        }
                    }
                    sort(transcript_id_set1.begin(), transcript_id_set1.end());
                    sort(transcript_id_set2.begin(), transcript_id_set2.end());
                    std::vector<size_t> intersected_transcripts;

                    std::set_intersection(transcript_id_set1.begin(), transcript_id_set1.end(),
                                        transcript_id_set2.begin(), transcript_id_set2.end(),
                                        std::back_inserter(intersected_transcripts));
                    std::vector<Hit> hits;
                    hits.reserve(intersected_transcripts.size());
                    for (const auto & transcript_id : intersected_transcripts){
                        Hit hit = Hit(transcript_id, 1.0);
                        if (opt.debug){
                            hit.set_mapped_kmer(jsonObject1["m"][toBase64(transcript_id)], jsonObject2["m"][toBase64(transcript_id)]);
                        }
                        hits.push_back(hit);
                    }


                    // Create a Read with its combined hits
                    if (hits.size() > 0){
                        Read read = Read(num_kmers1 + num_kmers2 + 2 * opt.short_reads_kmer_size - 2, hits,ShortReadType);
                        if (opt.debug){
                            read.set_mapped_kmer(num_mapped_kmers1, num_kmers1,num_mapped_kmers2, num_kmers2);
                            read.set_id(read_id1);
                        }
                        reads.push_back(read);
                        num_short_reads += 1;
                    }
                    

                    // Move to the next line in both files
                    if (!std::getline(decompressedStream_R1, line1)) eof1 = true;
                    if (!std::getline(decompressedStream_R2, line2)) eof2 = true;

                } else if (read_id1 < read_id2) {
                    // Move to the next line only in file_pair1
                    if (!std::getline(decompressedStream_R1, line1)) eof1 = true;
                } else {
                    // read_id2 < read_id1, move to the next line only in file_pair2
                    if (!std::getline(decompressedStream_R2, line2)) eof2 = true;
                }

            } catch (const nlohmann::json::exception& e) {
                std::cerr << "JSON error: " << e.what() << std::endl;
                // Advance the files in case of parsing error to avoid an infinite loop
                if (!std::getline(decompressedStream_R1, line1)) eof1 = true;
                if (!std::getline(decompressedStream_R2, line2)) eof2 = true;
            } catch (const std::invalid_argument& ia) {
                std::cerr << "Base64 conversion error: " << ia.what() << std::endl;
                // Advance the files in case of conversion error to avoid an infinite loop
                if (!std::getline(decompressedStream_R1, line1)) eof1 = true;
                if (!std::getline(decompressedStream_R2, line2)) eof2 = true;
            }
        }    
    }
}
void AlignmentReader::parallelProcessAlignments(
    std::vector<Transcript>& transcripts,
    std::vector<Read>& reads,
    ProgramOptions& opt,ReadType readType){

    size_t num_buckets = 1 << 10;
    size_t bucketPerThread = num_buckets / opt.threads;
    size_t additionalBuckets = num_buckets % opt.threads;

    if (readType == ShortReadType){
        std::vector<std::thread> threads;
        std::vector<std::vector<Read>> all_reads(opt.threads);
        std::vector<std::vector<size_t>> all_bucket_id(opt.threads);
    
        size_t start = 0;
        for (int i = 0; i < opt.threads; ++i) {
            size_t end = start + bucketPerThread + (i < additionalBuckets ? 1 : 0);
            if (start < num_buckets) {
                for (size_t id = start; id < end; ++id) {
                    all_bucket_id[i].push_back(id);
                }
                threads.emplace_back(
                    &AlignmentReader::processPairedShortReadsAlignments, // Member function pointer
                    this,                                  // Instance of the class
                    std::ref(transcripts),                   // Pass by reference using std::ref
                    std::ref(all_reads[i]),
                    std::ref(all_bucket_id[i]),
                    std::ref(opt)
                );
            }
            start = end;
        }
        for (auto& t : threads) {
            t.join();
        }
        for (auto& thread_reads : all_reads) {
            reads.insert(end(reads), begin(thread_reads), end(thread_reads));
        }
        
        opt.num_short_reads = reads.size();
        program_log_message(Info, ("Mapped "+std::to_string(opt.num_short_reads)+ " properly paried short read pairs.").c_str());
    }else if (readType == LongReadType){
        for (int chunk_id = 0; chunk_id < opt.num_long_reads_chunks; chunk_id++){
            std::vector<std::thread> threads;
            std::vector<std::vector<Read>> all_reads(opt.threads);
            std::vector<std::vector<size_t>> all_bucket_id(opt.threads);
        
            size_t start = 0;
            for (int i = 0; i < opt.threads; ++i) {
                size_t end = start + bucketPerThread + (i < additionalBuckets ? 1 : 0);
                if (start < num_buckets) {
                    for (size_t id = start; id < end; ++id) {
                        all_bucket_id[i].push_back(id);
                    }

                    threads.emplace_back(
                        &AlignmentReader::processLongReadsAlignments, // Member function pointer
                        this,                                  // Instance of the class
                        std::ref(transcripts),                   // Pass by reference using std::ref
                        std::ref(all_reads[i]),
                        std::ref(all_bucket_id[i]),
                        chunk_id,
                        std::ref(opt)
                    );
                }
                start = end;
            }
            for (auto& t : threads) {
                t.join();
            }
            for (auto& thread_reads : all_reads) {
                reads.insert(end(reads), begin(thread_reads), end(thread_reads));
            }
    
        }
        opt.num_long_reads = reads.size();
        program_log_message(Info, ("Mapped "+std::to_string(opt.num_long_reads)+ " long reads.").c_str());

    }
    opt.long_reads_weight_sum = (double) opt.num_long_reads;
    opt.short_reads_weight_sum = (double) opt.num_short_reads;
 }