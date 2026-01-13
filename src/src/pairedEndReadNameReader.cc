// #include <pairedEndReadNameReader.hh>
// #include <iostream>
// #include <thread>
// #include <mutex>

// std::mutex mtx;

// PairedEndReadNameReader::PairedEndReadNameReader(const std::string& filename, int numThreads)
//     : filename(filename), numThreads(numThreads) {}

// PairedEndReadNameReader::~PairedEndReadNameReader() {}

// void PairedEndReadNameReader::readSequenceNames() {
//     gzFile gzfile = gzopen(filename.c_str(), "rb");
//     if (!gzfile) {
//         std::cerr << "Could not open file " << filename << std::endl;
//         return;
//     }

//     gzseek(gzfile, 0, SEEK_END);
//     long fileSize = gzseek(gzfile, 0, SEEK_END);
//     gzrewind(gzfile);
    
//     std::vector<std::thread> threads;
//     long chunkSize = fileSize / numThreads;

//     for (int i = 0; i < numThreads; ++i) {
//         long start = i * chunkSize;
//         long end = (i == numThreads - 1) ? fileSize : start + chunkSize;
//         threads.emplace_back(&PairedEndReadNameReader::processChunk, this, start, end);
//     }

//     for (auto &t : threads) {
//         t.join();
//     }

//     gzclose(gzfile);
// }

// void PairedEndReadNameReader::processChunk(long start, long end) {
//     gzFile gzfile = gzopen(filename.c_str(), "rb");
//     gzseek(gzfile, start, SEEK_SET);
//     char buffer[256];
//     long current = start;

//     while (current < end && gzgets(gzfile, buffer, sizeof(buffer)) != Z_NULL) {
//         std::string line(buffer);
//         if (line[0] == '@') {
//             std::lock_guard<std::mutex> lock(mtx);
//             size_t pos = line.find_first_of("\n\r");
//             sequenceNames.push_back(line.substr(0, pos));
//         }
//         // Skip the other 3 lines in the record
//         gzgets(gzfile, buffer, sizeof(buffer)); // Sequence
//         gzgets(gzfile, buffer, sizeof(buffer)); // '+'
//         gzgets(gzfile, buffer, sizeof(buffer)); // Quality
//         current = gztell(gzfile);
//     }

//     gzclose(gzfile);
// }

// const std::vector<std::string>& PairedEndReadNameReader::getSequenceNames() const {
//     return sequenceNames;
// }