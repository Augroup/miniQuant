// #ifndef PAIRED_END_READ_NAME_READER_H
// #define PAIRED_END_READ_NAME_READER_H

// #include <string>
// #include <vector>
// #include <zlib.h>

// class PairedEndReadNameReader {
// public:
//     PairedEndReadNameReader(const std::string& filename, int numThreads);
//     ~PairedEndReadNameReader();

//     void readSequenceNames();
//     const std::vector<std::string>& getSequenceNames() const;

// private:
//     std::string filename;
//     int numThreads;
//     std::vector<std::string> sequenceNames;

//     void processChunk(long start, long end);
// };

// #endif // PAIRED_END_READ_NAME_READER_H