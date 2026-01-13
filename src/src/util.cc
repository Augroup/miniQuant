#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <util.hh>
#include <unordered_map>
#include <lz4frame.h>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdexcept>
#include <filesystem.hpp>
namespace fs = ghc::filesystem;

bool isFile(const std::string fileName){
    std::ifstream infile(fileName);
    return infile.good();
}
bool has_suffix(const std::string &str, const std::string &suffix){
    return str.size() >= suffix.size() &&
           str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}
bool isReferenceSequenceFile(const std::string &fileName){
    return has_suffix(fileName,".fasta") || has_suffix(fileName,".fa");
}
bool isReadSequenceFile(const std::string &fileName){
    return has_suffix(fileName,".fasta") || has_suffix(fileName,".fa") || \
    has_suffix(fileName,".fasta.gz") || has_suffix(fileName,".fa.gz") || \
    has_suffix(fileName,".fastq") || has_suffix(fileName,".fq") || \
    has_suffix(fileName,".fastq.gz") || has_suffix(fileName,".fq.gz");
}
int smallestPowerOfTwoLargerOrEqual(int n) {
    if (n < 1) return 0;  // No power of two is less than or equal to zero

    int result = 1;

    // Shift left until the result is at least n
    while (result < n) {
        result <<= 1;  // Equivalent to multiplying result by 2
    }

    return result;
}

// enum ProgramMessageLevel {
//     Info = 0,
//     Warning = 1,
//     Error = 2,
//     UnrecoverableError = 3};
void program_log_message(ProgramMessageLevel level, const char *message){
    // Get the current time as a time_point
    auto now = std::chrono::system_clock::now();
    
    // Convert the time_point to a time_t, which represents time in seconds since epoch
    std::time_t currentTime = std::chrono::system_clock::to_time_t(now);
    
    // Convert to local time and format as string
    std::tm localTime = *std::localtime(&currentTime);
    switch (level)
    {
    case Info:
    {
        std::cout << "[INFO][" << std::put_time(&localTime, "%Y-%m-%d %H:%M:%S") << "] " << message << "\n";
    }
    break;
    case Warning:
    {
    std::cout << "[WARNING] " << message <<"\n";
    }
    break;
    case Error:
    {
        std::cout << "[ERROR] " << message <<"\n";
    }
    break;
    case UnrecoverableError:
    {
        std::cout << "[FATAL_ERROR] " << message <<"\n";
        exit(1);
    }
    };
}
static std::string base64_chars = 
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    "abcdefghijklmnopqrstuvwxyz"
    "0123456789+/";
std::unordered_map<char, int>createBase64IndexMap(){
    std::unordered_map<char, int> indexMap;
    for (int i = 0; i < base64_chars.size(); ++i) {
        indexMap[base64_chars[i]] = i;
    }

    return indexMap;
};


// Initialize the map only once
static std::unordered_map<char, int> base64_index = createBase64IndexMap();

size_t fromBase64(const std::string& base64){
    size_t num = 0;
    size_t power = 1; // Represents 64^i

    // Iterate from the end to the beginning to calculate the base 10 number
    for (auto it = base64.rbegin(); it != base64.rend(); ++it) {
        char c = *it;
        
        auto map_it = base64_index.find(c);
        if (map_it == base64_index.end()) {
            std::ostringstream error_msg;
            error_msg << "Intermediate files in output folder (-o) is corrupted. Make sure no other program use it. Exit now!";
            program_log_message(UnrecoverableError,error_msg.str().c_str());
            // throw std::invalid_argument("Invalid Base64 character found.");
        }
        
        num += map_it->second * power;
        power *= 64;
    }

    return num;
};

std::string toBase64(size_t num) {
    std::string base64;
    
    if (num == 0) {
        return "A"; // Base64 representation of zero
    }

    while (num > 0) {
        size_t remainder = num % 64;
        base64 += base64_chars[remainder];
        num /= 64;
    }

    // Base64 digits are reverse order compared to typical numerical systems
    std::reverse(base64.begin(), base64.end());
    
    return base64;
}

void decompressLZ4File(const std::string &inputFile, std::istringstream &outputStream) {
    std::ifstream input(inputFile, std::ios::binary);
    if (!input) {
        std::ostringstream error_msg;
        error_msg << "Intermediate files in output folder (-o) is corrupted. Make sure no other program use it. Exit now!";
        program_log_message(UnrecoverableError,error_msg.str().c_str());
    }

    std::vector<char> compressedData((std::istreambuf_iterator<char>(input)), std::istreambuf_iterator<char>());
    input.close();

    // Create a decompression context
    LZ4F_dctx *dctx;
    size_t decompressionOk = LZ4F_createDecompressionContext(&dctx, LZ4F_VERSION);
    if (LZ4F_isError(decompressionOk)) {
        std::ostringstream error_msg;
        error_msg << "Intermediate files in output folder (-o) is corrupted. Make sure no other program use it. Exit now!";
        program_log_message(UnrecoverableError,error_msg.str().c_str());
        // throw std::runtime_error("Failed to create decompression context.");
    }

    // A buffer for decompressed data
    std::vector<char> buffer(4096); // Adjust buffer size as needed
    std::stringstream decompressedContent;

    size_t offset = 0;
    size_t remaining = compressedData.size();
    const char *src = compressedData.data();

    while (remaining > 0) {
        size_t decompressedSize = buffer.size();
        size_t readSize = remaining;

        decompressionOk = LZ4F_decompress(dctx, buffer.data(), &decompressedSize, src + offset, &readSize, nullptr);
        if (LZ4F_isError(decompressionOk)) {
            LZ4F_freeDecompressionContext(dctx);
            std::ostringstream error_msg;
            error_msg << "Intermediate files in output folder (-o) is corrupted. Make sure no other program use it. Exit now!";
            program_log_message(UnrecoverableError,error_msg.str().c_str());
        }

        decompressedContent.write(buffer.data(), decompressedSize);

        offset += readSize;
        remaining -= readSize;
    }

    LZ4F_freeDecompressionContext(dctx);

    // Pass decompressed data to the output stream
    outputStream.str(decompressedContent.str());
}

void emptyAndCreateFolder(const std::string directory_path){
    try {
        // Check if the directory exists
        if (fs::exists(directory_path)) {
            // Remove the directory and its contents
            fs::remove_all(directory_path);
        }

        // Create the directory along with any necessary parent directories
        fs::create_directories(directory_path); 
    } catch (const std::exception& e) {
        std::ostringstream error_msg;
        error_msg << "Failed to create folder of " <<directory_path << ". Make sure you have permission to write. Exit now!";
        program_log_message(UnrecoverableError,error_msg.str().c_str());
    }
}
void removeFolder(const std::string directory_path){
    try {
        // Check if the directory exists
        if (fs::exists(directory_path)) {
            // Remove the directory and its contents
            fs::remove_all(directory_path);
        }
    } catch (const std::exception& e) {
        std::ostringstream error_msg;
        error_msg << "Failed to delete folder of " <<directory_path << ". Make sure you have permission to write. Exit now!";
        program_log_message(UnrecoverableError,error_msg.str().c_str());
    }
}

