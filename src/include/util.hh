#ifndef UTIL_H
#define UTIL_H

#include <chrono>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <unordered_map>
bool isFile(const std::string fileName);
bool has_suffix(const std::string &str, const std::string &suffix);
bool isReferenceSequenceFile(const std::string &fileName);
bool isReadSequenceFile(const std::string &fileName);
int smallestPowerOfTwoLargerOrEqual(int n);

enum ProgramMessageLevel {
    Info = 0,
    Warning = 1,
    Error = 2,
    UnrecoverableError = 3};
void program_log_message(ProgramMessageLevel level, const char *message);
// static std::unordered_map<char, int> base64_index;
// static std::string base64_chars;
// Helper methods
size_t fromBase64(const std::string& base64); // Base64 decoding
std::string toBase64(size_t num); // Base64 encoding

// Utility method to create the Base64 index map; marked static to signify it's shared across instances
std::unordered_map<char, int> createBase64IndexMap();
void decompressLZ4File(const std::string &inputFile, std::istringstream &outputStream);
void emptyAndCreateFolder(const std::string directory_path);
void removeFolder(const std::string directory_path);

#endif