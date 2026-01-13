#ifndef CAL_KVALUE_H
#define CAL_KVALUE_H
#include <annotationParser.hh>
#include <common.hh>
#include <sstream>
void parallelCalKvalue(ProgramOptions &opt);
void calKvalue(std::vector<Gene> &genes,size_t start, size_t end,ProgramOptions &opt);

#endif