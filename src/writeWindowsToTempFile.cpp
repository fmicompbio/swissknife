#include <iostream>
#include <fstream>
#include <string>
#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export(.writeWindowsToTempFileCPP)]]
CharacterVector writeWindowsToTempFileCPP(std::string chr, size_t w, std::string fname) {
    std::ofstream myfile;
    myfile.open(fname);
    if (myfile.is_open()) {
        for (size_t i = 0; i < chr.size() - w + 1; i++) {
            myfile << ">" << (i + 1) << std::endl << chr.substr(i, w) << std::endl;
        }
        myfile.close();
    } else {
        stop("Could not open file 'fname' for writing");
    }
    return fname;
}
