#include <Rcpp.h>
using namespace Rcpp;

//' @title [internal] Count frequency of differences between values in integer vectors.
//'
//' @description Given two ascendingly sorted integer vectors \code{x} and \code{p}, calculate
//'     and cound the differences between their elements that are greater than zero
//'     and less than \code{maxd}. The number of observed distances \code{d} are reported in \code{cnt[d]},
//'     and `maxd` corresponds to the \code{length(cnt)}.
//'
//' @author Michael Stadler
//'
//' @param x first \code{integer} vector.
//' @param p second \code{integer} vector. Distances are calculated from each element
//'     in \code{p} to each greater element in \code{x}.
//' @param cnt \code{NumericVector} to store the result in. The length of \code{cnt} defines
//'     the maximal distance that will be included in the analysis, and new counts
//'     will be added to the values of \code{cnt}.
//'
//' @return \code{numeric} vector \code{cnt}, where \code{cnt[d]} correspond to
//'     the number of observed distances \code{d}.
//'
// [[Rcpp::export]]
NumericVector calcAndCountDist(std::vector<int> x, std::vector<int> p, NumericVector cnt) {
    int pi, xi, pos, d, xiold = 0, maxd = cnt.size();
    for (pi=0; pi < p.size(); pi++) {
        pos = p[pi];
        // fast forward xi until x[xi] > pos
        while (xiold < x.size() && pos > x[xiold]) {
            xiold++;
        }
        for (xi = xiold; xi < x.size(); xi++) {
            d = x[xi] - pos + 1;
            if (d <= maxd) {
                cnt[d-1]++;
            } else {
                break;
            }
        }
    }
    return cnt;
}
