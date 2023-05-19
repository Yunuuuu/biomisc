#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
List rrho_hyper_overlap_cpp(CharacterVector sample1, CharacterVector sample2, int stepsize)
{
    int n = sample1.size();
    int n2 = sample2.size();

    // the total length we need to calculate in each sample
    int list1_len = ceil(n / stepsize);
    int list2_len = ceil(n2 / stepsize);

    // Pre-allocate results matrix
    IntegerMatrix counts_mat(list1_len, list2_len);
    NumericMatrix metrics(list1_len, list2_len);
    IntegerMatrix signs(list1_len, list2_len);

    int k = 0;
    for (int i = 0; i < list1_len; i++)
    {
        k += stepsize;
        CharacterVector list1 = sample1[Range(0, k - 1)];

        // inner loop
        int m = 0;
        int counts = 0;
        
        for (int j = 0; j < list2_len; j++)
        {
            m += stepsize;
            int total = list1.size();
            if (total > 0) 
            {
                list1 = setdiff(list1, sample2[Range(m - stepsize, m - 1)]);
                counts += total - list1.size();
            }

            if (counts > double(m) / n * k)
            // over-enrichment
            {
                signs(i, j) = 1;
                metrics(i, j) = R::phyper(counts - 1, m, n - m, k, false, true);
            }
            else
            // under-enrichment
            {
                signs(i, j) = -1;
                metrics(i, j) = R::phyper(counts, m, n - m, k, true, true);
            }
            counts_mat(i, j) = counts;
        }
    }
    return List::create(Named("counts") = counts_mat,
                        Named("metrics") = metrics,
                        Named("signs") = signs);
}
