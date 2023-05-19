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
    IntegerMatrix counts(list1_len, list2_len);
    NumericMatrix metrics(list1_len, list2_len);
    IntegerMatrix signs(list1_len, list2_len);
    for (int i = 0; i < list1_len; i++)
    {
        int k = (i + 1) * stepsize;
        CharacterVector list1 = sample1[Range(0, k - 1)];

        for (int j = 0; j < list2_len; j++)
        {
            int m = (j + 1) * stepsize;
            CharacterVector list2 = sample2[Range(0, m - 1)];
            int count = intersect(list1, list2).size();
            counts(i, j) = count;
            if (count > double(m) / n * k)
            // over-enrichment
            {
                signs(i, j) = 1;
                metrics(i, j) = R::phyper(count - 1, m, n - m, k, false, true);
            }
            else
            // under-enrichment
            {
                signs(i, j) = -1;
                metrics(i, j) = R::phyper(count, m, n - m, k, true, true);
            }
        }
    }
    return List::create(Named("counts") = counts,
                        Named("metrics") = metrics,
                        Named("signs") = signs);
}
