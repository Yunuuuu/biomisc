#include <Rcpp.h>
using namespace Rcpp;

List hyper_test_cpp(CharacterVector sample1, CharacterVector sample2, int n)
{
    List res(3);
    int count = intersect(sample1, sample2).size();
    res[0] = count;
    int m = sample1.size(), k = sample2.size();
    if (count <= m * k / n)
    // under-enrichment
    {
        res[1] = R::phyper(count, m, n - m, k, true, false);
        res[2] = -1;
    }
    else
    // over-enrichment
    {
        res[1] = R::phyper(count - 1, m, n - m, k, false, false);
        res[2] = 1;
    }

    return res;
}

// [[Rcpp::export]]
List calculate_hyper_overlap_cpp(CharacterVector sample1, CharacterVector sample2, int n, int stepsize)
{
    // the total length we need to calculate in each sample
    int list1_len = floor((sample1.size() - stepsize) / stepsize) + 1;
    int list2_len = floor((sample2.size() - stepsize) / stepsize) + 1;
    // Pre-allocate results matrix
    IntegerMatrix counts = no_init_matrix(list1_len, list2_len);
    NumericMatrix pvalue = no_init_matrix(list1_len, list2_len);
    IntegerMatrix signs = no_init_matrix(list1_len, list2_len);
    for (int i = 0; i < list1_len; i++)
    {
        CharacterVector list1 = sample1[Range(0, (i + 1) * stepsize - 1)];
        int m = list1.size();
        int non_m = n - m;
        for (int j = 0; j < list2_len; j++)
        {
            CharacterVector list2 = sample2[Range(0, (j + 1) * stepsize - 1)];
            int k = list2.size();
            int count = intersect(list1, list2).size();
            counts(i, j) = count;
            if (count <= m * k / n)
            // under-enrichment
            {
                pvalue(i, j) = R::phyper(count, m, non_m, k, true, false);
                signs(i, j) = -1;
            }
            else
            // over-enrichment
            {
                pvalue(i, j) = R::phyper(count - 1, m, non_m, k, false, false);
                signs(i, j) = 1;
            }
        }
    }
    return List::create(Named("counts") = counts,
                        Named("pvalue") = pvalue,
                        Named("signs") = signs);
}
