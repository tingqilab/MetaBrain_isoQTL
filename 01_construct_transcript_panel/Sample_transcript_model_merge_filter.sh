#!/bin/sh

## gffcompare merge
gffcompare -i ${gtf_list} -o sample_merged -T -D

### extract transcript supported by at least 2 independent samples
awk '
$3=="transcript" {
    if (match($0, /num_samples "([0-9]+)"/, a)) {
        ns = a[1];
        if (ns > 1) {
            match($0, /transcript_id "([^"]+)"/, t);
            keep[t[1]] = 1;
        }
    }
}
{
    if (match($0, /transcript_id "([^"]+)"/, t2)) {
        tid = t2[1];
        if (keep[tid] == 1) {
            print $0;
        }
    }
}' sample_merged.combined.gtf > sample_merged_two_samples_support.gtf