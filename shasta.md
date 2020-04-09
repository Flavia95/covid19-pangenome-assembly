# Shasta parameters tuning - covid-19-bh20-assembly

**A de novo assembly of pulled down RNA sequenced on a nanopore device.***

A little temporary hack to remove uracils which shasta doesn't like in its default configuration.

<pre>zcat covid_update.fastq.gz | paste - - - - | tr ' ' '_' | tr -d '@' | tr 'U' 'T' | awk 'length($2) > 1500 { print ">"$1; print $2; }' > covid_update.1.5kb.UtoT.fasta</pre>

Modifying a little the default parameters

<pre>shasta-Linux-0.4.0 --input covid_update.1.5kb.UtoT.fasta --Reads.minReadLength 3460 --MarkerGraph.minCoverage 6 --MarkerGraph.maxCoverage 5000</pre>

we arrive to have a linear contig of nearly 20kbps, not enough:

![](images/01_change_min_len_bandage.png)

The first BLAST match is <a href='https://www.ncbi.nlm.nih.gov/nucleotide/MT007544.1?report=genbank&log$=nuclalign&blast_rank=1&RID=8XU4NDS5016'>Severe acute respiratory syndrome coronavirus 2 isolate Australia/VIC01/2020, complete genome, MT007544.1</a> (query coverge 100%, identity 97.91%).

Forcing the tools to work at higher coverage (modifying its coverage thresholds)

<pre>shasta-Linux-0.4.0 --input covid_update.1.5kb.UtoT.fasta --Reads.minReadLength 3460 --MarkerGraph.minCoverage 10 --MarkerGraph.maxCoverage 5000 --MinHash.maxBucketSize 100 --MarkerGraph.lowCoverageThreshold 20 --MarkerGraph.highCoverageThreshold 2560 --MarkerGraph.edgeMarkerSkipThreshold 1000</pre>

we get a bigger conting (31.5 kbps), but also a little mess.

![](images/02_change_coverage_parameters_bandage.png)

We still need to tune the pruning step and take conficence of the impact of this parameters tuning.
