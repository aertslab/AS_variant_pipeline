# Allele specific variant pipeline


## Step 1. Generation of personalized genomes and chain files

Here, we create a diploid personalized genome per samples by inserting phasedi
variants from whole genomes (single nucleotide variants, small insertions and
deletions, and structural variants) to reference genome (hg38).

Then we create chain files to go from personalized genome to reference genome,
and vice versa.


## Step2. Mapping functional genomics data to personalized genomes

ATAC-seq reads are mapped to HAP1, HAP2 and reference genome (hg38)
individually.

Bam file are further processed to remove duplicate reads and reads with low
mapping quality (Q2).

Peak calling is performed in each alignment file seperately with MACS2.

Peak summit locations from HAP1 and HAP2 are lifted over to hg38 coordinates,
filtered for blacklisted regions. Then peak summits are extended 250bp to each
direction, and peak scores within each peak set is normalized (by dividing
each peak score to the total peak score).

Then, three peak sets (HAP1, HAP2 and hg38) are consolidated: all peaks are
combined and ranked by their normalized score.
Starting with the highest scoring peak, we filter out any peak that overlaps
with the highest scoring peak. Then we move onto the second peak, and repeat
this process until we have a non-overlapping peak set. This step ensures we
obtain an accurate peak set considering all three alignments.


## Step 3. Obtaining allele specific counts and predict ASCAVs

At this step, we aim to obtain allele specific counts.

We first start by comparing HAP1 and HAP2 bam files generated in the first
step. We evaluate each read in these bam files in order to identify its most likely
origin. We do this by looking first at mapping quality, then at the number of
mismatches. If a read has higher mapping quality for either of the haplotypes
it is assigned to that haplotype. If the mapping quality is the same, we look
at the number of mismatches and assigned the read to the haplotype for which
the read is aligned with least number of mismatches. If both metrics are
equivalent, we deem the read as commonly mapping to both haplotypes. This step
results in four bam files: `HAP1.unique.bam`, `HAP2.unique.bam`,
`HAP1.common.bam` and `HAP2.common.bam` (`HAP1.common.bam` and
`HAP2.common.bam` should essentially be identical).

We next check duplicate and ambiguously mapping reads. We mark the duplicate
reads in each of these bam files with Picard. Then, by lifting over the read
positions in `HAP1.common.bam` and `HAP2.common.bam` to reference genome, we
check for ambiguously mapping reads (if the reads are truly indistinguishable
between two haplotypes, they should correspond to the same location when
lifted over to the reference genome).

Next, we obtain phased heterozygous variants from whole genome sequence data
and overlap them with consolidated peaks obtained in the previous step. These
variants will be the ones that we will try to obtain allele specific counts.
SNV locations in peaks are lifted over to HAP1 and HAP2, and we count the
alleles at these variants using HAP1.bam and HAP2.bam files with
`samtools mpileup`. We also count the alleles using the `COMMON.bam` file.

Next, we combine allelic counts from HAP1, HAP2 and COMMON alignments, add
genomic allelic ratios from WGS data, and predict allele specific chromatin
accessibility variants using BaalCHIP.


## Step 4. Scoring ASB and control events with Cluster-Buster

At this step, we evaluate the impact of variants on transcription factor motif
binding by evaluating cluster-buster scores of ASCAVs over non-ASCAVs. We
perform this analysis in peak-wise and not variant-wise. Thus, we start by
filtering peaks that contain discordant ASCAVs (ie. the peak has multiple
ASCAVs with different winning haplotypes). Then the remaining peaks are lifted
over to HAP1 or HAP2 (depending on the winning haplotype) and sequences were
extracted from corresponding haplotype fastas. Next, we select control peaks
that are not predicted as ASCAVs from the BaalCHIP output, again filtering for
discordant variants. Peak positions are lifted over to HAP1 or HAP2, and
sequences were extracted as in previous step. For each peak, the sequences
from the reference genome were extracted as well. And finally, all sequences
are scored with cluster-buster using a motif collection of 22k motifs.

