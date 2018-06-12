#/!bin/bash

#haplotype1_input_bam_filename=/staging/leuven/stg_00002/lcb/zkalender/melanoma_bulk_AC/20.bam_files/20.personalized_reference_mapping/A375_hap1.sam
#haplotype2_input_bam_filename=/staging/leuven/stg_00002/lcb/zkalender/melanoma_bulk_AC/20.bam_files/20.personalized_reference_mapping/A375_hap2.sam

#samtools sort -n -O sam ${haplotype1_input_bam_filename} > /staging/leuven/stg_00002/lcb/ghuls/tmp/A375_hap1.sam
#samtools sort -n -O sam ${haplotype2_input_bam_filename} > /staging/leuven/stg_00002/lcb/ghuls/tmp/A375_hap2.sam

haplotype1_input_bam_filename=/staging/leuven/stg_00002/lcb/ghuls/tmp/A375_hap1.sam
haplotype2_input_bam_filename=/staging/leuven/stg_00002/lcb/ghuls/tmp/A375_hap2.sam

haplotype1_output_bam_filename=/staging/leuven/stg_00002/lcb/ghuls/tmp/A375_hap1.phased.sam
haplotype2_output_bam_filename=/staging/leuven/stg_00002/lcb/ghuls/tmp/A375_hap2.phased.sam

mawk \
        -v haplotype1_input_bam_filename="${haplotype1_input_bam_filename}" \
        -v haplotype2_input_bam_filename="${haplotype2_input_bam_filename}" \
        -v haplotype1_output_bam_filename="${haplotype1_output_bam_filename}" \
        -v haplotype2_output_bam_filename="${haplotype2_output_bam_filename}" \
        '
        BEGIN {
                OFS=FS="\t";

                # Print BAM header of haplotype 1 input BAM file to the haplotype 1 output BAM file.
                while ( ( (getline < haplotype1_input_bam_filename) > 0 ) && ( $0 ~ /^@/ ) ) {
                    haplotype1_input_bam_line_number += 1;

                    print $0 > haplotype1_output_bam_filename;
                }

                # Wrtie BAM header of haplotype 2 input BAM file to the haplotype 2 output BAM file.
                while ( ( (getline < haplotype2_input_bam_filename) > 0 ) && ( $0 ~ /^@/ ) ) {
                    haplotype2_input_bam_line_number += 1;

                    print $0 > haplotype2_output_bam_filename;
                }

                haplotype1_input_bam_eof = 0
                haplotype2_input_bam_eof = 0

                while (haplotype1_input_bam_eof == 0 && haplotype1_input_bam_eof == 0) {
                    if ( (getline < haplotype1_input_bam_filename) > 0 ) {
                        haplotype1_input_bam_line_number += 1;

                        haplotype1["qname"] = $1;
                        haplotype1["rname"] = $3;
                        haplotype1["pos"] = $4;
                        haplotype1["mapq"] = $5;
                        haplotype1["cigar"] = $6;
                        haplotype1["line"] = $0;
                    } else {
                        haplotype1_input_bam_eof = 1;
                        break;
                    }

                    if ( (getline < haplotype2_input_bam_filename) > 0 ) {
                        haplotype2_input_bam_line_number += 1;

                        haplotype2["qname"] = $1;
                        haplotype2["rname"] = $3;
                        haplotype2["pos"] = $4;
                        haplotype2["mapq"] = $5;
                        haplotype2["cigar"] = $6;
                        haplotype2["line"] = $0;
                    } else {
                        haplotype2_input_bam_eof = 1;
                        break;
                    }

                    if (haplotype1["qname"] == haplotype2["qname"]) {
                        if ( (haplotype1["mapq"] >= 4) && (haplotype2["mapq"] >= 4)) {
                            # Read maps uniquely both in haplotype 1 and 2.

                            if ( haplotype1["rname"] == haplotype2["rname"] ) {
                                # Read maps to the same chromosome.
#                                if ( haplotype1["pos"] == haplotype2["pos"] ) {
                                    # Read maps to the same position.
                                    if ( haplotype1["cigar"] == haplotype2["cigar"] ) {
                                        # Read has the same cigar string.

                                        # Write this read to the haplotype 1 and 2 output BAM file.
                                        print haplotype1["line"] > haplotype1_output_bam_filename;
                                        print haplotype2["line"] > haplotype2_output_bam_filename;
                                    } else {
                                        printf("Error: CIGAR string of read '%s' is different for haplotype 1 and 2 BAM files:", haplotype1["qname"]);
                                        printf("  - haplotype 1 (line number: %d): %s\n", haplotype1_input_bam_line_number, haplotype1["line"]);
                                        printf("  - haplotype 2 (line number: %d): %s\n", haplotype2_input_bam_line_number, haplotype2["line"]);
                                    }
#                                } else {
#                                    printf("Error: Position of read '%s' is different for haplotype 1 and 2 BAM files:", haplotype1["qname"]);
#                                    printf("  - haplotype 1 (line number: %d): %s\n", haplotype1_input_bam_line_number, haplotype1["line"]);
#                                    printf("  - haplotype 2 (line number: %d): %s\n", haplotype2_input_bam_line_number, haplotype2["line"]);
#                                }
                            } else {
                                    printf("Error: Chromosome of read '%s' is different for haplotype 1 and 2 BAM files:", haplotype1["qname"]);
                                    printf("  - haplotype 1 (line number: %d): %s\n", haplotype1_input_bam_line_number, haplotype1["line"]);
                                    printf("  - haplotype 2 (line number: %d): %s\n", haplotype2_input_bam_line_number, haplotype2["line"]);
                            }
                        } else if (haplotype1["mapq"] >= 4) {
                            # Read maps only uniquely in haplotype 1.

                            # Write this read to the haplotype 1 output BAM file.
                            print haplotype1["line"] > haplotype1_output_bam_filename;
                        } else if (haplotype2["mapq"] >= 4) {
                            # Read maps only uniquely in haplotype 2.

                            # Write this read to the haplotype 2 output BAM file.
                            print haplotype2["line"] > haplotype2_output_bam_filename;
                        }
                    } else {
                        print "Error: haplotype 1 and 2 BAM files are not sorted by read name:";
                        printf("  - haplotype 1 (line number: %d): %s\n", haplotype1_input_bam_line_number, haplotype1["line"]);
                        printf("  - haplotype 2 (line number: %d): %s\n", haplotype2_input_bam_line_number, haplotype2["line"]);
                    }
            }
        }
'
