#!/bin/bash
#
# Copyright (C): 2018-2019 - Gert Hulselmans
#
# Purpose: Create haplotype specific BAM files where each haplotype BAM file only
#          contains those reads that map better to that haplotype than to the other
#          haplotype or contains those reads that map exactly the same (same mapQ,
#          CIGAR and XM tag value) to both haplotypes.
#          To create the input haplotype BAM files, the reads for a sample should be
#          mapped with bowtie2 to haplotype 1 and haplotype 2 reference separately.



create_haplotype_specific_bam_files_from_bowtie2_haplotype_mapped_bam_files () {
    if [ ${#@} -ne 5 ] ; then
        printf 'Usage:  create_haplotype_specific_bam_files_from_bowtie2_haplotype_mapped_bam_files \\\n';
        printf '            haplotype1_input_bam_filename \\\n';
        printf '            haplotype2_input_bam_filename \\\n';
        printf '            haplotype1_output_bam_filename_prefix \\\n';
        printf '            haplotype2_output_bam_filename_prefix \\\n';
        printf '            pos_sorted|rn_sorted\n\n';
        printf 'Arguments:\n\n';
        printf '    * haplotype1_input_bam_filename:\n';
        printf '        input BAM file mapped with bowtie2 against haplotype 1 reference.\n';
        printf '    * haplotype2_input_bam_filename:\n';
        printf '        input BAM file mapped with bowtie2 against haplotype 2 reference.\n';
        printf '    * haplotype1_output_bam_filename_prefix:\n';
        printf '        prefix for output BAM file names for haplotype 1:\n';
        printf '          - ${haplotype1_output_bam_filename_prefix}.all.bam:\n';
        printf '            all reads that map to haplotype 1.\n';
        printf '          - ${haplotype1_output_bam_filename_prefix}.common.bam:\n';
        printf '            reads which map equally well to haplotype 1 as to haplotype 2.\n';
        printf '          - ${haplotype1_output_bam_filename_prefix}.unique.bam:\n';
        printf '            reads which map the best to haplotype 1.\n';
        printf '    * haplotype2_output_bam_filename_prefix:\n';
        printf '        prefix for output BAM file names for haplotype 2:\n';
        printf '          - ${haplotype2_output_bam_filename_prefix}.all.bam:\n';
        printf '            all reads that map to haplotype 2.\n';
        printf '          - ${haplotype2_output_bam_filename_prefix}.common.bam:\n';
        printf '            reads which map equally well to haplotype 1 as to haplotype 2.\n';
        printf '          - ${haplotype2_output_bam_filename_prefix}.unique.bam:\n';
        printf '            reads which map the best to haplotype 2.\n';
        printf '    * pos_sorted|rn_sorted:\n';
        printf '        Sort haplotype output BAM files by chromosomal position or by read name.\n\n';
        printf 'Purpose:\n\n';
        printf '    Create haplotype specific BAM files where each haplotype BAM file only\n';
        printf '    contains those reads that map better to that haplotype than to the other\n';
        printf '    haplotype or contains those reads that map exactly the same (same mapQ,\n';
        printf '    CIGAR and XM tag value) to both haplotypes.\n';
        printf '    To create the input haplotype BAM files, the reads for a sample should be\n';
        printf '    mapped with bowtie2 to haplotype 1 and haplotype 2 reference separately.\n\n';
        return 1;
    fi

    local haplotype1_input_bam_filename="${1}";
    local haplotype2_input_bam_filename="${2}";

    local haplotype1_output_bam_filename_prefix="${3%.bam}";
    local haplotype2_output_bam_filename_prefix="${4%.bam}";

    local sorted_type="${5}";
    local samtools_sort_order_argument="";

    case "${sorted_type}" in
        pos_sorted)
            # Sorted by genomic coordinate.
            samtools_sort_order_argument="";
            ;;
        rn_sorted)
            # Sorted by read name.
            samtools_sort_order_argument=" -n ";
            ;;
        *)
            printf 'Error: Argument 5 "%s" is not "pos_sorted" (position sorted) or "rn_sorted" (read sorted).\n' "${sorted_type}";
            return 1;
            ;;
    esac

    local haplotype1_all_output_bam_filename="${haplotype1_output_bam_filename_prefix}.all.bam";
    local haplotype2_all_output_bam_filename="${haplotype2_output_bam_filename_prefix}.all.bam";
    local haplotype1_common_output_bam_filename="${haplotype1_output_bam_filename_prefix}.common.bam";
    local haplotype2_common_output_bam_filename="${haplotype2_output_bam_filename_prefix}.common.bam";
    local haplotype1_unique_output_bam_filename="${haplotype1_output_bam_filename_prefix}.unique.bam";
    local haplotype2_unique_output_bam_filename="${haplotype2_output_bam_filename_prefix}.unique.bam";

    local samtools_threads=4;

    if [ ! -f "${haplotype1_input_bam_filename}" ] ; then
        printf 'Error: Haplotype 1 BAM file "%s" does not exist.\n.' "${haplotype1_input_bam_filename}";
        return 1;
    fi

    if [ ! -f "${haplotype2_input_bam_filename}" ] ; then
        printf 'Error: Haplotype 2 BAM file "%s" does not exist.\n.' "${haplotype2_input_bam_filename}";
        return 1;
    fi


    mawk \
            -v haplotype1_input_bam_filename="${haplotype1_input_bam_filename}" \
            -v haplotype2_input_bam_filename="${haplotype2_input_bam_filename}" \
            -v haplotype1_all_output_bam_filename="${haplotype1_all_output_bam_filename}" \
            -v haplotype2_all_output_bam_filename="${haplotype2_all_output_bam_filename}" \
            -v haplotype1_common_output_bam_filename="${haplotype1_common_output_bam_filename}" \
            -v haplotype2_common_output_bam_filename="${haplotype2_common_output_bam_filename}" \
            -v haplotype1_unique_output_bam_filename="${haplotype1_unique_output_bam_filename}" \
            -v haplotype2_unique_output_bam_filename="${haplotype2_unique_output_bam_filename}" \
            -v samtools_sort_order_argument="${samtools_sort_order_argument}" \
            -v samtools_threads="${samtools_threads}" \
            '
            BEGIN {
                    # Construct command to read haplotype1_input_bam_filename and haplotype2_input_bam_filename,
                    # sort them by read name and convert them to SAM.
                    read_haplotype1_input_bam_filename_cmd = "samtools sort -@ " samtools_threads " -n -O SAM " haplotype1_input_bam_filename;
                    read_haplotype2_input_bam_filename_cmd = "samtools sort -@ " samtools_threads " -n -O SAM " haplotype2_input_bam_filename;

                    # Construct command to write various haplotype 1 and haplotype 2 related BAM output files.
                    write_haplotype1_all_output_bam_filename_cmd = "samtools sort -@ " samtools_threads samtools_sort_order_argument " -O BAM -o " haplotype1_all_output_bam_filename;
                    write_haplotype2_all_output_bam_filename_cmd = "samtools sort -@ " samtools_threads samtools_sort_order_argument " -O BAM -o " haplotype2_all_output_bam_filename;
                    write_haplotype1_common_output_bam_filename_cmd = "samtools sort -@ " samtools_threads samtools_sort_order_argument " -O BAM -o " haplotype1_common_output_bam_filename;
                    write_haplotype2_common_output_bam_filename_cmd = "samtools sort -@ " samtools_threads samtools_sort_order_argument " -O BAM -o " haplotype2_common_output_bam_filename;
                    write_haplotype1_unique_output_bam_filename_cmd = "samtools sort -@ " samtools_threads samtools_sort_order_argument " -O BAM -o " haplotype1_unique_output_bam_filename;
                    write_haplotype2_unique_output_bam_filename_cmd = "samtools sort -@ " samtools_threads samtools_sort_order_argument " -O BAM -o " haplotype2_unique_output_bam_filename;

                    OFS=FS="\t";

                    # Write BAM header of haplotype 1 input BAM file to the haplotype 1 output BAM files.
                    while ( ( (read_haplotype1_input_bam_filename_cmd | getline) > 0 ) && ( $0 ~ /^@/ ) ) {
                        haplotype1_input_sam_line_number += 1;

                        print $0 | write_haplotype1_all_output_bam_filename_cmd;
                        print $0 | write_haplotype1_common_output_bam_filename_cmd;
                        print $0 | write_haplotype1_unique_output_bam_filename_cmd;
                    }

                    # Save first read of haplotype 1.
                    haplotype1_first_read = $0;

                    # Write BAM header of haplotype 2 input BAM file to the haplotype 2 output BAM files.
                    while ( ( (read_haplotype2_input_bam_filename_cmd | getline) > 0 ) && ( $0 ~ /^@/ ) ) {
                        haplotype2_input_sam_line_number += 1;

                        print $0 | write_haplotype2_all_output_bam_filename_cmd;
                        print $0 | write_haplotype2_common_output_bam_filename_cmd;
                        print $0 | write_haplotype2_unique_output_bam_filename_cmd;
                    }

                    # Save first read of haplotype 2.
                    haplotype2_first_read = $0;

                    haplotype1_input_sam_eof = 0;
                    haplotype2_input_sam_eof = 0;

                    # Initialize variables that keep track of how many reads are in each category to zero.
                    stats["haplotype1__input"] = 0;
                    stats["haplotype1__kept"] = 0;
                    stats["haplotype1__common"] = 0;
                    stats["haplotype1__unique"] = 0;
                    stats["haplotype1__mapq_greater_or_equal_to_2"] = 0;
                    stats["haplotype1__mapq_lower_than_2"] = 0;
                    stats["haplotype2__input"] = 0;
                    stats["haplotype2__kept"] = 0;
                    stats["haplotype2__common"] = 0;
                    stats["haplotype2__unique"] = 0;
                    stats["haplotype2__mapq_greater_or_equal_to_2"] = 0;
                    stats["haplotype2__mapq_lower_than_2"] = 0;
                    stats["haplotype1__mapq_greater_than_haplotype2"] = 0;
                    stats["haplotype2__mapq_greater_than_haplotype1"] = 0;
                    stats["haplotype1_haplotype2__same_mapq"] = 0;
                    stats["haplotype1_haplotype2__same_mapq__same_cigar"] = 0;
                    stats["haplotype1_haplotype2__same_mapq__same_cigar__haplotype1_less_XM"] = 0;
                    stats["haplotype1_haplotype2__same_mapq__same_cigar__haplotype2_less_XM"] = 0;
                    stats["haplotype1_haplotype2__same_mapq__same_cigar__same_XM"] = 0;
                    stats["haplotype1_haplotype2__same_mapq__different_cigar"] = 0;
                    stats["haplotype1_haplotype2__same_mapq__different_cigar__haplotype1_more_cigar_Ms"] = 0;
                    stats["haplotype1_haplotype2__same_mapq__different_cigar__haplotype2_more_cigar_Ms"] = 0;
                    stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__haplotype1_less_cigar_nbr_M_patterns"] = 0;
                    stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__haplotype2_less_cigar_nbr_M_patterns"] = 0;
                    stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns"] = 0;
                    stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns__haplotype1_less_XM"] = 0;
                    stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns__haplotype2_less_XM"] = 0;
                    stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns__same_XM"] = 0;

                    while ( haplotype1_input_sam_eof == 0 && haplotype2_input_sam_eof == 0 ) {
                        # Try to read a line from haplotype 1 BAM file.
                        if (haplotype1_first_read != "") {
                            # First read of haplotype 1 was cached when reading past the BAM header.
                            $0 = haplotype1_first_read;
                            haplotype1_first_read = "";
                        } else if ( (read_haplotype1_input_bam_filename_cmd | getline) > 0 ) {
                            # Read next read of haplotype 1 and store in $0.
                            haplotype1_input_sam_line_number += 1
                        } else {
                            # Reached end of haplotype 1 BAM file.
                            haplotype1_input_sam_eof = 1;
                            break;
                        }

                        if (haplotype1_input_sam_eof == 0) {
                            stats["haplotype1__input"] += 1;

                            # Extract read info.
                            haplotype1["qname"] = $1;
                            haplotype1["rname"] = $3;
                            haplotype1["pos"] = $4;
                            haplotype1["mapq"] = $5;
                            haplotype1["cigar"] = $6;
                            haplotype1["line"] = $0;

                            # Set high value for the XM tag (number of mismatches), in case it is not found.
                            haplotype1["XM"] = 1000000;

                            # Search XM tag.
                            for (i = 12; i <= NF; i++) {
                                if ($i ~ /^XM:i:/) {
                                    # Extract associated number for the XM tag.
                                    haplotype1["XM"] = int(substr($i, 6));

                                    # Stop the loop if we found the tag.
                                    break;
                                }
                            }
                        }


                        # Try to read a line from haplotype 2 BAM file.
                        if (haplotype2_first_read != "") {
                            # First read of haplotype 2 was cached when reading past the BAM header.
                            $0 = haplotype2_first_read;
                            haplotype2_first_read = "";
                        } else if ( (read_haplotype2_input_bam_filename_cmd | getline) > 0 ) {
                            # Read next read of haplotype 2 and store in $0.
                            haplotype2_input_sam_line_number += 1
                        } else {
                            # Reached end of haplotype 2 BAM file.
                            haplotype2_input_sam_eof = 1;
                            break;
                        }

                        if (haplotype2_input_sam_eof == 0) {
                            stats["haplotype2__input"] += 1;

                            # Extract read info.
                            haplotype2["qname"] = $1;
                            haplotype2["rname"] = $3;
                            haplotype2["pos"] = $4;
                            haplotype2["mapq"] = $5;
                            haplotype2["cigar"] = $6;
                            haplotype2["line"] = $0;

                            # Set high value for the XM tag (number of mismatches), in case it is not found.
                            haplotype2["XM"] = 1000000;

                            # Search XM tag.
                            for (i = 12; i <= NF; i++) {
                                if ($i ~ /^XM:i:/) {
                                    # Extract associated number for the XM tag.
                                    haplotype2["XM"] = int(substr($i, 6));

                                    # Stop the loop if we found the tag.
                                    break;
                                }
                            }
                        }


                        if ( haplotype1["qname"] == haplotype2["qname"] ) {
                            # Read names match.

                            if ( (haplotype1["mapq"] >= 2) && (haplotype2["mapq"] >= 2) ) {
                                # Read maps uniquely both in haplotype 1 and 2.
                                stats["haplotype1__mapq_greater_or_equal_to_2"] += 1;
                                stats["haplotype2__mapq_greater_or_equal_to_2"] += 1;

                                if ( haplotype1["mapq"] == haplotype2["mapq"] ) {
                                    # Read maps with same mapping quality to haplotype 1 and 2.
                                    stats["haplotype1_haplotype2__same_mapq"] += 1;

                                    if ( haplotype1["cigar"] == haplotype2["cigar"] ) {
                                        # Read has the same cigar string in haplotype 1 and 2.
                                        stats["haplotype1_haplotype2__same_mapq__same_cigar"] += 1;

                                        if ( haplotype1["XM"] < haplotype2["XM"] ) {
                                            # Read has the less mismatches according to XM tag value in haplotype 1.
                                            stats["haplotype1_haplotype2__same_mapq__same_cigar__haplotype1_less_XM"] += 1;
                                            stats["haplotype1__kept"] += 1;
                                            stats["haplotype1__unique"] += 1;

                                            # Write this read to the haplotype 1 all and unique output BAM files.
                                            print haplotype1["line"] | write_haplotype1_all_output_bam_filename_cmd;
                                            print haplotype1["line"] | write_haplotype1_unique_output_bam_filename_cmd;
                                        } else if ( haplotype1["XM"] > haplotype2["XM"] ) {
                                            # Read has the less mismatches according to XM tag value in haplotype 2.
                                            stats["haplotype1_haplotype2__same_mapq__same_cigar__haplotype2_less_XM"] += 1;
                                            stats["haplotype2__kept"] += 1;
                                            stats["haplotype2__unique"] += 1;

                                            # Write this read to the haplotype 2 all and unique output BAM files.
                                            print haplotype2["line"] | write_haplotype2_all_output_bam_filename_cmd;
                                            print haplotype2["line"] | write_haplotype2_unique_output_bam_filename_cmd;
                                        } else {
                                            # Read has the same XM tag value (number of mismatches) in haplotype 1 and 2.
                                            stats["haplotype1_haplotype2__same_mapq__same_cigar__same_XM"] += 1;
                                            stats["haplotype1__kept"] += 1;
                                            stats["haplotype2__kept"] += 1;
                                            stats["haplotype1__common"] += 1;
                                            stats["haplotype2__common"] += 1;

                                            # Write this read to both haplotype 1 and 2 all and common output BAM files.
                                            print haplotype1["line"] | write_haplotype1_all_output_bam_filename_cmd;
                                            print haplotype2["line"] | write_haplotype2_all_output_bam_filename_cmd;
                                            print haplotype1["line"] | write_haplotype1_common_output_bam_filename_cmd;
                                            print haplotype2["line"] | write_haplotype2_common_output_bam_filename_cmd;
                                        }
                                    } else {
                                        # Read has different CIGAR string in haplotype 1 and 2.
                                        # Check which one has the most bases that match their respective reference to
                                        # decide to which haplotype the read most likely belongs too.
                                        stats["haplotype1_haplotype2__same_mapq__different_cigar"] += 1;

                                        haplotype1["cigar_Ms"] = 0;
                                        haplotype2["cigar_Ms"] = 0;

                                        haplotype1["cigar_nbr_M_patterns"] = 0;
                                        haplotype2["cigar_nbr_M_patterns"] = 0;

                                        cigar_search_start = 1;

                                        # Extract all "<number>M" parts from CIGAR string for haplotype 1 and calculate
                                        # total number of matches.
                                        while ( match(substr(haplotype1["cigar"], cigar_search_start), /[0-9]+M/) != 0 ) {
                                            # Extract number from matched "<number>M" part of CIGAR string and add to
                                            # the total number of matches for haplotype 1.
                                            haplotype1["cigar_Ms"] += int(substr(haplotype1["cigar"], cigar_search_start + RSTART - 1, RLENGTH - 1));

                                            # Set new start position to search for next "<number>M" pattern in CIGAR
                                            # string to after current found "<number>M" part of CIGAR string of
                                            # haplotype 1.
                                            cigar_search_start = cigar_search_start + RSTART + RLENGTH - 1;

                                            # Keep track of the number of "<number>M" patterns found in CIGAR string
                                            # of haplotype 1.
                                            haplotype1["cigar_nbr_M_patterns"] += 1;
                                        }

                                        cigar_search_start = 1;

                                        # Extract all "<number>M" parts from CIGAR string for haplotype 2 and calculate
                                        # total number of matches.
                                        while ( match(substr(haplotype2["cigar"], cigar_search_start), /[0-9]+M/) != 0 ) {
                                            # Extract number from matched "<number>M" part of CIGAR string and add to
                                            # the total number of matches for haplotype 2.
                                            haplotype2["cigar_Ms"] += int(substr(haplotype2["cigar"], cigar_search_start + RSTART - 1, RLENGTH - 1));

                                            # Set new start position to search for next "<number>M" pattern in CIGAR
                                            # string to after current found "<number>M" part of CIGAR string of
                                            # haplotype 2.
                                            cigar_search_start = cigar_search_start + RSTART + RLENGTH - 1;

                                            # Keep track of the number of "<number>M" patterns found in CIGAR string
                                            # of haplotype 2.
                                            haplotype2["cigar_nbr_M_patterns"] += 1
                                        }

                                        if ( haplotype1["cigar_Ms"] > haplotype2["cigar_Ms"] ) {
                                            # Read maps with more matches to haplotype 1.
                                            stats["haplotype1_haplotype2__same_mapq__different_cigar__haplotype1_more_cigar_Ms"] += 1;
                                            stats["haplotype1__kept"] += 1;
                                            stats["haplotype1__unique"] += 1;

                                            # Write this read to the haplotype 1 all and unique output BAM files.
                                            print haplotype1["line"] | write_haplotype1_all_output_bam_filename_cmd;
                                            print haplotype1["line"] | write_haplotype1_unique_output_bam_filename_cmd;
                                        } else if ( haplotype1["cigar_Ms"] < haplotype2["cigar_Ms"] ) {
                                            # Read maps with more matches to haplotype 2.
                                            stats["haplotype1_haplotype2__same_mapq__different_cigar__haplotype2_more_cigar_Ms"] += 1;
                                            stats["haplotype2__kept"] += 1;
                                            stats["haplotype2__unique"] += 1;

                                            # Write this read to the haplotype 2 all and unique output BAM files.
                                            print haplotype2["line"] | write_haplotype2_all_output_bam_filename_cmd;
                                            print haplotype2["line"] | write_haplotype2_unique_output_bam_filename_cmd;
                                        } else if ( haplotype1["cigar_nbr_M_patterns"] < haplotype2["cigar_nbr_M_patterns"] ) {
                                            # Read maps with same number of matches but the least number of M patterns
                                            # (so less insertions and deletions) in haplotype 1.
                                            stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__haplotype1_less_cigar_nbr_M_patterns"] += 1;
                                            stats["haplotype1__kept"] += 1;
                                            stats["haplotype1__unique"] += 1;

                                            # Write this read to the haplotype 1 all and unique output BAM files.
                                            print haplotype1["line"] | write_haplotype1_all_output_bam_filename_cmd;
                                            print haplotype1["line"] | write_haplotype1_unique_output_bam_filename_cmd;
                                        } else if ( haplotype1["cigar_nbr_M_patterns"] > haplotype2["cigar_nbr_M_patterns"] ) {
                                            # Read maps with same number of matches but the least number of M patterns
                                            # (so less insertions and deletions) in haplotype 2.
                                            stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__haplotype2_less_cigar_nbr_M_patterns"] += 1;
                                            stats["haplotype2__kept"] += 1;
                                            stats["haplotype2__unique"] += 1;

                                            # Write this read to the haplotype 2 all and unique output BAM files.
                                            print haplotype2["line"] | write_haplotype2_all_output_bam_filename_cmd;
                                            print haplotype2["line"] | write_haplotype2_unique_output_bam_filename_cmd;
                                        } else if ( haplotype1["XM"] < haplotype2["XM"] ) {
                                            # Read maps with same number of matches and the same number of M patterns
                                            # but with less mismatches (according to XM tag) in haplotype 1.
                                            stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns__haplotype1_less_XM"] += 1;
                                            stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns"] += 1;
                                            stats["haplotype1__kept"] += 1;
                                            stats["haplotype1__unique"] += 1;

                                            # Write this read to the haplotype 1 all and unique output BAM files.
                                            print haplotype1["line"] | write_haplotype1_all_output_bam_filename_cmd;
                                            print haplotype1["line"] | write_haplotype1_unique_output_bam_filename_cmd;
                                        } else if ( haplotype1["XM"] > haplotype2["XM"] ) {
                                            # Read maps with same number of matches and the same number of M patterns
                                            # but with less mismatches (according to XM tag) in haplotype 2.
                                            stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns__haplotype2_less_XM"] += 1;
                                            stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns"] += 1;
                                            stats["haplotype2__kept"] += 1;
                                            stats["haplotype2__unique"] += 1;

                                            # Write this read to the haplotype 1 all and unique output BAM files.
                                            print haplotype1["line"] | write_haplotype1_all_output_bam_filename_cmd;
                                            print haplotype1["line"] | write_haplotype1_unique_output_bam_filename_cmd;
                                        } else {
                                            # Read maps with same number of matches and the same number of M patterns
                                            # (but CIGAR strings) and same number of mismatches (according to XM tag)
                                            # in both haplotypes ==> ambiguous.
                                            stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns__same_XM"] += 1
                                            stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns"] += 1;
                                            stats["haplotype1__kept"] += 1;
                                            stats["haplotype2__kept"] += 1;
                                            stats["haplotype1__common"] += 1;
                                            stats["haplotype2__common"] += 1;

                                            # Write this read to the haplotype 1 and 2 all and common output BAM files.
                                            print haplotype1["line"] | write_haplotype1_all_output_bam_filename_cmd;
                                            print haplotype2["line"] | write_haplotype2_all_output_bam_filename_cmd;
                                            print haplotype1["line"] | write_haplotype1_common_output_bam_filename_cmd;
                                            print haplotype2["line"] | write_haplotype2_common_output_bam_filename_cmd;
                                        }

                                    }
                                } else if ( haplotype1["mapq"] > haplotype2["mapq"] ) {
                                    # Read maps with higher mapping quality to haplotype 1.
                                    stats["haplotype1__mapq_greater_than_haplotype2"] += 1;
                                    stats["haplotype1__kept"] += 1;
                                    stats["haplotype1__unique"] += 1;

                                    # Write this read to the haplotype 1 all and unique output BAM files.
                                    print haplotype1["line"] | write_haplotype1_all_output_bam_filename_cmd;
                                    print haplotype1["line"] | write_haplotype1_unique_output_bam_filename_cmd;
                                } else {
                                    # Read maps with higher mapping quality to haplotype 2.
                                    stats["haplotype2__mapq_greater_than_haplotype1"] += 1;
                                    stats["haplotype2__kept"] += 1;
                                    stats["haplotype2__unique"] += 1;

                                    # Write this read to the haplotype 2 all and unique output BAM files.
                                    print haplotype2["line"] | write_haplotype2_all_output_bam_filename_cmd;
                                    print haplotype2["line"] | write_haplotype2_unique_output_bam_filename_cmd;
                                }
                            } else if (haplotype1["mapq"] >= 2) {
                                # Read maps only uniquely in haplotype 1.
                                stats["haplotype1__mapq_greater_or_equal_to_2"] += 1;
                                stats["haplotype2__mapq_lower_than_2"] += 1;
                                stats["haplotype1__kept"] += 1;
                                stats["haplotype1__unique"] += 1;

                                # Write this read to the haplotype 1 all and unique output BAM files.
                                print haplotype1["line"] | write_haplotype1_all_output_bam_filename_cmd;
                                print haplotype1["line"] | write_haplotype1_unique_output_bam_filename_cmd;
                            } else if (haplotype2["mapq"] >= 2) {
                                # Read maps only uniquely in haplotype 2.
                                stats["haplotype2__mapq_greater_or_equal_to_2"] += 1;
                                stats["haplotype1__mapq_lower_than_2"] += 1;
                                stats["haplotype2__kept"] += 1;
                                stats["haplotype2__unique"] += 1;

                                # Write this read to the haplotype 2 all and unique output BAM files.
                                print haplotype2["line"] | write_haplotype2_all_output_bam_filename_cmd;
                                print haplotype2["line"] | write_haplotype2_unique_output_bam_filename_cmd;
                            }
                        } else {
                            print "Error: haplotype 1 and 2 BAM files are not sorted by read name:";
                            printf("  - haplotype 1 (line number: %d): %s\n", haplotype1_input_sam_line_number, haplotype1["line"]);
                            printf("  - haplotype 2 (line number: %d): %s\n", haplotype2_input_sam_line_number, haplotype2["line"]);
                        }
                    }

                    # Close file handles.
                    close(read_haplotype1_input_bam_filename_cmd);
                    close(read_haplotype2_input_bam_filename_cmd);
                    close(write_haplotype1_all_output_read_name_sorted_bam_filename_cmd);
                    close(write_haplotype2_all_output_read_name_sorted_bam_filename_cmd);
                    close(write_haplotype1_common_output_bam_filename_cmd);
                    close(write_haplotype2_common_output_bam_filename_cmd);
                    close(write_haplotype1_unique_output_bam_filename_cmd);
                    close(write_haplotype2_unique_output_bam_filename_cmd);

                    # Print some statistics.
                    print "Statistics:";
                    print "haplotype1__input\t" stats["haplotype1__input"];
                    print "haplotype1__kept\t" stats["haplotype1__kept"];
                    print "haplotype1__common\t" stats["haplotype1__common"];
                    print "haplotype1__unique\t" stats["haplotype1__unique"];
                    print "haplotype1__mapq_greater_or_equal_to_2\t" stats["haplotype1__mapq_greater_or_equal_to_2"];
                    print "haplotype1__mapq_lower_than_2\t" stats["haplotype1__mapq_lower_than_2"];
                    print "haplotype2__input\t" stats["haplotype2__input"];
                    print "haplotype2__kept\t" stats["haplotype2__kept"];
                    print "haplotype2__common\t" stats["haplotype2__common"];
                    print "haplotype2__unique\t" stats["haplotype2__unique"];
                    print "haplotype2__mapq_greater_or_equal_to_2\t" stats["haplotype2__mapq_greater_or_equal_to_2"];
                    print "haplotype2__mapq_lower_than_2\t" stats["haplotype2__mapq_lower_than_2"];
                    print "haplotype1__mapq_greater_than_haplotype2\t" stats["haplotype1__mapq_greater_than_haplotype2"];
                    print "haplotype2__mapq_greater_than_haplotype1\t" stats["haplotype2__mapq_greater_than_haplotype1"];
                    print "haplotype1_haplotype2__same_mapq\t" stats["haplotype1_haplotype2__same_mapq"];
                    print "haplotype1_haplotype2__same_mapq__same_cigar\t" stats["haplotype1_haplotype2__same_mapq__same_cigar"];
                    print "haplotype1_haplotype2__same_mapq__same_cigar__haplotype1_less_XM\t" stats["haplotype1_haplotype2__same_mapq__same_cigar__haplotype1_less_XM"];
                    print "haplotype1_haplotype2__same_mapq__same_cigar__haplotype2_less_XM\t" stats["haplotype1_haplotype2__same_mapq__same_cigar__haplotype2_less_XM"];
                    print "haplotype1_haplotype2__same_mapq__same_cigar__same_XM\t" stats["haplotype1_haplotype2__same_mapq__same_cigar__same_XM"];
                    print "haplotype1_haplotype2__same_mapq__different_cigar\t" stats["haplotype1_haplotype2__same_mapq__different_cigar"];
                    print "haplotype1_haplotype2__same_mapq__different_cigar__haplotype1_more_cigar_Ms\t" stats["haplotype1_haplotype2__same_mapq__different_cigar__haplotype1_more_cigar_Ms"];
                    print "haplotype1_haplotype2__same_mapq__different_cigar__haplotype2_more_cigar_Ms\t" stats["haplotype1_haplotype2__same_mapq__different_cigar__haplotype2_more_cigar_Ms"];
                    print "haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__haplotype1_less_cigar_nbr_M_patterns\t" stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__haplotype1_less_cigar_nbr_M_patterns"];
                    print "haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__haplotype2_less_cigar_nbr_M_patterns\t" stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__haplotype2_less_cigar_nbr_M_patterns"];
                    print "haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns\t" stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns"];
                    print "haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns__haplotype1_less_XM\t" stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns__haplotype1_less_XM"];
                    print "haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns__haplotype2_less_XM\t" stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns__haplotype2_less_XM"];
                    print "haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns__same_XM (ambiguous)\t" stats["haplotype1_haplotype2__same_mapq__different_cigar__same_cigar_Ms__same_cigar_nbr_M_patterns__same_XM"];
            }
    '

    if [ "${sorted_type}" = "pos_sorted" ] ; then
        printf 'Indexing BAM files...';

        samtools index "${haplotype1_all_output_bam_filename}";
        samtools index "${haplotype2_all_output_bam_filename}";
        samtools index "${haplotype1_common_output_bam_filename}";
        samtools index "${haplotype2_common_output_bam_filename}";
        samtools index "${haplotype1_unique_output_bam_filename}";
        samtools index "${haplotype2_unique_output_bam_filename}";
    fi
}



create_haplotype_specific_bam_files_from_bowtie2_haplotype_mapped_bam_files "${@}";
