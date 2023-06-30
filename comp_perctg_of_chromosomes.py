import pysam

# Original bam file
original_bam_file = "ont_ill_spades_contigs2ref.sorted.bam"

# Chromosome number
chromosome_count = 15

# Selecting the 15 longest chromosomes in the bam file and writing to the new bam file
with pysam.AlignmentFile(original_bam_file, "rb") as original_bam:
    chromosome_lengths = {}

    for contig in original_bam.header.references:
        chromosome_lengths[contig] = original_bam.get_reference_length(contig)

    sorted_chromosomes = sorted(chromosome_lengths, key=chromosome_lengths.get, reverse=True)[:chromosome_count]

    # Creating new bam file and calculating percent complete while writing
    with pysam.AlignmentFile("new_contigs.bam", "wb", header=original_bam.header) as new_bam:
        total_completed_percentages = 0.0
        completed_percentages = []
        for contig in sorted_chromosomes:
            total_bases = 0
            completed_bases = 0

            for read in original_bam.fetch(contig):
                total_bases += read.reference_length
                if read.query_alignment_start is not None and read.query_alignment_end is not None:
                    completed_bases += read.query_alignment_end - read.query_alignment_start

            if total_bases > 0:
                completed_percentage = (completed_bases / total_bases) * 100
            else:
                completed_percentage = 0.0

            completed_percentages.append(completed_percentage)
            total_completed_percentages += completed_percentage

            print(f"Completion percentage ({contig}): {completed_percentage:.2f}%")

        average_completed_percentage = total_completed_percentages / len(sorted_chromosomes)
        print(f"\nAverage completion percentage: {average_completed_percentage:.2f}%")
