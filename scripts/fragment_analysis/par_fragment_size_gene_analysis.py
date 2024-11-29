import pysam
import numpy as np
import pandas as pd
import argparse
import HTSeq
from collections import defaultdict
from HTSeq import GFF_Reader, GenomicArrayOfSets
from multiprocessing import Pool, cpu_count
import matplotlib.pyplot as plt


def load_gene_annotations(gtf_file):
    """
    Load gene annotations from a GTF file into a GenomicArrayOfSets.
    """
    gtf_reader = GFF_Reader(gtf_file)
    ga = GenomicArrayOfSets("auto", stranded=False)
    gene_lengths = defaultdict(int)

    for feature in gtf_reader:
        if feature.type == "exon":
            ga[feature.iv] += feature.attr["gene_id"]
            gene_lengths[feature.attr["gene_id"]] += feature.iv.end - feature.iv.start

    return ga, gene_lengths


def process_region(region, bam_file, ga):
    """
    Process a specific region of the BAM file and extract fragment sizes mapped to genes,
    ensuring at least 50% overlap of the fragment with the gene.
    """
    bam = pysam.AlignmentFile(bam_file, "rb")
    gene_fragment_sizes = defaultdict(list)

    for read in bam.fetch(region=region):
        if read.is_proper_pair and not read.is_unmapped and not read.mate_is_unmapped:
            fragment_size = abs(read.template_length)
            fragment_start = read.reference_start
            fragment_end = fragment_start + fragment_size

            # Validate the interval
            if fragment_start >= fragment_end:  # Invalid interval
                continue

            try:
                fragment_interval = HTSeq.GenomicInterval(
                    read.reference_name, fragment_start, fragment_end
                )
            except ValueError:  # Catch intervals outside valid ranges
                continue

            # Check overlap with genes
            for iv, gene_ids in ga[fragment_interval].steps():
                overlap_length = min(fragment_interval.end, iv.end) - max(fragment_interval.start, iv.start)
                if overlap_length > 0 and (overlap_length / fragment_size) >= 0.5:
                    for gene_id in gene_ids:
                        gene_fragment_sizes[gene_id].append(fragment_size)

    bam.close()
    return gene_fragment_sizes



def merge_fragment_sizes(results):
    """
    Merge fragment size results from multiple processes.
    """
    merged = defaultdict(list)
    for result in results:
        for gene_id, sizes in result.items():
            merged[gene_id].extend(sizes)
    return merged


def categorize_size(size):
    """
    Categorize fragment size into ranges.
    """
    if 100 <= size < 200:
        return "100-200"
    elif 200 <= size < 500:
        return "200-500"
    elif 500 <= size < 1000:
        return "500-1000"
    elif size >= 1000:
        return ">1000"
    return "<100"


def summarize_fragment_sizes(gene_fragment_sizes, gene_lengths, sample_name, total_fragments):
    """
    Summarize fragment sizes per gene for a single BAM file, including normalization by gene length,
    total fragments, calculating normalized counts within size ranges, and calculating percentages
    for each size range based on normalized counts.
    """
    gene_stats = []
    for gene_id, sizes in gene_fragment_sizes.items():
        mean_size = np.mean(sizes)
        median_size = np.median(sizes)
        std_dev = np.std(sizes)
        num_fragments = len(sizes)

        # Normalization
        gene_length = gene_lengths.get(gene_id, 1)  # Avoid division by zero
        normalized_count = num_fragments / (total_fragments * gene_length)

        # Size range counts
        size_ranges = pd.Series(sizes).apply(categorize_size).value_counts()

        # Normalize counts within each range
        normalized_ranges = {
            range: (size_ranges.get(range, 0) / (gene_length * total_fragments))
            for range in ["100-200", "200-500", "500-1000", ">1000"]
        }

        # Total normalized counts for the gene
        total_normalized = sum(normalized_ranges.values())

        # Calculate percentages based on normalized counts
        percentage_normalized_ranges = {
            f"%_normalized_{range}": (normalized_ranges[range] / total_normalized * 100)
            if total_normalized > 0 else 0
            for range in ["100-200", "200-500", "500-1000", ">1000"]
        }

        gene_stats.append({
            "sample_name": sample_name,
            "gene_id": gene_id,
            "mean_fragment_size": mean_size,
            "median_fragment_size": median_size,
            "std_dev_fragment_size": std_dev,
            "num_fragments": num_fragments,
            "normalized_count": normalized_count,
            **normalized_ranges,
            **percentage_normalized_ranges
        })

    return pd.DataFrame(gene_stats)


def process_bam_file(args):
    bam_file, sample_name, ga, gene_lengths, output_prefix, processes = args
    print(f"Processing sample {sample_name} from BAM file: {bam_file}")

    # Process BAM file in parallel by regions
    bam = pysam.AlignmentFile(bam_file, "rb")
    regions = bam.references
    total_fragments = sum([1 for read in bam.fetch() if read.is_proper_pair])
    bam.close()

    with Pool(processes=processes) as pool:
        results = pool.starmap(process_region, [(region, bam_file, ga) for region in regions])

    # Merge fragment sizes and summarize
    gene_fragment_sizes = merge_fragment_sizes(results)
    sample_df = summarize_fragment_sizes(gene_fragment_sizes, gene_lengths, sample_name, total_fragments)

    # Save individual sample results
    sample_df.to_csv(f"{output_prefix}_{sample_name}_gene_fragment_sizes.csv", index=False)
    print(f"Results saved for sample {sample_name} to {output_prefix}_{sample_name}_gene_fragment_sizes.csv")

    return sample_df


def main():
    parser = argparse.ArgumentParser(description="Analyze fragment sizes across multiple BAM files.")
    parser.add_argument("-l", "--bam_list", required=True, help="Text file with BAM file paths and sample names (tab-separated).")
    parser.add_argument("-g", "--gtf", required=True, help="GTF file with gene annotations.")
    parser.add_argument("-o", "--output", required=True, help="Output prefix for results.")
    parser.add_argument("-p", "--processes", type=int, default=cpu_count(), help="Number of parallel processes (default: all available cores).")

    args = parser.parse_args()

    # Step 1: Load gene annotations
    print("Loading gene annotations...")
    ga, gene_lengths = load_gene_annotations(args.gtf)

    # Step 2: Read BAM file list
    bam_files = []
    sample_names = []
    with open(args.bam_list, "r") as f:
        for line in f:
            bam_file, sample_name = line.strip().split("\t")
            bam_files.append(bam_file)
            sample_names.append(sample_name)

    # Step 3: Process each BAM file
    all_results = []
    for bam_file, sample_name in zip(bam_files, sample_names):
        sample_result = process_bam_file((bam_file, sample_name, ga, gene_lengths, args.output, args.processes))
        all_results.append(sample_result)

    # Step 4: Merge results from all samples
    print("Merging results from all samples...")
    merged_df = pd.concat(all_results)
    merged_df.to_csv(f"{args.output}_all_samples_gene_fragment_sizes.csv", index=False)
    print(f"Merged results saved to {args.output}_all_samples_gene_fragment_sizes.csv")


if __name__ == "__main__":
    main()
