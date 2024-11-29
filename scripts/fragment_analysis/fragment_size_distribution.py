import argparse
import os
import pysam
import pandas as pd
import matplotlib.pyplot as plt

def read_bam_list(file_path):
    """
    Read a list of BAM file paths from a text file.

    Args:
        file_path (str): Path to the text file.

    Returns:
        list: List of BAM file paths.
    """
    with open(file_path, "r") as f:
        return [line.strip() for line in f if line.strip()]

def extract_fragment_sizes(bam_file):
    """
    Extract fragment sizes from a BAM file.

    Args:
        bam_file (str): Path to the BAM file.

    Returns:
        list: List of fragment sizes (absolute values of TLEN).
    """
    fragment_sizes = []
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for read in bam.fetch():
            if read.is_proper_pair and not read.is_unmapped and not read.mate_is_unmapped:
                fragment_sizes.append(abs(read.template_length))
    return fragment_sizes

def categorize_fragment_size(size):
    """
    Categorize fragment size into predefined ranges.

    Args:
        size (int): Fragment size.

    Returns:
        str: Fragment size range category.
    """
    if size <= 200:
        return "100-200"
    elif size <= 500:
        return "200-500"
    elif size <= 1000:
        return "500-1000"
    else:
        return ">1000"

def summarize_fragment_sizes_by_size(fragment_sizes, group_name, bam_file):
    """
    Summarize fragment sizes by individual size and categorize into ranges.

    Args:
        fragment_sizes (list): List of fragment sizes.
        group_name (str): Name of the group (e.g., 'Plasma' or 'Tumor').
        bam_file (str): Name of the BAM file.

    Returns:
        pd.DataFrame: Summary DataFrame with counts and ranges for each fragment size.
    """
    summary = pd.Series(fragment_sizes).value_counts().reset_index()
    summary.columns = ["Fragment Size", "Count"]
    summary["Range"] = summary["Fragment Size"].apply(categorize_fragment_size)
    summary["Group"] = group_name
    summary["BAM File"] = os.path.basename(bam_file)
    return summary

def analyze_bams(bam_files, group_name):
    """
    Analyze a group of BAM files for fragment size distribution.

    Args:
        bam_files (list): List of BAM file paths.
        group_name (str): Name of the group (e.g., 'Plasma' or 'Tumor').

    Returns:
        pd.DataFrame: Combined summary DataFrame for all BAM files in the group.
    """
    summaries = []
    for bam_file in bam_files:
        print(f"Processing {bam_file}...")
        fragment_sizes = extract_fragment_sizes(bam_file)
        summary = summarize_fragment_sizes_by_size(fragment_sizes, group_name, bam_file)
        summaries.append(summary)
    return pd.concat(summaries, ignore_index=True)

def plot_fragment_size_distribution(df, output_file):
    """
    Plot fragment size distribution as line plots for each group.

    Args:
        df (pd.DataFrame): DataFrame containing fragment size summaries.
        output_file (str): Path to save the plot.
    """
    pivot_table = df.pivot_table(
        index="Fragment Size", columns="Group", values="Count", aggfunc="sum", fill_value=0
    )
    pivot_table.plot(kind="line", figsize=(12, 6), colormap="viridis")
    plt.title("Fragment Size Distribution by Group")
    plt.xlabel("Fragment Size (bp)")
    plt.ylabel("Frequency")
    plt.grid(axis="y", linestyle="--", alpha=0.7)
    plt.legend(title="Group")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

def main():
    parser = argparse.ArgumentParser(description="Analyze fragment size distributions from BAM files.")
    parser.add_argument(
        "--plasma_bam_list", type=str, required=True,
        help="Path to a text file containing the list of BAM files for the plasma cohort."
    )
    parser.add_argument(
        "--tumor_bam_list", type=str, required=True,
        help="Path to a text file containing the list of BAM files for the tumor cohort."
    )
    parser.add_argument(
        "--output_csv", type=str, required=True,
        help="Path to save the output CSV file."
    )
    parser.add_argument(
        "--output_plot", type=str, required=True,
        help="Path to save the output plot."
    )
    args = parser.parse_args()

    # Read BAM file paths
    plasma_bams = read_bam_list(args.plasma_bam_list)
    tumor_bams = read_bam_list(args.tumor_bam_list)

    # Analyze both groups
    plasma_df = analyze_bams(plasma_bams, "Plasma")
    tumor_df = analyze_bams(tumor_bams, "Tumor")

    # Combine data and save as CSV
    combined_df = pd.concat([plasma_df, tumor_df], ignore_index=True)
    combined_df.to_csv(args.output_csv, index=False)
    print(f"Results saved to {args.output_csv}")

    # Plot and save the distribution
    plot_fragment_size_distribution(combined_df, args.output_plot)
    print(f"Plot saved to {args.output_plot}")

if __name__ == "__main__":
    main()
