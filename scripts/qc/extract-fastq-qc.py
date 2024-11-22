import json
import sys

sys.stderr = open(snakemake.log[0], "a")
json_file  = open(snakemake.input.rep_json, 'r')
data       = json.load(json_file)
data_sum   = data['summary']['before_filtering']
json_file.close()

reads1       = snakemake.params.read_filenames['r1']
reads2       = snakemake.params.read_filenames['r2']
total_reads  = str(data_sum['total_reads'])
total_r1     = str(data['read1_before_filtering']['total_reads'])
total_r2     = str(data['read2_before_filtering']['total_reads'])
read_lengths = str(data_sum['read1_mean_length']) + "-" + str(data_sum['read2_mean_length'])
insert_peak  = str(data["insert_size"]["peak"])
gccont       = str(round(data_sum['gc_content'] * 100, 2))
bases20      = str(round(data_sum['q20_rate'] * 100, 2))
bases30      = str(round(data_sum['q30_rate'] * 100, 2))
lowqual      = str(round(data['filtering_result']['low_quality_reads'] / data_sum['total_reads'] * 100, 2))
dups         = str(round(data['duplication']['rate'] * 100, 2))

report = f"""## Sequencing Report
File names:
    - Read1: {reads1}
    - Read2: {reads2}
Sample: {snakemake.params.sample_name}
File type: {snakemake.params.file_type}
Encoding: {snakemake.params.encoding}
Total reads: {total_reads}
Read1: {total_r1}
Read2: {total_r2}
Mean read length: {read_lengths}
Insert size peak: {insert_peak}
GC content [%]: {gccont}
Basecalls with Q>20 [%]: {bases20}
Basecalls with Q>30 [%]: {bases30}
Low quality reads [%]: {lowqual}
Duplicated reads [%]: {dups}
"""

out_rep = open(snakemake.output.stats_json, "w")
out_rep.writelines(report)
out_rep.close()

