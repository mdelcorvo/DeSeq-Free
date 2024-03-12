import pandas as pd

#### Define input and output tables ####
s1_report = snakemake.input.report
s2_stats = snakemake.input.stats
s3_sum_stats = snakemake.input.sum_stats
s4_mosdepth = snakemake.input.mosdepth
output = snakemake.output.summary_report


#### Import tables and clean fields when needed ####
s1 = pd.read_table(s1_report, delimiter=':',names=["", "Value"])
s2 = pd.read_table(s2_stats, names=["Pairs statistics", "Value"])
s3 = pd.read_table(s3_sum_stats, names=["Pairs statistics", "Value", "Percentage"])
s4 = pd.read_table(s4_mosdepth)

#### Write sheets to Excel object ####
writer = pd.ExcelWriter(output, engine='xlsxwriter')

s1.to_excel(writer, sheet_name='QC sequencing', index=False)
s3.to_excel(writer, sheet_name='Proximity-ligation summary', index=False)
s2.to_excel(writer, sheet_name='Proximity-ligation assessment', index=False)
s4.to_excel(writer, sheet_name='Coverage summary', index=False)


worksheet = writer.sheets['QC sequencing']
worksheet.set_column(0, 1, 30)

worksheet1 = writer.sheets['Proximity-ligation summary']
worksheet1.set_column(0, 0, 40)
worksheet1.set_column(1, 2, 20)

worksheet2 = writer.sheets['Proximity-ligation assessment']
worksheet2.set_column(0, 1, 30)

writer.close()
