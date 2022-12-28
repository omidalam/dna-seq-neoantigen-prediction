import sys
import pandas as pd

sys.stderr = open(snakemake.log[0], "w")

data = pd.read_csv(snakemake.input.tsv, sep="\t")
xlWriter = pd.ExcelWriter(snakemake.output.xlsx, engine='xlsxwriter')
data.to_excel(xlWriter, sheet_name='Sheet1')
xlWriter.save()