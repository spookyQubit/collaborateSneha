import pandas as pd

intron_coordintes_file = "Galaxy_intron_gtf.txt"
dmso_file = "Galaxy_dmso_bedgraph.txt"

dmso = pd.read_csv(dmso_file, sep='\t')
print(dmso.head())





print(dmso[(dmso['start'] > 10000)
           & (dmso['end'] < 10061)
           & (dmso['chr'] == 'chr1')]['dmso_count'].sum())