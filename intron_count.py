import pandas as pd

intron_coordintes_file = "data/Galaxy_intron_gtf.txt"
intron_coordintes_scores_file = "Galaxy_intron_scores_gtf.txt"
dmso_file = "data/Galaxy_dmso_bedgraph.txt"

dmso = pd.read_csv(dmso_file, sep='\t')


#def get_intron_sum(chr_intron, start_intron, end_intron):
#    return dmso[(dmso['start'] >= start_intron)
#                & (dmso['end'] <= end_intron)
#                & (dmso['chr'] == chr_intron)]['dmso_count'].sum()


def get_intron_sum(df_chr, start_intron, end_intron):
    return df_chr[(df_chr['start'] >= start_intron)
                & (df_chr['end'] <= end_intron)]['dmso_count'].sum()


with open(intron_coordintes_file) as f1, \
        open(intron_coordintes_scores_file, 'w') as f2:
    lines = f1.read().splitlines()
    i = 0
    for line in lines:
        chr_intron, bed_, intron_, \
        start_intron, end_intron, zero, sign_, do_, \
        gene_name = line.split('\t')

        intron_sum = get_intron_sum(dmso[dmso['chr'] == chr_intron],
                                    int(start_intron),
                                    int(end_intron))
        line_to_write = '\t'.join([chr_intron, bed_, intron_, start_intron, end_intron, zero, sign_, do_, gene_name, str(intron_sum)])
        line_to_write += '\n'

        f2.write(line_to_write)
        i += 1
        #if i==5:
        #    break
        if i%1000 == 0:
            print("Computed {} lines".format(i))


