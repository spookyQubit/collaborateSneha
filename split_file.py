import os
import pandas as pd

g_data_dir = "./data"

g_treatment_parent_file = os.path.join(g_data_dir, "Galaxy_si_bedgraph.txt")
g_treatment_template_file = os.path.join(g_data_dir, "si_data/Galaxy_si_bedgraph_%s.txt")

g_intron_parent_file = os.path.join(g_data_dir, "Galaxy_intron_gtf.txt")
g_intron_template_file = os.path.join(g_data_dir, "intron_data/Galaxy_intron_gtf_%s.txt")

g_intron_scores_template_file = os.path.join(g_data_dir, "intron_si_scores_data/Galaxy_intron_scores_gtf_%s.txt")


def split_intron_files_by_chromosome(parent_file):
    df = pd.read_csv(parent_file, sep='\t')
    chromosomes = df['chr'].unique()
    files_per_chromosome = {}
    for chromosome in chromosomes:
        file_to_write_to = g_intron_template_file % chromosome
        if not os.path.isfile(file_to_write_to):
            df[df['chr']==chromosome].to_csv(file_to_write_to, sep='\t')
        files_per_chromosome[chromosome] = file_to_write_to
    return files_per_chromosome


def split_treatment_files_by_chromosome(parent_file, chromosomes_in_intron_data):
    df = pd.read_csv(parent_file, sep='\t')
    files_per_chromosome = {}
    for chromosome in chromosomes_in_intron_data:
        file_to_write_to = g_treatment_template_file % chromosome
        if not os.path.isfile(file_to_write_to):
            df[df['chr']==chromosome].to_csv(file_to_write_to, sep='\t')
        files_per_chromosome[chromosome] = file_to_write_to
    return files_per_chromosome


def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def get_intron_sum(df_treatment_per_chr, start_idx, end_idx):

    df_overlap = df_treatment_per_chr[((df_treatment_per_chr['end'] >= start_idx) & (df_treatment_per_chr['end'] <= end_idx))
                    | ((df_treatment_per_chr['start'] >= start_idx) & (df_treatment_per_chr['start'] <= end_idx))
                    | ((df_treatment_per_chr['start'] <= start_idx) & (df_treatment_per_chr['end'] >= end_idx))]

    count_average = 0
    length_of_introns_covered = 0
    for _, row in df_overlap.iterrows():
        overlap_window = getOverlap([start_idx, end_idx], [row['start'], row['end']])
        count_average += (row['count'] * overlap_window) / (row['end'] - row['start'])
        length_of_introns_covered += overlap_window

    if length_of_introns_covered == 0:
        return 0
    else:
        return count_average/length_of_introns_covered


def write_intron_scores_data(intron_data_file, treatment_data_file, file_to_write):

    print("intron_data_file = {}, treatment_data_file = {}".format(intron_data_file, treatment_data_file))

    df_treatment = pd.read_csv(treatment_data_file, sep='\t')
    with open(intron_data_file) as f1, \
            open(file_to_write, 'w') as f2:

        next(f1)  # skip the first line because it is header

        lines = f1.read().splitlines()
        i = 0
        for line in lines:
            idx_, chr_intron, bed_, intron_, \
            start_intron, end_intron, zero_, sign_, do_, \
            gene_name = line.split('\t')
            intron_sum = get_intron_sum(df_treatment,
                                        int(start_intron),
                                        int(end_intron))

            line_to_write = '\t'.join(
                [chr_intron, bed_, intron_, start_intron, end_intron, zero_, sign_, do_, gene_name, str(intron_sum)])
            line_to_write += '\n'

            f2.write(line_to_write)
            i += 1
            if i % 1000 == 0:
                print("Computed {} lines".format(i))


def main():
    intron_data_files_per_chromosome = split_intron_files_by_chromosome(g_intron_parent_file)
    chromosomes_in_intron_data = [k for k,v in intron_data_files_per_chromosome.items()]
    treatment_data_files_per_chromosome = split_treatment_files_by_chromosome(g_treatment_parent_file, chromosomes_in_intron_data)

    files_per_chromosome = {}
    for chr in intron_data_files_per_chromosome:
        files_per_chromosome[chr] = (intron_data_files_per_chromosome[chr], treatment_data_files_per_chromosome[chr])

    for chr in files_per_chromosome:
        intron_data_file, treatment_data_file = files_per_chromosome[chr]
        write_intron_scores_data(intron_data_file, treatment_data_file, g_intron_scores_template_file%chr)


if __name__=="__main__":
    main()