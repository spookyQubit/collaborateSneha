import os
import pandas as pd

g_data_dir = "./data_exon"

g_treatment_parent_file = os.path.join(g_data_dir, "Galaxy_dmso_bedgraph.txt")
g_treatment_template_file = os.path.join(g_data_dir, "dmso_data/Galaxy_dmso_bedgraph_%s.txt")

g_feature_parent_file = os.path.join(g_data_dir, "Galaxy_exon_gtf.txt")
g_feature_template_file = os.path.join(g_data_dir, "exon_data/Galaxy_exon_gtf_%s.txt")

g_feature_scores_template_file = os.path.join(g_data_dir, "exon_dmso_scores_data/Galaxy_exon_scores_gtf_%s.txt")


def split_feature_files_by_chromosome(parent_file):
    df = pd.read_csv(parent_file, sep='\t')
    chromosomes = df['chr'].unique()
    files_per_chromosome = {}
    for chromosome in chromosomes:
        file_to_write_to = g_feature_template_file % chromosome
        if not os.path.isfile(file_to_write_to):
            df[df['chr']==chromosome].to_csv(file_to_write_to, sep='\t')
        files_per_chromosome[chromosome] = file_to_write_to
    return files_per_chromosome


def split_treatment_files_by_chromosome(parent_file, chromosomes_in_feature_data):
    df = pd.read_csv(parent_file, sep='\t')
    files_per_chromosome = {}
    for chromosome in chromosomes_in_feature_data:
        file_to_write_to = g_treatment_template_file % chromosome
        if not os.path.isfile(file_to_write_to):
            df[df['chr']==chromosome].to_csv(file_to_write_to, sep='\t')
        files_per_chromosome[chromosome] = file_to_write_to
    return files_per_chromosome


def getOverlap(a, b):
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def get_feature_sum(df_treatment_per_chr, start_idx, end_idx):

    df_overlap = df_treatment_per_chr[((df_treatment_per_chr['end'] >= start_idx) & (df_treatment_per_chr['end'] <= end_idx))
                    | ((df_treatment_per_chr['start'] >= start_idx) & (df_treatment_per_chr['start'] <= end_idx))
                    | ((df_treatment_per_chr['start'] <= start_idx) & (df_treatment_per_chr['end'] >= end_idx))]

    count_average = 0
    length_of_features_covered = 0
    for _, row in df_overlap.iterrows():
        overlap_window = getOverlap([start_idx, end_idx], [row['start'], row['end']])
        count_average += (row['count'] * overlap_window) / (row['end'] - row['start'])
        length_of_features_covered += overlap_window

    if length_of_features_covered == 0:
        return 0
    else:
        return count_average/length_of_features_covered


def write_feature_scores_data(feature_data_file, treatment_data_file, file_to_write):

    print("feature_data_file = {}, treatment_data_file = {}".format(feature_data_file, treatment_data_file))

    df_treatment = pd.read_csv(treatment_data_file, sep='\t')
    with open(feature_data_file) as f1, \
            open(file_to_write, 'w') as f2:

        next(f1)  # skip the first line because it is header

        lines = f1.read().splitlines()
        i = 0
        for line in lines:

            line_data = line.split('\t')
            start_feature = int(line_data[4])
            end_feature = int(line_data[5])

            feature_sum = get_feature_sum(df_treatment,
                                          start_feature,
                                          end_feature)

            line_data.append(str(feature_sum))

            line_to_write = '\t'.join(line_data)
            line_to_write += '\n'

            f2.write(line_to_write)
            i += 1
            if i % 1000 == 0:
                print("Computed {} lines".format(i))


def main():
    feature_data_files_per_chromosome = split_feature_files_by_chromosome(g_feature_parent_file)
    chromosomes_in_feature_data = [k for k,v in feature_data_files_per_chromosome.items()]
    treatment_data_files_per_chromosome = split_treatment_files_by_chromosome(g_treatment_parent_file, chromosomes_in_feature_data)

    files_per_chromosome = {}
    for chr in feature_data_files_per_chromosome:
        files_per_chromosome[chr] = (feature_data_files_per_chromosome[chr], treatment_data_files_per_chromosome[chr])

    for chr in files_per_chromosome:
        feature_data_file, treatment_data_file = files_per_chromosome[chr]
        write_feature_scores_data(feature_data_file, treatment_data_file, g_feature_scores_template_file%chr)


if __name__=="__main__":
    main()
