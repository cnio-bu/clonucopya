import seaborn as sns
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

if __name__ == '__main__':

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--samples", nargs='+', action='append')
    input_parser.add_argument("--results_mutations", action='store')
    input_parser.add_argument("--pandrugs_dir", action='store')
    input_parser.add_argument("--heatmap_dir", action='store')
    input_parser.add_argument("--out_file_gd", action='store')
    input_parser.add_argument("--out_file_vars", action='store')
    args = input_parser.parse_args()

    # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/pandrugs'
    pandrugs_dir = args.pandrugs_dir

    # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/heatmap'
    heatmap_dir = args.heatmap_dir

    # i.e: /home/cleon/PycharmProjects/MSc-Project/out/mapscape/933124/results_mutations.tsv file
    results_mutations_tsv = pd.read_csv(args.results_mutations, sep=',')

    # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/heatmap/933124/heatmap_all_compared_aml_patient.png'
    out_file_gd = args.out_file_gd

    # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/heatmap/933124/heatmap_all_vars_compared_aml_patient.png'
    out_file_vars = args.out_file_vars

    arr_dataframes = []
    arr_dataframes_variants = []

    # ['tumor', 'relapse']
    arr_samples = np.array(args.samples)
    # loop for samples in the results_mutations.tsv file
    for idx_sample, row_sample in np.ndenumerate(arr_samples):
        # filter by sample
        df_sample = pd.DataFrame(results_mutations_tsv.loc[results_mutations_tsv['sample_id'] == row_sample])

        sample_id = row_sample

        # array of unique clones of the sample
        arr_sample_clones = np.array(df_sample['clone_id'].unique())

        clones_genes_vaf = {}
        genes_vaf = {}

        clones_variants_vaf = {}
        variants_vaf = {}

        # loop for sample's clones in the results_mutations.tsv file
        for idx_smp_clone, row_smp_clone in np.ndenumerate(arr_sample_clones):
            # i.e: /home/cleon/vep_data/933124/tumor/ensemble_tumor.0.tsv
            cluster_id = row_smp_clone

            genes_vaf = {}
            num_drugs = 0

            variants_vaf = {}
            num_variants = 0

            # i.e: /home/cleon/PycharmProjects/MSc-Project/out/pandrugs/933124/tumor/0
            path_computation_tsv = pandrugs_dir + '/' + str(sample_id) + '/' + str(cluster_id) + '/computation.tsv'
            computation_tsv = pd.read_csv(path_computation_tsv, sep='\t')
            url_computation = ''
            computation_id = ''
            for idx_comp_tsv, row_comp_tsv in computation_tsv.iterrows():
                url_computation = row_comp_tsv['url_computation']
                url_computation_words = url_computation.split('/')
                computation_id = url_computation_words[len(url_computation_words) - 1]

            # loop for the gene-drug.ann.tsv
            path_gd_ann_tsv = heatmap_dir + '/' + str(sample_id) + '/' + str(cluster_id) + '/gene-drug.ann.tsv'
            if os.path.exists(path_gd_ann_tsv):
                gd_ann_tsv = pd.read_csv(path_gd_ann_tsv, sep='\t')

                for idx_gd_ann_tsv, row_gd_ann_tsv in gd_ann_tsv.iterrows():
                    drug = row_gd_ann_tsv['drug']
                    gene_vaf = row_gd_ann_tsv['VAF']
                    genes_vaf[drug] = gene_vaf
                    num_drugs = num_drugs + 1

                clones_genes_vaf[cluster_id] = genes_vaf
            else:
                clones_genes_vaf[cluster_id] = {}

            # loop for the xxx-vscore.ann.tsv
            path_var_ann_tsv = heatmap_dir + '/' + str(sample_id) + '/' + str(cluster_id) + '/' + str(computation_id) \
                               + '-vscore.ann.tsv'
            if os.path.exists(path_var_ann_tsv):
                var_ann_tsv = pd.read_csv(path_var_ann_tsv, sep='\t')

                for idx_var_ann_tsv, row_var_ann_tsv in var_ann_tsv.iterrows():
                    variant_id = row_var_ann_tsv['ID']
                    variant_vaf = row_var_ann_tsv['VAF']
                    variants_vaf[variant_id] = variant_vaf
                    num_variants = num_variants + 1

                clones_variants_vaf[cluster_id] = variants_vaf
            else:
                clones_variants_vaf[cluster_id] = {}

        # end of loop for sample's clones in the results_mutations.tsv file
        # print(clones_genes_vaf)
        # print(clones_variants_vaf)

        # heatmap for gene-drug/clones association
        df_heatmap = pd.DataFrame(clones_genes_vaf)
        # print(df_heatmap)
        df_heatmap_norm = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in df_heatmap.items()]))
        df_heatmap_norm = df_heatmap_norm.fillna(0)
        # print(df_heatmap_norm)
        arr_dataframes.append(df_heatmap_norm.copy(True))

        # plot a heatmap with annotation
        sns.set(font_scale=0.5)
        p = sns.heatmap(df_heatmap_norm, annot=True, annot_kws={"size": 7}, cmap="BuGn",
                        cbar_kws={'label': 'VAF'})
        p.set_xlabel("Clones")
        p.set_ylabel("Drugs")
        p.set_title(str(sample_id).upper())

        # Saving the Seaborn Figure:
        heatmap_plots_dir = heatmap_dir + '/' + str(sample_id) + '/plots'
        if not os.path.exists(heatmap_plots_dir):
            os.makedirs(heatmap_plots_dir)

        path_heatmap_png = heatmap_plots_dir + '/heatmap_{}.png'.format(sample_id)
        plt.savefig(path_heatmap_png)

        # heatmap for variants/clones association
        df_heatmap_variants = pd.DataFrame(clones_variants_vaf)
        # print(df_heatmap_variants)
        df_heatmap_variants_norm = pd.DataFrame(
            dict([(k, pd.Series(v)) for k, v in df_heatmap_variants.items()]))
        df_heatmap_variants_norm = df_heatmap_variants_norm.fillna(0)
        # print(df_heatmap_variants_norm)
        arr_dataframes_variants.append(df_heatmap_variants_norm.copy(True))

        # plot a heatmap with annotation
        sns.set(font_scale=0.5)
        p_variants = sns.heatmap(df_heatmap_variants_norm, annot=True, annot_kws={"size": 7}, cmap="BuGn",
                                 cbar_kws={'label': 'VAF'})
        p_variants.set_xlabel("Clones")
        p_variants.set_ylabel("Variant")
        p_variants.set_title(str(sample_id).upper())

        # Saving the Seaborn Figure:
        heatmap_variants_plots_dir = heatmap_dir + '/' + str(sample_id) + '/plots'
        if not os.path.exists(heatmap_variants_plots_dir):
            os.makedirs(heatmap_variants_plots_dir)

        path_heatmap_variants_png = heatmap_variants_plots_dir + '/heatmap_variants_{}.png'.format(sample_id)
        plt.savefig(path_heatmap_variants_png)

    # end of loop for samples in the results_mutations.tsv file
    # TODO: currently for only two samples (in working progress)
    # Exporting the Seaborn Figure for gene-drug/clones association:
    fig, (ax, ax2) = plt.subplots(ncols=2, figsize=(10, 6))
    fig.subplots_adjust(wspace=0.01)
    sns.heatmap(arr_dataframes[0], annot=True, annot_kws={"size": 7}, cmap="BuGn", ax=ax, cbar=False)
    # fig.colorbar(ax.collections[0], ax=ax, location="left", use_gridspec=False, pad=0.2)
    sns.heatmap(arr_dataframes[1], annot=True, annot_kws={"size": 7}, cmap="BuGn", ax=ax2, cbar=True, cbar_kws={'label': 'VAF'}, yticklabels=False)
    # fig.colorbar(ax2.collections[0], ax=ax2, location="right", use_gridspec=False, pad=0.2)
    # ax2.yaxis.tick_right()
    # ax2.tick_params(rotation=0)
    # ax2.set_ylim(0, 1)
    ax.set(title='Tumor')
    ax.set_ylabel("Drugs")
    ax.set_xlabel("Clones")
    ax2.set(title='Relapse')
    ax2.set_xlabel("Clones")

    # plt.show()
    plt.savefig(out_file_gd)

    # Exporting the Seaborn Figure for variant/clones association:
    fig_variants, (ax_variants, ax2_variants) = plt.subplots(ncols=2, figsize=(10, 10))
    fig_variants.subplots_adjust(wspace=0.01)
    sns.heatmap(arr_dataframes_variants[0], annot=True, annot_kws={"size": 7}, cmap="BuGn", ax=ax_variants, cbar=False)
    # fig.colorbar(ax.collections[0], ax=ax, location="left", use_gridspec=False, pad=0.2)
    sns.heatmap(arr_dataframes_variants[1], annot=True, annot_kws={"size": 7}, cmap="BuGn", ax=ax2_variants, cbar=True,
                cbar_kws={'label': 'VAF'}, yticklabels=False)
    # fig.colorbar(ax2.collections[0], ax=ax2, location="right", use_gridspec=False, pad=0.2)
    # ax2.yaxis.tick_right()
    # ax2.tick_params(rotation=0)
    # ax2.set_ylim(0, 1)
    ax_variants.set(title='Tumor')
    ax_variants.set_ylabel("Variants")
    ax_variants.set_xlabel("Clones")
    ax2_variants.set(title='Relapse')
    ax2_variants.set_xlabel("Clones")

    # plt.show()
    plt.savefig(out_file_vars)
