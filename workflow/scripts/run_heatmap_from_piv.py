import seaborn as sns
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import json
from matplotlib.patches import Rectangle

if __name__ == '__main__':

    # get the script input params
    input_parser = argparse.ArgumentParser()
    input_parser.add_argument("--prevalence", action='store')
    input_parser.add_argument("--patient", action='store')
    input_parser.add_argument("--samples", nargs='+', action='append')
    input_parser.add_argument("--pandrugs_dir", action='store')
    input_parser.add_argument("--heatmap_dir", action='store')
    input_parser.add_argument("--out_file_gd", action='store')
    input_parser.add_argument("--out_file_vars", action='store')
    args = input_parser.parse_args()

    # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/pandrugs'
    pandrugs_dir = args.pandrugs_dir

    # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/heatmap'
    heatmap_dir = args.heatmap_dir

    # i.e.: /home/cleon/PycharmProjects/MSc-Project/out/pyclone-vi/UPN933124/output_933124.tsv
    prevalence_tsv = pd.read_csv(args.prevalence, sep='\t')

    # i.e.: 'UPN933124'
    patient_id = args.patient

    # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/heatmap/UPN933124/heatmap_all_compared.png'
    out_file_gd = args.out_file_gd

    # i.e: '/home/cleon/PycharmProjects/MSc-Project/out/heatmap/UPN933124/heatmap_all_vars_compared.png'
    out_file_vars = args.out_file_vars

    # ['tumor', 'relapse']
    arr_samples = np.array(args.samples)

    arr_dataframes = []
    arr_dataframes_variants = []

    # loop for samples in the prevalence_tsv file
    for idx_sample, row_sample in np.ndenumerate(arr_samples):
        # print(f'row_sample='.format(row_sample))

        # filter by sample
        df_sample = pd.DataFrame(prevalence_tsv.loc[prevalence_tsv['sample_id'] == row_sample])

        # array of unique clones of the sample
        arr_sample_clones = np.array(df_sample['cluster_id'].unique())

        clones_genes_vaf = {}
        genes_vaf = {}

        clones_variants_vaf = {}
        variants_vaf = {}

        # loop for sample's clones in the prevalence_tsv file
        for idx_smp_clone, row_smp_clone in np.ndenumerate(arr_sample_clones):
            # i.e: /home/cleon/vep_data/UPN933124/tumor/ensemble_tumor.0.tsv
            sample_id = row_sample
            cluster_id = row_smp_clone

            genes_vaf = {}
            num_drugs = 0

            variants_vaf = {}
            num_variants = 0

            # i.e: /home/cleon/PycharmProjects/MSc-Project/out/pandrugs/UPN933124/tumor/0
            path_computation_tsv = pandrugs_dir + '/' + str(patient_id) + '/' + str(sample_id) + '/' + str(cluster_id) + '/computation.tsv'
            if os.path.exists(path_computation_tsv):
                computation_tsv = pd.read_csv(path_computation_tsv, sep='\t')
                url_computation = ''
                computation_id = ''
                for idx_comp_tsv, row_comp_tsv in computation_tsv.iterrows():
                    url_computation = row_comp_tsv['url_computation']
                    url_computation_words = url_computation.split('/')
                    computation_id = url_computation_words[len(url_computation_words) - 1]

                # loop for the gene-drug.ann.tsv
                path_gd_ann_tsv = heatmap_dir + '/' + str(patient_id) + '/' + str(sample_id) + '/' + str(cluster_id) + '/gene-drug.ann.tsv'
                if os.path.exists(path_gd_ann_tsv):
                    gd_ann_tsv = pd.read_csv(path_gd_ann_tsv, sep='\t')

                    for idx_gd_ann_tsv, row_gd_ann_tsv in gd_ann_tsv.iterrows():
                        gd_gene_obj_json_str = str(row_gd_ann_tsv['gene']).replace('[', '').replace(']', '').replace("'", "\"")
                        gd_gene_json_object = json.loads(gd_gene_obj_json_str)
                        gd_gene_symbol = gd_gene_json_object["geneSymbol"]

                        str_drug = row_gd_ann_tsv['drug']
                        str_status = row_gd_ann_tsv['status']

                        if str_status == "APPROVED":
                            drug = f'{gd_gene_symbol}  -  {str_drug} [{str_status}]'
                        else:
                            drug = f'{gd_gene_symbol}  -  {str_drug}'

                        gene_vaf = row_gd_ann_tsv['VAF']
                        genes_vaf[drug] = gene_vaf
                        num_drugs = num_drugs + 1

                    clones_genes_vaf[cluster_id] = genes_vaf
                else:
                    clones_genes_vaf[cluster_id] = {}

                # loop for the xxx-vscore.ann.tsv
                path_var_ann_tsv = heatmap_dir + '/' + str(patient_id) + '/' + str(sample_id) + '/' + str(cluster_id) + '/' + str(computation_id) \
                                   + '-vscore.ann.tsv'
                if os.path.exists(path_var_ann_tsv):
                    var_ann_tsv = pd.read_csv(path_var_ann_tsv, sep='\t')

                    for idx_var_ann_tsv, row_var_ann_tsv in var_ann_tsv.iterrows():
                        # variant_id = row_var_ann_tsv['gene_hgnc']
                        variant_id = row_var_ann_tsv['ID']
                        variant_vaf = row_var_ann_tsv['VAF']
                        variants_vaf[variant_id] = variant_vaf
                        num_variants = num_variants + 1

                    clones_variants_vaf[cluster_id] = variants_vaf
                else:
                    clones_variants_vaf[cluster_id] = {}
        # print(clones_genes_vaf)
        # print(clones_variants_vaf)

        # heatmap for gene-drug/clones association
        df_heatmap = pd.DataFrame(clones_genes_vaf)
        if df_heatmap.size > 0:
            # print(df_heatmap)
            df_heatmap_norm = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in df_heatmap.items()]))
            df_heatmap_norm = df_heatmap_norm.fillna(0)
            # print(df_heatmap_norm)
            arr_dataframes.append(df_heatmap_norm.copy(True))

            # plot a heatmap with annotation
            # sns.set(font_scale=0.2)
            msk_hd = arr_dataframes[0] == 0
            # sns.set(rc={'figure.figsize': (14, 15)})
            with sns.axes_style("white"):
                f, p_ax = plt.subplots(figsize=(25, 15))
                p_ax = sns.heatmap(df_heatmap_norm, linewidth=0.1, linecolor='gray', annot=True, annot_kws={"size": 12}, cmap="Reds", mask=msk_hd, cbar_kws={'label': 'VAF'})
                p_ax.set_xlabel("Clones")
                p_ax.set_ylabel("Drugs")
                p_ax.set_title(str(sample_id).upper())
                # p_ax.add_patch(Rectangle((2, 48), -3, 2, fill=False, edgecolor='green', lw=1))

            # Saving the Seaborn Figure:
            heatmap_plots_dir = heatmap_dir + '/' + str(patient_id) + '/' + str(sample_id) + '/plots'
            if not os.path.exists(heatmap_plots_dir):
                os.makedirs(heatmap_plots_dir)

            path_heatmap_png = heatmap_plots_dir + '/heatmap_{}.png'.format(sample_id)
            plt.savefig(path_heatmap_png)

            # heatmap for variants/clones association
            df_heatmap_variants = pd.DataFrame(clones_variants_vaf)
            # print(df_heatmap_variants)
            df_heatmap_variants_norm = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in df_heatmap_variants.items()]))
            df_heatmap_variants_norm = df_heatmap_variants_norm.fillna(0)
            # print(df_heatmap_variants_norm)
            arr_dataframes_variants.append(df_heatmap_variants_norm.copy(True))

            # plot a heatmap with annotation
            # sns.set(font_scale=0.2)
            sns.set(rc={'figure.figsize': (15, 8)})
            p_variants = sns.heatmap(df_heatmap_variants_norm, annot=True, annot_kws={"size": 7}, cmap="Reds", cbar_kws={'label': 'VAF'})
            p_variants.set_xlabel("Clones")
            p_variants.set_ylabel("Genes")
            p_variants.set_title(str(sample_id).upper())

            # Saving the Seaborn Figure:
            heatmap_variants_plots_dir = heatmap_dir + '/' + str(patient_id) + '/' + str(sample_id) + '/plots'
            if not os.path.exists(heatmap_variants_plots_dir):
                os.makedirs(heatmap_variants_plots_dir)

            path_heatmap_variants_png = heatmap_variants_plots_dir + '/heatmap_genes_{}.png'.format(sample_id)
            plt.savefig(path_heatmap_variants_png)

    # TODO: currently for only two samples (in working progress)
    # Exporting the Seaborn Figure for gene-drug/clones association:
    sns.set(font_scale=1)
    fig, (ax, ax2) = plt.subplots(ncols=2, figsize=(30, 15))
    fig.subplots_adjust(wspace=0.01)
    # sns.set_theme(context='paper', style='ticks', font_scale=1.5)

    if np.size(arr_dataframes) == 2:
        msk1 = arr_dataframes[0] == 0
        msk2 = arr_dataframes[1] == 0

        sns.heatmap(arr_dataframes[0], annot=True, linewidth=0.1, annot_kws={"fontsize": "medium"}, mask=msk1, cmap="Reds", ax=ax, cbar=False)
        # fig.colorbar(ax.collections[0], ax=ax, location="left", use_gridspec=False, pad=0.2)
        sns.heatmap(arr_dataframes[1], annot=True, linewidth=0.1, annot_kws={"fontsize": "medium"}, mask=msk2, cmap="Reds", ax=ax2, cbar=True, cbar_kws={'label': 'VAF'}, yticklabels=False)
        # fig.colorbar(ax2.collections[0], ax=ax2, location="right", use_gridspec=False, pad=0.2)
        # ax2.yaxis.tick_right()
        # ax2.tick_params(rotation=0)
        # ax2.set_ylim(0, 1)
        ax.set(title='Tumor')
        ax.set_ylabel("Drugs")
        ax.set_xlabel("Clones")
        ax2.set(title='Relapse')
        ax2.set_xlabel("Clones")

        # highlighting approved drugs
        for text in ax.texts:
            # text.set_size(14)
            if str(text.get_text()).__contains__('APPROVED'):
                # text.set_size(18)
                text.set_weight('bold')
                # text.set_style('italic')

        # plt.show()
        plt.savefig(out_file_gd)

        sns.set(font_scale=1.2)

        # mask = np.zeros_like(arr_dataframes_variants[0])
        mask1 = arr_dataframes_variants[0] == 0
        mask2 = arr_dataframes_variants[1] == 0

        # Exporting the Seaborn Figure for gene/clones association:
        with sns.axes_style("white"):
            fig_variants, (ax_variants, ax2_variants) = plt.subplots(ncols=2, figsize=(20, 15))
            fig_variants.subplots_adjust(wspace=0.01)
            pv_1 = sns.heatmap(arr_dataframes_variants[0], annot=False, annot_kws={"size": 6}, cmap="Reds", mask=mask1, ax=ax_variants, cbar=False, yticklabels=False)
            # pv_1.invert_yaxis()
            # make frame visible
            for _, spine in pv_1.spines.items():
                spine.set_visible(True)
            # fig.colorbar(ax.collections[0], ax=ax, location="left", use_gridspec=False, pad=0.2)
            pv_2 = sns.heatmap(arr_dataframes_variants[1], annot=False, annot_kws={"size": 6}, cmap="Reds", mask=mask2, ax=ax2_variants, cbar=True,
                               cbar_kws={'label': 'VAF'}, yticklabels=False)
            # pv_2.invert_yaxis()
            # make frame visible
            for _, spine in pv_2.spines.items():
                spine.set_visible(True)

        # fig.colorbar(ax2.collections[0], ax=ax2, location="right", use_gridspec=False, pad=0.2)
        # ax2.yaxis.tick_right()
        # ax2.tick_params(rotation=0)
        # ax2.set_ylim(0, 1)
        ax_variants.set(title='Tumor')
        ax_variants.set_ylabel("Variants")
        ax_variants.set_xlabel("Clones")

        ax2_variants.set(title='Relapse')
        ax2_variants.set_xlabel("Clones")
        # ax2_variants.add_patch(Rectangle((3, 4), 1, 1, fill=True, edgecolor='blue', lw=1))

        # plt.show()
        plt.savefig(out_file_vars)
