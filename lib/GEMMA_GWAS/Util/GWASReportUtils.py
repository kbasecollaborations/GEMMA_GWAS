import os
import shutil
import logging
import pandas as pd
from pprint import pprint as pp

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.snp2geneClient import snp2gene


class GWASReportUtils:
    def __init__(self, config):
        self.scratch = config["scratch"]
        self.callback_url = config["SDK_CALLBACK_URL"]
        self.htmldir = os.path.join(self.scratch, 'mhplot')

        # modules
        self.dfu = DataFileUtil(self.callback_url)
        self.snp2gene = snp2gene(self.callback_url)

        if os.path.isdir(os.path.join(self.scratch,'mhplot')):
            shutil.rmtree(os.path.join(self.scratch,'mhplot'))
        shutil.copytree('/kb/module/lib/GEMMA_GWAS/Util/Report/mhplot/', os.path.join(self.scratch,'mhplot'))

    def label_dot_snps(self, row):
        if row['rs'] == '.':
            row['rs'] = 'Chr' + str(row['chr']) + '_' + str(row['ps'])

        return row

    def filter_gemma_results(self, assoc_files):
        filtered_files = dict.fromkeys(assoc_files.keys(), '')
        self.filtered_result_dir = os.path.join(self.scratch, 'filtered_gemma')

        if not os.path.isdir(self.filtered_result_dir):
            os.mkdir(self.filtered_result_dir)

        for pheno, assoc_file in assoc_files.items():
            df = pd.read_csv(assoc_file, sep='\t')

            chrs = list(set(list(df['chr'])))

            chrs_no_scaff = [int(x) for x in chrs if 'scaffold' not in str(x)]
            chrs_no_scaff = list(set(chrs_no_scaff))
            chrs_no_scaff.sort()
            scaffold_chr = chrs_no_scaff[-1] + 1

            df['chr'] = df['chr'].apply(lambda x: x if 'scaffold' not in str(x) else scaffold_chr)
            df['p_wald'] = df['p_wald'].apply(lambda x: float(x))
            # TODO: test for snp = '.'
            df = df.apply(self.label_dot_snps, axis=1)

            df2 = df[['rs', 'chr', 'ps', 'p_wald']]
            df2.columns = ['SNP', 'CHR', 'BP', 'P']

            df3 = df2[df2['P'] <= 0.05]
            df3 = df3.sort_values(by='P', ascending=True)

            filtered_file = os.path.join(self.filtered_result_dir, pheno+'_filtered.txt')

            df3.to_csv(filtered_file, sep='\t', index=False)

            filtered_files[pheno] = filtered_file

        return filtered_files

    def mk_html_report(self, assoc_files, stats):
        with open(os.path.join(self.htmldir, 'pheno.js'), 'w') as f:
            f.write("var inputs = {\n")
            phenolen = len(assoc_files)
            i = 0

            for pheno, file in assoc_files.items():
                if i == (phenolen -1):
                    f.write('\t\'' + pheno + '\': \'' + file + '\'\n}\n')
                else:
                    f.write('\t\'' + pheno + '\': \'' + file + '\',\n')
                i+=1

            f.write('var ind = ' + stats['individuals'] + ';\n')
            f.close()

        html_return = {
            'path': self.htmldir,
            'name': "index.html",
            'description': 'Manhattan plot of GEMMA GWAS association tests'
        }

        return html_return

    def mk_assoc_object(self, assoc_results, params):
        var_ref = params['variation']
        trait_ref = params['trait_matrix']

        for pheno, file in assoc_results.items():
            df = pd.read_csv(file, sep='\t')
            df_tuple = list(df.itertuples(index=False, name=None))



    def mk_output(self, params, assoc_results, assoc_stats):
        if params['model'] is 0:
            filtered_files = self.filter_gemma_results(assoc_results)
            html_out = self.mk_html_report(filtered_files, assoc_stats)
            assoc_obj = self.mk_assoc_object(filtered_files, params)
            exit()
        elif params['model'] is 1:
            raise NotImplementedError('Only univariate analysis are supported right now.')
        else:
            raise ValueError('GEMMA linear mixed model is not set.')

        return reportobj
