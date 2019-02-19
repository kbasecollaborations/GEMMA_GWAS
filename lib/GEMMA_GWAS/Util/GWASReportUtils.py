import os
import json
import csv
import shutil
import logging
import uuid

from installed_clients.DataFileUtilClient import DataFileUtil

class GWASReportUtils:
    def __init__(self, config):
        self.config = config
        self.scratch = config["scratch"]
        self.callback_url = config["SDK_CALLBACK_URL"]
        self.dfu = DataFileUtil(self.callback_url)
        shutil.copytree('/kb/module/lib/GEMMA_GWAS/Util/Report/mhplot/', os.path.join(self.scratch,'mhplot'))
        self.htmldir = os.path.join(self.scratch,'mhplot')

    def _filter_assoc_results(self, assoc_file):
        tsv_unfiltered = csv.reader(open(assoc_file, 'r', newline=''), delimiter='\t')
        # skip old csv headers
        next(tsv_unfiltered, None)

        # gemma assoc results columns:
        # 0 - chr - chromosome numbers
        # 1 - rs - snp ids
        # 2 - ps - base pair positions on the chromosome
        # 3 - n_miss - number of missing individuals for a given snp
        # 4 - allele1 - minor allele
        # 5 - allele0 - major allele
        # 6 - af - allele frequency
        # 7 - beta - beta estimates
        # 8 - se - standard errors for beta
        # 9 - logl_H1 -
        # 10 - l_remle - remle estimates for lambda
        # 11 - l_mle -
        # 12 - p_wald - p value from the Wald test
        # 13 - p_lrt - likelihood ratio test
        # 14 - p_score - p score test

        tsv_sorted = sorted(tsv_unfiltered, key=lambda col: float(col[12]))

        tsv_filtered_headers = "SNP\tCHR\tBP\tP\tAF\n"
        filtered_tsv_file = os.path.join(self.htmldir, 'snpdata.tsv')
        assoc_entry_limit = 5000
        assoc_details = []

        with open(filtered_tsv_file,'w') as tsv_filtered:
            tsv_filtered.write(tsv_filtered_headers)

            k = 0

            if len(tsv_sorted) > assoc_entry_limit:
                for snp in tsv_sorted:
                    if k < assoc_entry_limit:
                        tsv_filtered.write(snp[1]+"\t"+snp[0]+"\t"+snp[2]+"\t"+snp[13]+"\n")
                        k += 1
                    assoc_details.append((snp[1], snp[0], int(snp[2]), float(snp[13]), float(0.0)))
            else:
                for snp in tsv_sorted:
                    tsv_filtered.write(snp[1]+"\t"+snp[0]+"\t"+snp[2]+"\t"+snp[13]+"\n")
                    assoc_details.append((snp[1], snp[0], int(snp[2]), float(snp[13]), float(0.0)))

            tsv_filtered.close()

        return assoc_details

    def _mk_html_report(self, assoc_file):
        assoc_results = self._filter_assoc_results(assoc_file)

        logging.info("\n\n\nfiltered:\n")
        os.system("wc -l "+os.path.join(self.htmldir, 'snpdata.tsv'))
        logging.info("\n\n\nunfiltered:\n")
        os.system("wc -l " + assoc_file)
        logging.info("\n\n")

        html_return = {
            'path': self.htmldir,
            'name': "index.html",
            'description': 'Manhattan plot of GEMMA GWAS association tests'
        }

        return html_return, assoc_results

    def _save_assoc_obj(self, params, assoc_details_list):
        assoc_details = {
            'traits': "trait description",
            'association_results': assoc_details_list
        }

        assoc = {
            'description': 'test description',
            'variation_id': params['variation'],
            'trait_ref': params['trait_matrix'],
            'association_details': [assoc_details]
        }

        if 'assoc_obj_name' in params:
            assoc_obj_name = params['assoc_obj_name']
        else:
            assoc_obj_name = 'assoc_obj_'+str(uuid.uuid4())

        assoc_obj = self.dfu.save_objects({
            'id': self.dfu.ws_name_to_id(params['workspace_name']),
            'objects': [{
                'type': 'KBaseGwasData.Associations',
                'data': assoc,
                'name': assoc_obj_name
            }]
        })[0]

        assoc_obj_ref = str(assoc_obj[6]) + "/" + str(assoc_obj[0]) + "/" + str(assoc_obj[4])

        return assoc_obj_ref

    def mk_output(self, params, assoc_file):
        html_info, assoc_details_list = self._mk_html_report(assoc_file)
        assoc_obj = self._save_assoc_obj(params, assoc_details_list)

        reportobj = {
            'message': "The variation object: " + str(params['variation']) + "\nThe association object:" +
                       str(assoc_obj),
            'objects_created': [{'ref': assoc_obj, 'description': 'Association GWAS object from GEMMA algorithm.'}],
            'direct_html': None,
            'direct_html_link_index': 0,
            'html_links': [html_info],
            'file_links': [],
            'report_object_name': 'GEMMA_GWAS_report_' + str(uuid.uuid4()),
            'workspace_name': params['workspace_name']
        }

        return reportobj
