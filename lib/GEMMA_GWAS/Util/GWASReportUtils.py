import os
import json
import csv
import shutil
import uuid

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.snp2geneClient import snp2gene

class GWASReportUtils:
    def __init__(self, config):
        self.config = config
        self.scratch = config["scratch"]
        self.callback_url = config["SDK_CALLBACK_URL"]
        self.dfu = DataFileUtil(self.callback_url)
        self.snp2gene = snp2gene(self.callback_url)
        if os.path.isdir(os.path.join(self.scratch,'mhplot')):
            shutil.rmtree(os.path.join(self.scratch,'mhplot'))
        shutil.copytree('/kb/module/lib/GEMMA_GWAS/Util/Report/mhplot/', os.path.join(self.scratch,'mhplot'))
        self.htmldir = os.path.join(self.scratch,'mhplot')

    def _filter_assoc_results(self, trait_info, var_ref, model):
        if trait_info['gemma'] == '':
            exit(trait_info['id']+" has an empty gemma element.")

        if not trait_info['gemma'] == 'fail':
            tsv_unfiltered = csv.reader(open(trait_info['gemma'], 'r', newline=''), delimiter='\t')
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

            tsv_filtered_headers = "SNP\tCHR\tBP\tP\tPOS\n"

            if model is 0:
                filtered_tsv_file = os.path.join(self.htmldir, 'snpdata'+str(trait_info['id'])+'.tsv')
            elif model is 1:
                filtered_tsv_file = os.path.join(self.htmldir, 'snpdata-multi.tsv')
            else:
                raise ValueError('GEMMA Linear Mixed Model not set from UI.')

            assoc_entry_limit = 5000
            assoc_details = []

            # correct lengths for scatter plot
            assembly = self.dfu.get_objects({'object_refs': [var_ref]})['data'][0]
            # TODO: change to assembly_ref, change data model
            assembly_ref = assembly['data']['assemby_ref']
            genome_ref = assembly['data']['genome_ref']

            assembly_obj = self.dfu.get_objects({'object_refs': [assembly_ref]})['data'][0]
            contigs = assembly_obj['data']['contigs']
            contig_ids = []

            for contig in contigs:
                if contig.upper().startswith('CHR'):
                    try:
                        int(contig[3:])
                        contig_ids.append(int(contig[3:]))
                    except ValueError:
                        pass
                else:
                    try:
                        int(contig)
                        contig_ids.append(int(contig))
                    except ValueError:
                        pass

            contig_ids.sort()
            contig_baselengths = {}
            prev_len = 0

            from pprint import pprint as pp
            pp(contig_ids)
            exit()

            for id in contig_ids:
                try:
                    contig_baselengths[id] = prev_len
                    prev_len += contigs[str(id)]['length']
                except KeyError:
                    contig_baselengths[id] = prev_len
                    prev_len += contigs['Chr'+str(id)]['length']

            with open(filtered_tsv_file,'w') as tsv_filtered:
                tsv_filtered.write(tsv_filtered_headers)

                k = 0

                if len(tsv_sorted) > assoc_entry_limit:
                    for snp in tsv_sorted:
                        if k < assoc_entry_limit:
                            tsv_filtered.write(snp[1]+"\t"+snp[0]+"\t"+snp[2]+"\t"+snp[13]+"\t" \
                                               + str((contig_baselengths[int(snp[0])]+int(snp[2])))+"\n")
                            k += 1
                        assoc_details.append((snp[1], snp[0], int(snp[2]), float(snp[13]), float(snp[6])))
                else:
                    for snp in tsv_sorted:
                        tsv_filtered.write(snp[1]+"\t"+snp[0]+"\t"+snp[2]+"\t"+snp[13]+"\t" \
                                           + str((contig_baselengths[int(snp[0])]+int(snp[2])))+"\n")
                        assoc_details.append((snp[1], snp[0], int(snp[2]), float(snp[13]), float(snp[6])))

                tsv_filtered.close()

            # annotate gwas result file
            self.snp2gene.annotate_gwas_results({
                'genome_obj': genome_ref,
                'gwas_result_file': filtered_tsv_file
            })

            return assoc_details
        else:
            return False

    def _mk_html_report(self, trait_info, var_ref, model):
        assoc_results = self._filter_assoc_results(trait_info, var_ref, model)

        """ 
        TODO: change output here, actually parse # variants for FDR and Bonferoni
        logging.info("\n\n\nfiltered:\n")
        if model is 0:
            os.system("wc -l "+os.path.join(self.htmldir, 'snpdata' + str(trait_info['id']) + '.tsv'))
        elif model is 1:
            os.system("wc -l " + os.path.join(self.htmldir, 'snpdata-multi.tsv'))
        else:
            raise ValueError('GEMMA Linear Mixed Model not set from UI.')
        logging.info("\n\n\nunfiltered:\n")
        os.system("wc -l " + trait_info['gemma'])
        logging.info("\n\n")
        """

        html_return = {
            'path': self.htmldir,
            'name': "index.html",
            'description': 'Manhattan plot of GEMMA GWAS association tests'
        }

        return html_return, assoc_results

    def _save_assoc_obj(self, params, assoc_results, assoc_details_list):
        assoc_details = []
        if params['model'] is 0:
            for x in range(0, len(assoc_details_list)):
                assoc_details_entry = {
                    'traits': assoc_results[x]['id'],
                    'association_results': assoc_details_list[x]
                }
                assoc_details.append(assoc_details_entry)
        elif params['model'] is 1:
            assoc_details_entry = {
                'traits': assoc_results['multi']['id'],
                'association_results': assoc_details_list[0]
            }
            assoc_details.append(assoc_details_entry)

        assoc = {
            'description': 'test description',
            'variation_id': params['variation'],
            'trait_ref': params['trait_matrix'],
            'association_details': assoc_details
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

    def mk_output(self, params, assoc_results):
        assoc_details = []
        html_info = []
        js_pheno_inputs = []
        file_links = []
        failed_phenos = []

        if params['model'] is 0:
            for x in range(0, len(assoc_results)):
                html_info_entry, assoc_details_entry = self._mk_html_report(assoc_results[x],
                                                                            params['variation'], params['model'])
                if assoc_details_entry:
                    assoc_details.append(assoc_details_entry)
                else:
                    failed_phenos.append(assoc_results[x]['id'])
                html_info.append(html_info_entry)
                js_pheno_inputs.append('snpdata'+str(assoc_results[x]['id'])+'.tsv')
                file_links.append({
                    'path': os.path.join(self.htmldir,'snpdata'+str(assoc_results[x]['id'])+'.tsv'),
                    'name': 'GEMMA Association results for phenotype '+str(assoc_results[x]['id'])+'.'
                })
        elif params['model'] is 1:
            html_info_entry, assoc_details_entry = self._mk_html_report(assoc_results['multi'],
                                                                        params['variation'], params['model'])
            if assoc_details_entry:
                assoc_details.append(assoc_details_entry)
            else:
                failed_phenos.append(assoc_results['multi']['id'])

            html_info.append(html_info_entry)
            js_pheno_inputs.append('snpdata-multi.tsv')
            file_links.append({
                'path': os.path.join(self.htmldir, 'snpdata-multi.tsv'),
                'name': 'GEMMA Association results for phenotypes: ' + str(assoc_results['multi']['id']) + '.'
            })
        else:
            raise ValueError('GEMMA linear mixed model is not set.')

        with open(os.path.join(self.htmldir, 'pheno.js'), 'w') as f:
            f.write("var inputs = ['")
            for x in range(0, len(js_pheno_inputs)):
                if x is (len(js_pheno_inputs)-1):
                    f.write(js_pheno_inputs[x] + "'];\n")
                else:
                    f.write(js_pheno_inputs[x]+"', '")
            f.close()

        assoc_obj = self._save_assoc_obj(params, assoc_results, assoc_details)

        if failed_phenos:
            reportmsg = "The variation object: " + str(params['variation']) + "\n"+ \
                        "The association object: " + str(assoc_obj) + "\n" \
                        "Some phenotypes failed association testing:" + "\n"
            for pheno in failed_phenos:
                reportmsg += pheno
        else:
            reportmsg = "The variation object: " + str(params['variation']) + "\n" + \
                        "The association object: " + str(assoc_obj) + "\n"

        reportobj = {
            'message': reportmsg,
            'objects_created': [{'ref': assoc_obj, 'description': 'Association GWAS object from GEMMA algorithm.'}],
            'direct_html': None,
            'direct_html_link_index': 0,
            'html_links': html_info,
            'file_links': file_links,
            'report_object_name': 'GEMMA_GWAS_report_' + str(uuid.uuid4()),
            'workspace_name': params['workspace_name']
        }

        return reportobj
