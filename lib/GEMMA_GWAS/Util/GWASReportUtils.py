import os
import pandas as pd
import shutil
import uuid
import logging
import hashlib
from pprint import pprint as pp

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
        self.assembly_info = ''
        self.map = {}

    def find_global_contig_length(self, row):
        list_contigs = self.assembly_info['contigs'].keys()
        list_contigs = [x for x in list_contigs if x.lower().startswith('chr')]

        fail_safe_contig_baselength = self.assembly_info['contigs'][list_contigs[-1]]

        chr = str(row['CHR'])
        assembly_id = self.get_real_contig_id({'chr': chr})

        try:
            length = int(row['POS']) + int(self.contig_baselengths[assembly_id])
        except KeyError:
            length = fail_safe_contig_baselength

        return length


    def filter_gemma_results(self, assoc_info, params):
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
        logging.info('Filering and scrubbing GEMMA results')
        if os.path.exists(params['variation']):
            var_obj = self.dfu.get_objects({'object_refs': ['26322/20/1']})['data'][0]
        else:
            var_obj = self.dfu.get_objects({'object_refs': [params['variation']]})['data'][0]

        genome_ref = var_obj['data']['genome_ref']

        self.state = assoc_info

        for pheno in assoc_info:
            # read gemma results
            gemma_results = pd.read_csv(assoc_info[pheno]['gemma']['file'], sep='\t')
            # sort by lowest p-value
            gemma_results = gemma_results.sort_values(by='p_wald')

            # relabel contig to match assembly
            gemma_results['chr'] = gemma_results.apply(self.get_real_contig_id, axis=1)
            relabeled_gemma_results = os.path.join(self.scratch, 'output',
                                                   'relabeled_' + os.path.basename(self.state[pheno]['gemma']['file']))
            gemma_results.to_csv(path_or_buf=relabeled_gemma_results, sep='\t', index=False)
            self.state[pheno]['gemma']['file'] = relabeled_gemma_results
            self.state[pheno]['gemma']['md5'] = hashlib.md5(open(relabeled_gemma_results, 'rb').read()).hexdigest()

            # filter gemma_results, save to filtered file
            filtered_gemma_file = os.path.join(self.htmldir, pheno + '_filtered_results.txt')
            gemma_results = gemma_results[['chr', 'rs', 'ps', 'p_wald']]
            gemma_results.columns = ['CHR', 'SNP', 'POS', 'P']
            gemma_results['BP'] = gemma_results.apply(self.find_global_contig_length, axis=1)
            gemma_results[:5000].to_csv(path_or_buf=filtered_gemma_file,
                                        sep='\t', index=False)
            self.state[pheno]['gemma']['filtered_file'] = filtered_gemma_file
            self.state[pheno]['gemma']['filtered_md5'] = hashlib.md5(open(filtered_gemma_file, 'rb').read()).hexdigest()
            logging.info(f'Filering and scrubbing GEMMA results done for phenotype {pheno}')

            logging.info(f'Annotating GEMMA results for phenotype {pheno}')
            annotated_file = self.snp2gene.annotate_gwas_results({
                'genome_obj': genome_ref,
                'gwas_result_file': filtered_gemma_file
            })
            logging.info(f'Finished GEMMA result annotation for phenotype {pheno}')

            if os.path.exists(annotated_file['snp_to_gene_list']):
                self.state[pheno]['gemma']['annotated_file'] = annotated_file['snp_to_gene_list']
                self.state[pheno]['gemma']['annotated_md5'] = \
                    hashlib.md5(open(annotated_file['snp_to_gene_list'], 'rb').read()).hexdigest()
            else:
                raise FileNotFoundError(f'Annotated snp2gene results not found at: {annotated_file}')

        return True

    def mk_mapping_file(self, contigs):
        logging.info('Creating mapping files for plink -> assembly ids')
        # map assembly contigs to plink integers
        map_file = "assembly\tplink\n"
        for contig in contigs:
            plink = int(''.join(filter(str.isdigit, contig)))
            map_file += contig + '\t' + str(plink) + '\n'

        map_file_out = os.path.join(self.scratch, 'assembly_to_plink_map.txt')

        with open(map_file_out, 'w') as map:
            map.write(map_file)
            map.close()

        if os.path.exists(map_file_out):
            self.map['file'] = map_file_out
            self.map['md5'] = hashlib.md5(open(map_file_out, 'rb').read()).hexdigest()
            self.map['shock'] = self.dfu.file_to_shock({ 'file_path': map_file_out })
        else:
            raise FileNotFoundError(f'Map file not created: {map_file_out}')

        logging.info('Mapping file created')

        return self.map['file']

    def get_assembly_info(self, params):
        logging.info('Pulling assembly data on contigs')
        # if testing, variation is a flat file
        # so let's hardcode the reference for the variation object
        if os.path.exists(params['variation']):
            var_obj = self.dfu.get_objects({'object_refs': ['26322/20/1']})['data'][0]
        else:
            var_obj = self.dfu.get_objects({'object_refs': [params['variation']]})['data'][0]

        assembly_ref = var_obj['data']['assemby_ref']
        genome_ref = var_obj['data']['genome_ref']

        assembly_obj = self.dfu.get_objects({'object_refs': [assembly_ref]})['data'][0]

        contigs = assembly_obj['data']['contigs']

        contig_ids = list(contigs.keys())
        contig_ids.sort()

        list_contigs = contigs.keys()
        list_contigs = [x for x in list_contigs if x.lower().startswith('chr')]
        # list_contigs = [x.replace('r0', 'r') for x in list_contigs]

        self.mk_mapping_file(list_contigs)

        contig_baselengths = {}
        prev_len = 0

        for id in contig_ids:
            contig_baselengths[id] = prev_len
            prev_len += contigs[id]['length']

        self.contig_baselengths = contig_baselengths

        assembly_info = {
            'ref': assembly_ref,
            'genome_ref': genome_ref,
            'contigs': contig_baselengths
        }

        logging.info('Assembly data collected')

        return assembly_info

    def get_real_contig_id(self, row):
        chr = str(row['chr'])
        map = {}
    
        with open(self.map['file'], 'r') as f:
            next(f)  # ignore tsv header row
            for line in f:
                sp_line = line.strip('\n').split('\t')
                map[sp_line[1]] = sp_line[0]
        try:
            return str(map[chr])
        except KeyError:
            return str(chr)

    def save_assoc_obj(self, params):
        logging.info('Creating association KBase object')
        if os.path.exists(params['variation']):
            var_ref = '26322/20/1'
        else:
            var_ref = params['variation']

        trait_ref = params['trait_matrix']

        assoc_details = []

        for pheno in self.state:
            df_filtered = pd.read_csv(self.state[pheno]['gemma']['file'], sep='\t')
            df_filtered = df_filtered[['chr', 'rs', 'ps', 'p_wald']]
            df_filtered.sort_values(by='p_wald')
            if len(df_filtered.index) > 5000:
                df_filtered = df_filtered[:5000]

            df_filtered['contig_id'] = df_filtered.apply(self.get_real_contig_id, axis=1)
            df_filtered['pve'] = 0.0
            df_filtered = df_filtered[['contig_id', 'rs', 'ps', 'p_wald', 'pve']]
            df_filtered.columns = ['contig_id', 'variant_id', 'position', 'pvalue', 'pve']

            results = []

            for i, row in df_filtered.iterrows():
                if row['variant_id'] == '.':
                    var_id = str(row['contig_id']) + '_' + str(row['position'])
                else:
                    var_id = str(row['variant_id'])

                results.append((str(row['contig_id']), var_id, int(row['position']),
                                float(row['pvalue']), float(row['pve'])))

            a_details = {
                'traits': f'Association analysis of trait values of phenotype: {pheno}',
                'association_results': results
            }

            assoc_details.append(a_details)

        associations_obj = {
            'description' : 'GEMMA Linear mixed models Univariate association analysis on phenotypes',
            'trait_ref': trait_ref,
            'variation_id': var_ref,
            'association_details': assoc_details
        }

        if 'assoc_obj_name' in params:
            assoc_obj_name = params['assoc_obj_name']
        else:
            assoc_obj_name = 'assoc_obj_' + str(uuid.uuid4())

        assoc_obj_saved = self.dfu.save_objects({
            'id': self.dfu.ws_name_to_id(params['workspace_name']),
            'objects': [{
                'type': 'KBaseGwasData.Associations',
                'data': associations_obj,
                'name': assoc_obj_name
            }]
        })[0]

        assoc_obj_ref = str(assoc_obj_saved[6]) + "/" + str(assoc_obj_saved[0]) + "/" + str(assoc_obj_saved[4])
        logging.info(f'Association object {assoc_obj_ref} created.')

        return assoc_obj_ref

    def create_report_msg(self, plink_info, vcf_info):
        reportmsg = "Internal files associated with this pipeline and their md5 values:\n"
        if vcf_info is not None:
            reportmsg += f"\tVCF: {os.path.basename(vcf_info['file'])}\n"
            reportmsg += f"\t\tmd5: {vcf_info['md5']}\n"

        if plink_info is not None:
            reportmsg += "\tPlink bed: " + os.path.basename(plink_info['plink_bed']['file']) + "\n"
            reportmsg += f"\t\tmd5: {plink_info['plink_bed']['md5']}\n"
            reportmsg += "\tPlink bim: " + os.path.basename(plink_info['plink_bim']['file']) + "\n"
            reportmsg += f"\t\tmd5: {plink_info['plink_bim']['md5']}\n"
        reportmsg += "\nSome statistics, separated by phenotypes run, about the GEMMA analysis(es):\n"

        for pheno in self.state:
            reportmsg += f"\tphenotype : {pheno}"
            reportmsg += f"\t\tgemma statistics: {self.state[pheno]['gemma']['stats']['individuals']} individuals,  "
            reportmsg += f"{self.state[pheno]['gemma']['stats']['pve']} pve, {self.state[pheno]['gemma']['stats']['se_pve']} se_pve\n"
            reportmsg += "\t\tgemma results file:" + os.path.basename(self.state[pheno]['gemma']['file']) + "\n"
            reportmsg += "\t\tgemma results md5:" + self.state[pheno]['gemma']['md5'] + "\n"
            reportmsg += "\t\tkinship matrix file:" + os.path.basename(self.state[pheno]['kinship']['file']) + "\n"
            reportmsg += "\t\tkinship matrix md5:" + self.state[pheno]['kinship']['md5'] + "\n"
            reportmsg += "\t\tphenotype value file:" + os.path.basename(self.state[pheno]['fam']['file']) + "\n"
            reportmsg += "\t\tphenotype value md5:" + self.state[pheno]['fam']['md5'] + "\n"

        return reportmsg

    def create_file_links(self):
        logging.info('Creating File output links.')
        file_links = []

        for pheno in self.state:
            file_links.append({
                'path': self.state[pheno]['gemma']['file'],
                'name': 'Raw_GEMMA_results.txt'
            })
            file_links.append({
                'path': self.state[pheno]['gemma']['annotated_file'],
                'name': 'Filtered_GEMMA_results.txt'
            })
            #file_links.append({
            #    'path': self.state[pheno]['kinship']['file'],
            #    'name': f'Kinship matrix file generated for phenotype: {pheno}'
            #})
            file_links.append({
                 'path': self.state[pheno]['fam']['file'],
                 'name': f'{pheno}.fam'
            })

        logging.info('File output links created.')

        return file_links

    def create_html_outputs(self):
        logging.info('Creating HTML outputs.')
        phenojs = 'var inputs = ['
        i = 0

        for pheno in self.state:
            if i == 0:
                phenojs += '"' + os.path.basename(self.state[pheno]['gemma']['annotated_file']) + '"'
            else:
                phenojs += ',"' + os.path.basename(self.state[pheno]['gemma']['annotated_file']) + '"'
            i += 1

        phenojs += '];\n'

        phenojs += 'var contig_lengths = {\n'

        list_contigs = self.assembly_info['contigs'].keys()
        list_contigs = [x for x in list_contigs if x.lower().startswith('chr')]

        for contig in list_contigs:
            phenojs += '\t'+str(contig)+': '+str(self.contig_baselengths[contig])+',\n'

        phenojs = phenojs[:-2] + '\n};'

        with open(os.path.join(self.htmldir, 'pheno.js'), 'w') as jsf:
            jsf.write(phenojs)
            jsf.close()

        html_dir = [{
            'path': self.htmldir,
            'name': "index.html",
            'description': 'Manhattan plot of GEMMA GWAS association tests'
        }]
        logging.info("Finished creating HTML outputs")

        return html_dir

    def mk_output(self, params, assoc_info):
        plink_info = assoc_info.pop('plink', None)

        vcf_info = assoc_info.pop('vcf', None)

        self.assembly_info = self.get_assembly_info(params)
        self.filter_gemma_results(assoc_info, params)
        obj_ref = self.save_assoc_obj(params)

        logging.info("Report object creation started")

        reportobj = {
            'message': self.create_report_msg(plink_info, vcf_info),
            'objects_created': [{'ref': obj_ref, 'description': 'Association GWAS object from GEMMA algorithm.'}],
            'direct_html': None,
            'direct_html_link_index': 0,
            'html_links': self.create_html_outputs(),
            'file_links': self.create_file_links(),
            'report_object_name': 'GEMMA_GWAS_report_' + str(uuid.uuid4()),
            'workspace_name': params['workspace_name']
        }

        logging.info("Report object created.")

        pp(self.state)

        return reportobj
