import subprocess
import os
import operator
import shutil
import csv
import logging

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class AssociationUtils:
    def __init__(self, config, varfiles):
        self.dfu = DataFileUtil(config['SDK_CALLBACK_URL'])
        # TODO: input variable for workspace url
        self.wsc = Workspace("https://appdev.kbase.us/services/ws")
        self.scratch = config["scratch"]
        self.plink_pref = 'plink_vcf'
        self._process_varfiles(varfiles)

    def _process_varfiles(self, files):
        if not isinstance(files, (list,)):
            # get extension from single file
            file_name, file_ext = os.path.splitext(files)
            if file_ext == '.vcf':
                if os.path.exists(files):
                    self.varfile = files
                else:
                    raise IOError('Variation file provided does not exist or is not readable.')
            else:
                raise ValueError('Only supporting VCF variation as input right now.')
        else:
            raise ValueError('Not accepting list of files as variation input yet!')

    def _mk_plink_bin_uni(self):
        plinkvars = ['--make-bed', '--vcf', self.varfile, '--allow-extra-chr', '--out', self.plink_pref]
        plinkcmd = ['plink']

        try:
            proc = subprocess.Popen(plinkcmd, cwd=self.scratch)
            proc.wait()
        except Exception as e:
            exit(e)

        for arg in plinkvars:
            plinkcmd.append(arg)

        plink_bed = os.path.join(self.scratch, plink_base_prefix + str(x) + '.bed')
        plink_bim = os.path.join(self.scratch, plink_base_prefix + str(x) + '.bim')
        plink_fam = os.path.join(self.scratch, plink_base_prefix + str(x) + '.fam')

        if not os.path.exists(plink_bed):
            raise FileNotFoundError('Plink bed doesn\'t exist')

        if not os.path.exists(plink_bim):
            raise FileNotFoundError('Plink bim doesn\'t exist')

        if os.path.exists(plink_fam):
            self.plink_fam_template = shutil.move(plink_fam, os.path.join(self.scratch, 'plink_fam_template.fam'))

        return True

    def _mk_phenos_from_trait_matrix_uni(self, trait_matrix_ref):
        logging.info("Parsing trait matrix.")

        trait_matrix_obj = self.dfu.get_objects({'object_refs': [trait_matrix_ref]})['data'][0]
        ws_type = trait_matrix_obj['info'][2]
        if 'KBaseMatrices.TraitMatrix' in ws_type:
            fam_dir = os.mkdir(os.path.join(self.scratch, 'fams'))
            phenotypes = trait_matrix_obj['data']['data']['row_ids']
            fids = trait_matrix_obj['data']['data']['col_ids']
            vi = 0

            for pheno in phenotypes:
                phenodict = {}
                values = trait_matrix_obj['data']['data']['values'][vi]

                for x in range(0, len(fids)):
                    phenodict[fids[x]] = values[x]

                exit(phenodict)
        else:
            raise ValueError('Cannot write data to VCF; invalid WS type (' + ws_type +
                             ').  Supported types is KBaseMatrices.TraitMatrix')

    def process_gemma_out(self, output):
        stats = {}
        for line in output:
            # ## number of analyzed individuals = XX
            print(line)
            if 'number of analyzed individuals' in line:
                spline = line.split(' ')
                stats['individuals'] = spline[6]
        return stats

    def run_assoc_exp(self, params):
        if params['model'] is 0:
            # univariate analysis
            plink = self._mk_plink_bin_uni(phenotype)
            phenotype = self._mk_phenos_from_trait_matrix_uni(params['trait_matrix'])
            kinmatrix = self._mk_centered_kinship(plink)
            gemma, gemma_output = self.run_gemma_assoc_uni(kinmatrix)
        else:
            raise ValueError('Only univariate analysis are supported right now.')

        stats = self.process_gemma_out(gemma_output)

        return gemma, stats
