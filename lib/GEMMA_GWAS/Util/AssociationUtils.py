import subprocess
import os
import operator
import shutil
import csv
import logging
from pprint import pprint as pp
import pandas as pd

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class AssociationUtils:
    def __init__(self, config, varfiles, plink_prefix = 'plink_vcf'):
        self.dfu = DataFileUtil(config['SDK_CALLBACK_URL'])
        # TODO: input variable for workspace url
        self.wsc = Workspace("https://appdev.kbase.us/services/ws")
        self.scratch = config["scratch"]
        self.plink_pref = plink_prefix
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

        for arg in plinkvars:
            plinkcmd.append(arg)

        try:
            proc = subprocess.Popen(plinkcmd, cwd=self.scratch)
            proc.wait()
        except Exception as e:
            exit(e)

        plink_bed = os.path.join(self.scratch, self.plink_pref + '.bed')
        plink_bim = os.path.join(self.scratch, self.plink_pref + '.bim')
        plink_fam = os.path.join(self.scratch, self.plink_pref + '.fam')

        if not os.path.exists(plink_bed):
            raise FileNotFoundError('Plink bed doesn\'t exist')

        if not os.path.exists(plink_bim):
            raise FileNotFoundError('Plink bim doesn\'t exist')

        if os.path.exists(plink_fam):
            self.plink_fam_template = shutil.move(plink_fam, os.path.join(self.scratch, 'plink_fam_template.fam'))
        else:
            raise FileNotFoundError('Plink fam doesn\'t exist')

        return self.plink_pref

    def _mk_phenos_from_trait_matrix_uni(self, trait_matrix_ref):
        logging.info("Parsing trait matrix.")

        trait_matrix_obj = self.dfu.get_objects({'object_refs': [trait_matrix_ref]})['data'][0]
        ws_type = trait_matrix_obj['info'][2]

        if 'KBaseMatrices.TraitMatrix' in ws_type:
            phenotypes = trait_matrix_obj['data']['data']['row_ids']
            fids = trait_matrix_obj['data']['data']['col_ids']
            phenosdict = {}

            # create dictionary of phenotype values
            # first level key is phenotype name, with dictionary value
            # second level key is sample id, with values of phenotype value
            for pheno in phenotypes:
                phenodict = {}
                values = trait_matrix_obj['data']['data']['values'][phenotypes.index(pheno)]
                for x in range(0, len(fids)):
                    phenodict[fids[x]] = values[x]
                phenosdict[pheno] = phenodict
        else:
            raise ValueError('Cannot write data to VCF; invalid WS type (' + ws_type +
                             ').  Supported types is KBaseMatrices.TraitMatrix')

        return phenosdict

    def mk_fam_files_from_phenos(self, phenovals):
        self.fam_directory = os.path.join(self.scratch, 'fams')

        if not os.path.isdir(self.fam_directory):
            os.mkdir(os.path.join(self.scratch, 'fams'))

        fam_template = pd.DataFrame.from_csv(self.plink_fam_template, sep=' ', header=None, index_col=False)
        fam_template = fam_template.replace(-9, 'NA')
        fids = list(fam_template[0])
        fam_files = []

        for pheno, values in phenovals.items():
            new_fam = fam_template

            for sample, val in values.items():
                if sample in fids:
                    new_fam.loc[new_fam[0] == sample] = [sample, sample, 0, 0, 0, val]

            new_fam_path = os.path.join(self.fam_directory, pheno+'.fam')
            new_fam.fillna(value='NA', inplace=True)
            new_fam.to_csv(new_fam_path, sep=' ', header=None, index=False)
            fam_files.append(new_fam_path)

        return fam_files

    def mk_centered_kinship(self, plink_prefix):

    def run_gemma_assoc_uni(self, kinship, phenovals):
        return 0

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
            plink = self._mk_plink_bin_uni()
            phenovals = self._mk_phenos_from_trait_matrix_uni(params['trait_matrix'])
            famfiles = self.mk_fam_files_from_phenos(phenovals)
            kinmatrix = self.mk_centered_kinship(plink)
            gemma, gemma_output = self.run_gemma_assoc_uni(kinmatrix, famfiles)
        else:
            raise ValueError('Only univariate analysis are supported right now.')

        stats = self.process_gemma_out(gemma_output)

        return gemma, stats
