import subprocess
import os
import operator
import shutil
import csv
import logging
from pprint import pprint as pp
import pandas as pd
import numpy as np
import hashlib

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class AssociationUtils:
    def __init__(self, config, varfile, plink_prefix = 'plink_vcf'):
        self.dfu = DataFileUtil(config['SDK_CALLBACK_URL'])
        # TODO: input variable for workspace url
        self.wsc = Workspace("https://appdev.kbase.us/services/ws")
        self.scratch = config["scratch"]
        self.plink_pref = plink_prefix
        self.state = {}
        self._process_varfiles(varfile)

    def _process_varfiles(self, file):
        file_name, file_ext = os.path.splitext(file)
        if file_ext == '.vcf' or file_ext == '.gz':
            if os.path.exists(file):
                self.varfile = file
                self.state['vcf'] = {}
                self.state['vcf']['file'] = self.varfile
                self.state['vcf']['md5'] = hashlib.md5(open(self.varfile, 'rb').read()).hexdigest()
            else:
                raise IOError('Variation file provided does not exist or is not readable.')
        else:
            raise ValueError('Only supporting VCF variation as input right now.')

    def _mk_plink_bin_uni(self):
        logging.info("Generating PLINK binary files.")
        plinkvars = ['--make-bed', '--vcf', self.varfile, '--allow-no-sex','--allow-extra-chr', '--out', self.plink_pref]
        plinkcmd = ['plink']

        for arg in plinkvars:
            plinkcmd.append(arg)

        try:
            proc = subprocess.Popen(plinkcmd, cwd=self.scratch, stdout=subprocess.PIPE, close_fds=True)
            # do not use proc.wait(): https://docs.python.org/3/library/subprocess.html#subprocess.Popen.wait
            # use proc.communicate() instead:
            # and close_fds=True
            out, err = proc.communicate()
        except Exception as e:
            exit(e)

        plink_bed = os.path.join(self.scratch, self.plink_pref + '.bed')
        plink_bim = os.path.join(self.scratch, self.plink_pref + '.bim')
        plink_fam = os.path.join(self.scratch, self.plink_pref + '.fam')

        self.state['plink'], self.state['plink']['plink_bed'], self.state['plink']['plink_bim'], \
            self.state['plink']['plink_fam'] = {}, {}, {}, {}

        if not os.path.exists(plink_bed):
            raise FileNotFoundError('Plink bed doesn\'t exist')
        else:
            self.state['plink']['plink_bed']['file'] = plink_bed
            self.state['plink']['plink_bed']['md5'] = hashlib.md5(open(plink_bed, 'rb').read()).hexdigest()

        if not os.path.exists(plink_bim):
            raise FileNotFoundError('Plink bim doesn\'t exist')
        else:
            self.state['plink']['plink_bim']['file'] = plink_bim
            self.state['plink']['plink_bim']['md5'] = hashlib.md5(open(plink_bim, 'rb').read()).hexdigest()


        if os.path.exists(plink_fam):
            self.plink_fam_template = shutil.move(plink_fam, os.path.join(self.scratch, 'plink_fam_template.fam'))
            self.state['plink']['plink_fam']['file'] = self.plink_fam_template
            self.state['plink']['plink_fam']['md5'] = hashlib.md5(open(self.plink_fam_template, 'rb').read()).hexdigest()
        else:
            raise FileNotFoundError('Plink fam doesn\'t exist')

        logging.info("Done making PLINK binaries")

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

        logging.info('Done parsing trait matrix')

        print(phenosdict)
        return phenosdict

    def _mk_phenos_from_am(self, trait_matrix_ref, selected_traits):
        #TODO: Add object type check and few other checks for the structure of the object
        #config = self.cfg
        #self.scratch = config["scratch"]
        #self.callback_url = config["SDK_CALLBACK_URL"]
        #self.dfu = DataFileUtil(self.callback_url)
        trait_matrix_obj = self.dfu.get_objects({'object_refs': [trait_matrix_ref]})['data'][0]['data']
        attributes = trait_matrix_obj['attributes']
        instances = trait_matrix_obj['instances']



        phenosdict = {}
        for count, trait in enumerate(attributes):
            if (trait['attribute'] not in selected_traits):
                continue
                
            phenodict = {}
            trait_label = trait['attribute']
            for key in instances:
                phenodict[key] = instances[key][count]
            phenosdict[trait_label] = phenodict
        print(phenosdict)
        return (phenosdict)

    def mk_fam_files_from_phenos(self, phenovals):
        logging.info('Making .fam files from phenotype values')
        self.fam_directory = os.path.join(self.scratch, 'fams')

        if not os.path.isdir(self.fam_directory):
            os.mkdir(os.path.join(self.scratch, 'fams'))

        print ("Osdddddddddddddddddddddddddddddddddddddddddddddddddddddd")

        fam_template = pd.read_csv(self.plink_fam_template, sep=' ',
                                   names=['fid', 'iid', 'iidf', 'iidm', 'sex', 'value'], index_col=False,
                                   dtype={'fid': object, 'iid': object, 'iidf': np.int64,
                                          'iidm': np.int64, 'sex': np.int64, 'value': np.int64})
        fam_template = fam_template.replace(-9, 'NA')
        fam_files = []

        for pheno, values in phenovals.items():
            new_fam = fam_template  # copy template data frame
            # pheno is the phenotype name
            # values is a dictionary with keys of sampleid and values of phenotypic value
            for sample, value in values.items():
                if not value:
                    value = "NA"
                new_fam.loc[new_fam['fid'] == sample] = sample, sample, 0, 0, 0, value

            new_fam_path = os.path.join(self.fam_directory, pheno + '.fam')
            new_fam.fillna(value='NA', inplace=True)
            new_fam.to_csv(new_fam_path, sep=' ', header=None, index=False)
            self.state[pheno], self.state[pheno]['fam'] = {}, {}
            self.state[pheno]['fam']['file'] = new_fam_path
            self.state[pheno]['fam']['md5'] = hashlib.md5(open(new_fam_path, 'rb').read()).hexdigest()
            fam_files.append(new_fam_path)

        logging.info('Done making fam files from phenotype values')
        return fam_files

    def mk_centered_kinship_uni(self, phenovalues, famfiles):
        # univariate analysis
        self.kinship_base_prefix = 'kinship'
        logging.info("Generating kinship matrices")
        kinshipmatricies = []

        for pheno in phenovalues:
            kin_cmd = ['gemma', '-bfile', self.plink_pref, '-gk', '1', '-o', self.kinship_base_prefix + '_' + pheno]

            print (kin_cmd)
            famfile = os.path.join(self.fam_directory, pheno + '.fam')

            if famfile not in famfiles:
                raise FileNotFoundError(f'{famfile} does not exist in famfile list: {str(famfiles)}')

            self._stage_fam_file(famfile)

            try:
                proc = subprocess.Popen(kin_cmd, cwd=self.scratch, stdout=subprocess.PIPE, close_fds=True)
                # do not use proc.wait() : https://docs.python.org/3/library/subprocess.html#subprocess.Popen.wait
                # use proc.communicate() instead:
                # and close_fds=True
                out, err = proc.communicate()
            except Exception as e:
                logging.error('Centered kinship generation failed')
                raise ChildProcessError(e)

            kinship_file = os.path.join(self.scratch, 'output', self.kinship_base_prefix + '_' + pheno + '.cXX.txt')

            self.state[pheno]['kinship'] = {}
            self.state[pheno]['kinship']['file'] = kinship_file
            self.state[pheno]['kinship']['md5'] = hashlib.md5(open(kinship_file, 'rb').read()).hexdigest()

            if not os.path.exists(kinship_file):
                raise FileNotFoundError("Kinship file does not exist: " +
                    os.path.join(self.scratch, 'output', self.kinship_base_prefix + '_' + pheno + '.cXX.txt'))
            else:
                kinshipmatricies.append(os.path.join(self.scratch, 'output',
                    self.kinship_base_prefix + '_' + pheno + '.cXX.txt'))

        logging.info("Kinship matrices generation done")

        return kinshipmatricies

    def _stage_fam_file(self, famfile):
        if os.path.isfile(os.path.join(self.scratch, self.plink_pref + '.fam')):
            os.remove(os.path.join(self.scratch, self.plink_pref + '.fam'))

        shutil.copyfile(famfile, os.path.join(self.scratch, self.plink_pref + '.fam'))

        return True

    def run_gemma_assoc_uni(self, kinship_files, famfiles, phenotypes, plink_prefix):
        logging.info("Running GEMMA Univariate lmm association tests")
        self.assoc_base_file_prefix = 'gemma_assoc'
        gemma_stats = {}
        gemma_files = {}
        for pheno in phenotypes:
            famfile = os.path.join(self.fam_directory, pheno + '.fam')
            kinship = os.path.join(self.scratch, 'output', self.kinship_base_prefix + '_' + pheno + '.cXX.txt')

            if famfile not in famfiles:
                raise FileNotFoundError(f'{famfile} does not exist in famfile list: {str(famfiles)}')

            if not self._stage_fam_file(famfile):
                raise FileNotFoundError(f'Fam file: {famfile} does not exist.')

            if kinship not in kinship_files:
                raise FileNotFoundError(f'Kinship file: {kinship} does not exist.')
                
            assoc_args = ['-bfile', plink_prefix, '-k', kinship, '-lmm', '1', '-o',
                          self.assoc_base_file_prefix + '_' + pheno]
            assoc_cmd = ['gemma']

            for arg in assoc_args:
                assoc_cmd.append(arg)
          
            try:
                proc = subprocess.Popen(assoc_cmd, cwd=self.scratch, stdout=subprocess.PIPE, close_fds=True)
                # do not use proc.wait() : https://docs.python.org/3/library/subprocess.html#subprocess.Popen.wait
                # use proc.communicate() instead:
                # and close_fds=True
                out, err = proc.communicate()
            except Exception as e:
                logging.error('Unspecified subprocess execution error' + str(err.decode('UTF-8')))
                raise ChildProcessError(e)

            decode_out = out.decode('utf-8')
            decode_out = decode_out.split('\n')
            decode_out = self.process_gemma_out(decode_out)

            gemma_stats[pheno] = decode_out
            gemma_out_file = os.path.join(self.scratch, 'output',
                                          self.assoc_base_file_prefix + '_' + pheno + '.assoc.txt')
            gemma_files[pheno] = gemma_out_file

            self.state[pheno]['gemma'] = {}
            self.state[pheno]['gemma']['file'] = gemma_out_file
            self.state[pheno]['gemma']['md5'] = hashlib.md5(open(gemma_out_file, 'rb').read()).hexdigest()
            self.state[pheno]['gemma']['stats'] = decode_out

            if not os.path.exists(gemma_out_file):
                raise FileNotFoundError(f'GEMMA association output not found: {gemma_out_file}')

        logging.info("GEMMA association testing done")
        return gemma_files, gemma_stats

    def process_gemma_out(self, output):
        stats = {}
        for line in output:
            # number of analyzed individuals = XX
            if 'number of analyzed individuals' in line:
                spline = line.split(' ')
                stats['individuals'] = spline[6]

            # pve estimate =
            if 'pve estimate =' in line:
                spline = line.split('=')
                stats['pve'] = spline[1]

            # se(pve) =
            if 'se(pve)' in line:
                spline = line.split('=')
                stats['se_pve'] = spline[1]

        return stats

    def run_assoc_exp(self, params):
        if params['model'] is 0:
            # univariate analysis
            plink = self._mk_plink_bin_uni()
            #phenovals = self._mk_phenos_from_trait_matrix_uni(params['trait_matrix'])
            phenovals = self._mk_phenos_from_am(params['trait_matrix'], params['selected_traits'])
            famfiles = self.mk_fam_files_from_phenos(phenovals)
            kinmatricies = self.mk_centered_kinship_uni(phenovals, famfiles)
            gemma, gemma_output = self.run_gemma_assoc_uni(kinmatricies, famfiles, phenovals, plink)
        else:
            raise NotImplementedError('Only univariate analysis are supported right now.')

        return self.state
