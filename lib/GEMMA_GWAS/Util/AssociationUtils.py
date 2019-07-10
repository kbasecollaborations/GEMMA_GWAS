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
        #self.wsc = Workspace(config['SDK_CALLBACK_URL'])
        self.wsc = Workspace("https://appdev.kbase.us/services/ws")
        self.scratch = config["scratch"]
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

    def _get_iids(self, col_attr_map_ref):
        iids = {}

        col_attr_map_ids = self.dfu.get_objects({'object_refs': [col_attr_map_ref]})['data'][0]['data']['instances']

        for id_dict in col_attr_map_ids:
            iids[col_attr_map_ids[id_dict][0]] = col_attr_map_ids[id_dict][1]

        return iids

    def _get_fid_master_list(self, row_attribute_ref):
        row_attr_ids = self.dfu.get_objects({'object_refs': [row_attribute_ref]})['data'][0]['data']['instances']
        fids = list(row_attr_ids.keys())

        if not fids:
            raise ValueError("Error retrieving fids from sample attribute object")

        return fids

    def _mk_phenos_from_trait_matrix_uni(self, trait_matrix_ref):
        """
        This function makes 3 column space delimited phenotype files from the Kbase trait matrix object.

        This function is for univariate GWAS only, a single phenotype file is created for each phenotype.

        The columns in the file are: Family ID, Within-family ID, Phenotypic value.

        :param trait_matrix_ref:
            Kbase trait matrix reference. Data referenced must be of type KBaseMatrices.TraitMatrix
        :return: phenotype_files
            Returns a dictionary of phenotype type file locations, and phenotype labels

        TODO: fix iids, if they're null use fid instead DO NOT LEAVE BLANK will disrupt plink operation

        """

        trait_matrix_obj = self.dfu.get_objects({'object_refs': [trait_matrix_ref]})['data'][0]

        logging.info("Parsing trait matrix.")

        ws_type = trait_matrix_obj['info'][2]
        if 'KBaseMatrices.TraitMatrix' in ws_type:
            fids = trait_matrix_obj['data']['data']['col_ids']
            iids = self._get_iids(trait_matrix_obj['data']['col_attributemapping_ref'])
            fid_master_list = self._get_fid_master_list(trait_matrix_obj['data']['col_attributemapping_ref'])

            phenotypeids = trait_matrix_obj['data']['data']['row_ids']

            phenotype_files = {}

            for x in range(0, len(phenotypeids)):
                phenotype_files[x] = {}

                phenotype_files[x]['id'] = phenotypeids[x]
                phenotypevals = trait_matrix_obj['data']['data']['values'][x]

                for fid in fids:
                    fids_not_documented = []
                    if fid not in fid_master_list:
                        fids_not_documented.append(fids_not_documented)

                    if fids_not_documented:
                        fids_not_listed = ', '.join(fids_not_documented)
                        raise ValueError('Fids: ' + fids_not_listed +'. Are not listed in the' \
                                         ' sample attribute mapping meta information object.')

                single_pheno_file_path = os.path.join(self.scratch, 'pheno'+str(x)+'.txt')

                pheno_tsv_headers = "FID IID " + phenotypeids[x] + "\n"

                with open(single_pheno_file_path, 'w') as f:
                    f.write(pheno_tsv_headers)
                    for k in range(0, len(fids)):
                        if phenotypevals[k] == 1:
                            phenoval = '1'
                        elif phenotypevals[k] == 0:
                            phenoval = '0'
                        elif str(phenotypevals[k]).upper() == 'NONE':
                            phenoval = 'NA'
                        else:
                            phenoval = str(phenotypevals[k])

                        f.write(fids[k] + " " + fids[k] + " " + phenoval + "\n")
                        """
                        if iids[fids[k]]:
                            f.write(fids[k] + " " + iids[fids[k]] + " " + str(phenotypevals[k]) + "\n")
                        else:
                            f.write(fids[k] + " NA " + str(phenotypevals[k]) + "\n")
                        """
                    f.close()
                    phenotype_files[x]['file'] = single_pheno_file_path

            logging.info("Finished trait matrix parsing.")

            return phenotype_files
        else:
            raise ValueError('Cannot write data to VCF; invalid WS type (' + ws_type +
                             ').  Supported types is KBaseMatrices.TraitMatrix')

    def _mk_phenos_from_trait_matrix_multi(self, trait_matrix_ref):
        """
        This function makes 3 column space delimited phenotype files from the Kbase trait matrix object.

        This function is for multivariate GWAS only, a single phenotype file is created for each phenotype.

        The columns in the file are: Family ID, Within-family ID, phenotype-n1, phenotype-n2...

        :param trait_matrix_ref:
            Kbase trait matrix reference. Data referenced must be of type KBaseMatrices.TraitMatrix
        :return: phenotype_files
            Returns a dictionary of phenotype type file locations, and phenotype labels

        TODO: fix iids, if they're null use fid instead DO NOT LEAVE BLANK will disrupt plink operation

        """

        trait_matrix_obj = self.dfu.get_objects({'object_refs': [trait_matrix_ref]})['data'][0]

        logging.info("Parsing trait matrix.")

        ws_type = trait_matrix_obj['info'][2]
        if 'KBaseMatrices.TraitMatrix' in ws_type:
            fids = trait_matrix_obj['data']['data']['col_ids']
            iids = self._get_iids(trait_matrix_obj['data']['col_attributemapping_ref'])
            fid_master_list = self._get_fid_master_list(trait_matrix_obj['data']['col_attributemapping_ref'])

            phenotypeids = trait_matrix_obj['data']['data']['row_ids']

            phenotype_files = {'multi': {}}
            phenotype_files['multi']['id'] = ' '.join(phenotypeids)

            for fid in fids:
                fids_not_documented = []
                if fid not in fid_master_list:
                    fids_not_documented.append(fids_not_documented)

                if fids_not_documented:
                    fids_not_listed = ', '.join(fids_not_documented)
                    raise ValueError('Fids: ' + fids_not_listed + '. Are not listed in the' \
                                                                  ' sample attribute mapping meta information object.')

            phenotype_files['multi']['pheno'] = os.path.join(self.scratch, 'multi.pheno')
            phenotypevals = trait_matrix_obj['data']['data']['values']

            with open(phenotype_files['multi']['pheno'], 'w') as f:
                headerline = "FID IID "+phenotype_files['multi']['id']
                f.write(headerline+"\n")

                for k in range(0, len(fids)):
                    fline = ""
                    # This line is the beginning of the .fam file for plink
                    # phenotype files, instead we are using bimbam phenotype
                    # file format for gemma multivariate analysis
                    # fline = fids[k] + " " + fids[k] + " 0 0 0 "
                    fline = fids[k] + " " + fids[k]

                    for j in range(0, len(phenotypevals)):
                        if str(phenotypevals[j][k]).rstrip().lstrip().upper() is 'NONE':
                            fline += " NA"
                        elif phenotypevals[j][k] == None:
                            fline += " NA"
                        elif phenotypevals[j][k] == 1:
                            fline += " 1"
                        elif phenotypevals[j][k] == 0:
                            fline += " 0"
                        else:
                            fline += " "+str(phenotypevals[j][k])
                    f.write(fline+"\n")
                f.close()

            logging.info("Finished parsing trait matrix")

            return phenotype_files
        else:
            raise ValueError('Cannot write data to VCF; invalid WS type (' + ws_type +
                             ').  Supported types is KBaseMatrices.TraitMatrix')

    def _mk_centered_kinship(self, plink_prefixes):
        if 'multi' in plink_prefixes:
            logging.info("Generating kinship matrix")

            # multivariate analysis
            kin_cmd = ['/usr/bin/time', '-v', 'gemma', '-bfile', plink_prefixes['multi']['plink'], '-gk', '1', '-o',
                       'kinship-multi']
            try:
                proc = subprocess.Popen(kin_cmd, cwd=self.scratch)
                proc.wait()
            except Exception as e:
                exit(e)

            logging.info("Kinship matrix done")

            plink_prefixes['multi']['kinship'] = os.path.join(self.scratch, 'output', 'kinship-multi.cXX.txt')

            if not os.path.exists(plink_prefixes['multi']['kinship']):
                exit("Kinship file does not exist: " + plink_prefixes['multi']['kinship'])
        else:
            # univariate analysis
            kinship_base_prefix = 'kinship'
            logging.info("Generating kinship matrices")

            for x in range(0, len(plink_prefixes)):
                kin_cmd = ['/usr/bin/time', '-v', 'gemma', '-bfile', plink_prefixes[x]['plink'], '-gk', '1', '-o', kinship_base_prefix+str(x)]

                try:
                    proc = subprocess.Popen(kin_cmd, cwd=self.scratch)
                    proc.wait()
                except Exception as e:
                    exit(e)

                plink_prefixes[x]['kinship'] = os.path.join(self.scratch, 'output', kinship_base_prefix+str(x)+'.cXX.txt')

                if not os.path.exists(plink_prefixes[x]['kinship']):
                    exit("Kinship file does not exist: "+plink_prefixes[x]['kinship'])

            logging.info("Kinship matrices generation done")

        return plink_prefixes

    def _mk_standardized_kinship(self, assoc_info):
        logging.info("Generating kinship matrix.")
        kinship_cmd = ['/usr/bin/time', '-v', 'gemma', '-bfile', os.path.join(self.scratch, assoc_info['plink']), '-gk', '2', '-o', str(assoc_info['id'])]
        kinship_file = os.path.join(self.scratch, 'output', str(assoc_info['id'])+'.sXX.txt')

        try:
            proc = subprocess.Popen(kinship_cmd, cwd=self.scratch)
            proc.wait()
        except Exception as e:
            exit(e)

        logging.info("Kinship matrix generation done")

        if not os.path.exists(kinship_file):
            raise ValueError('Unable to create Standardized Kinship Matrix. Contact Administrator')

        return kinship_file

    def _check_pheno_case_control(self, phenotypefile):
        # phenofile can be assumed to be a 3 column space-delimted file since,
        # this method is only used in the univariate analysis

        phenovals = []
        with open(phenotypefile, 'r') as f:
            pheno = csv.reader(f, delimiter=' ')
            # skip headers
            next(pheno)
            for row in pheno:
                phenovals.append(row[2])

        phenoset = set(phenovals)

        if phenoset == {'0', '1', 'NA'} \
                or phenoset == {'0', '1'}:
                # TODO: make decision on case/control detection
                #or phenoset == {'0', '1',  '2', 'NA'} \
                #or phenoset == {'0', '1',  '2'}:
            return '--1'
        else:
            return False

    def _mk_plink_bin_uni(self, phenotypes):
        plink_base_prefix = 'plink_variation'
        plink_prefixes = []

        logging.info('Encoding with plink.')

        for x in range(0, len(phenotypes)):
            plinkvars = ['--make-bed', '--vcf', self.varfile, '--pheno', phenotypes[x]['file'],
                         '--allow-no-sex', '--allow-extra-chr', '--output-chr', 'chr26']

            cc_flag = self._check_pheno_case_control(phenotypes[x]['file'])
            if cc_flag:
                plinkvars.append(cc_flag)

            plinkvars.append('--out')
            plinkvars.append(plink_base_prefix+str(x))

            plinkcmd = ['/usr/bin/time', '-v', 'plink']

            for arg in plinkvars:
                plinkcmd.append(arg)

            try:
                proc = subprocess.Popen(plinkcmd, cwd=self.scratch)
                proc.wait()
            except Exception as e:
                exit(e)

            phenotypes[x]['plink'] = plink_base_prefix+str(x)
            plink_bed = os.path.join(self.scratch, plink_base_prefix+str(x)+'.bed')
            plink_bim = os.path.join(self.scratch, plink_base_prefix+str(x)+'.bim')
            plink_fam = os.path.join(self.scratch, plink_base_prefix+str(x)+'.fam')

            if not os.path.exists(plink_bed):
                print(f'\n{plink_prefixes}')
                raise IOError(f"Plink bed doesn't exist: {plink_bed}")

            if not os.path.exists(plink_bim):
                raise IOError(f"Plink bim doesn't exist: {plink_bim}")

            if not os.path.exists(plink_fam):
                raise IOError(f"Plink fam doesn't exist: {plink_fam}")

        logging.info('Plink encoding finished')

        return phenotypes

    def _mk_plink_bin_multi(self, phenotypes):
        """
        :param phenotypes:
        :return:
        """

        """
        For multivariate analysis need to use plink2 and plink2 flags, this plink1.9 call will no work
        and causes kinship matrix issues
        """

        plinkvars = ['--make-bed', '--vcf', self.varfile, '--pheno', phenotypes['multi']['pheno'], '--double-id',
                     '--vcf-require-gt', '--max-alleles', '2', '--allow-extra-chr', '--output-chr', 'chr26',
                     '--out', 'plink_multi']
        plinkcmd = ['/usr/bin/time', '-v', 'plink2']

        logging.info('Encoding with plink2.')

        for arg in plinkvars:
            plinkcmd.append(arg)

        try:
            proc = subprocess.Popen(plinkcmd, cwd=self.scratch)
            proc.wait()
        except Exception as e:
            exit(e)

        phenotypes['multi']['plink'] = 'plink_multi'
        plink_bed = os.path.join(self.scratch, 'plink_multi.bed')
        plink_bim = os.path.join(self.scratch, 'plink_multi.bim')

        if not os.path.exists(plink_bed):
            raise IOError('Plink bed doesn\'t exist')

        if not os.path.exists(plink_bim):
            raise IOError('Plink bim doesn\'t exist')

        logging.info('Plink2 encoding finished')

        return phenotypes

    def run_gemma_assoc_uni(self, kinmatrix):
        if 'multi' in kinmatrix:
            raise ValueError('Attempted to run a univariate gemma analysis on a multivariate dataset')

        logging.info('Running univariate association analysis with GEMMA')

        for x in range(0, len(kinmatrix)):
            assoc_base_file_prefix = 'gemma_assoc'
            assoc_args = ['-bfile', kinmatrix[x]['plink'], '-k', kinmatrix[x]['kinship'], '-lmm', '4', '-debug', '-o',
                          assoc_base_file_prefix + kinmatrix[x]['id']]
            assoc_cmd = ['/usr/bin/time', '-v', 'gemma']

            for arg in assoc_args:
                assoc_cmd.append(arg)

            assoc_results = kinmatrix

            try:
                #proc = subprocess.Popen(assoc_cmd, cwd=self.scratch, stdout=subprocess.PIPE)
                proc = subprocess.Popen(assoc_cmd, cwd=self.scratch)
                proc.wait()
                #out, err = proc.communicate()

                if proc.returncode is -2:
                    # brent error
                    newkinship = self._mk_standardized_kinship(kinmatrix[x])
                    new_assoc_cmd = ['/usr/bin/time', '-v', 'gemma','-bfile', kinmatrix[x]['plink'], '-k', newkinship, '-lmm', '4','-debug', '-o',
                                     assoc_base_file_prefix + kinmatrix[x]['id']]

                    try:
                        newproc = subprocess.Popen(new_assoc_cmd, cwd=self.scratch, stdout=subprocess.PIPE)
                        newproc.wait()
                        out, err = newproc.communicate()
                    except Exception as e:
                        exit(e)

                    if not newproc.returncode is -2:
                        assoc_results[x]['gemma'] = os.path.join(self.scratch, 'output', assoc_base_file_prefix +
                                                                 kinmatrix[x]['id']+'.assoc.txt')
                    else:
                        logging.error('Unspecified subprocess execution error' + str(newproc))
                        raise IOError('GEMMA Association failed.')
                else:
                    assoc_results[x]['gemma'] = os.path.join(self.scratch, 'output',assoc_base_file_prefix +
                                                             kinmatrix[x]['id']+'.assoc.txt')
            except Exception as e:
                logging.error('Unspecified subprocess execution error' + str(proc))
                exit(e)

        #decode_out = out.decode('utf-8')
        #decode_out = decode_out.split('\n')
        decode_out = 'abcdef'

        if len(kinmatrix) is 1:
            logging.info('GEMMA univariate association analysis complete')
        else:
            logging.info('GEMMA univariate association analyses complete')

        return assoc_results, decode_out

    def run_gemma_assoc_multi(self, kinmatrix):
        if 'multi' not in kinmatrix:
            raise ValueError('Attempted to run a multivariate gemma analysis on a univariate dataset')

        assoc_args = ['-bfile', kinmatrix['multi']['plink'], '-k', kinmatrix['multi']['kinship'],
                      '-lmm', '4', '-debug', '-o', 'gemma_multi_assoc']
        assoc_cmd = ['/usr/bin/time', '-v', 'gemma']

        logging.info('Running multivariate association analysis with GEMMA')

        for arg in assoc_args:
            assoc_cmd.append(arg)

        assoc_results = kinmatrix

        try:
            proc = subprocess.Popen(assoc_cmd, cwd=self.scratch, stdout=subprocess.PIPE)
            proc.wait()
            out, err = proc.communicate()
            if proc.returncode is -2:
                # brent error
                newkinship = self._mk_standardized_kinship(kinmatrix['multi'])
                new_assoc_cmd = ['gemma', '-bfile', kinmatrix['multi']['plink'], '-k', newkinship, '-lmm', '4', '-debug',
                                 '-o', 'gemma_multi_assoc']
                try:
                    newproc = subprocess.Popen(new_assoc_cmd, cwd=self.scratch, stdout=subprocess.PIPE)
                    newproc.wait()
                    out, err = newproc.communicate()
                except Exception as e:
                    exit(e)
                if not newproc.returncode is -2:
                    assoc_results[x]['gemma'] = os.path.join(self.scratch, 'output', 'gemma_multi_assoc.assoc.txt')
                else:
                    logging.error('Failed to run gemma association:\n' + str(newproc))
                    raise RuntimeError('GEMMA Association failed.')
            else:
                assoc_results['multi']['gemma'] = os.path.join(self.scratch, 'output', 'gemma_multi_assoc.assoc.txt')
        except Exception as e:
            exit(e)

        logging.info('GEMMA multivariate association analyses complete')

        decode_out = out.decode('utf-8')
        decode_out = decode_out.split('\n')

        return assoc_results, decode_out

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
            phenotype = self._mk_phenos_from_trait_matrix_uni(params['trait_matrix'])
            plink = self._mk_plink_bin_uni(phenotype)
            kinmatrix = self._mk_centered_kinship(plink)
            gemma, gemma_output = self.run_gemma_assoc_uni(kinmatrix)
        else:
            # mutlivariate analysis
            phenotype = self._mk_phenos_from_trait_matrix_multi(params['trait_matrix'])
            plink = self._mk_plink_bin_multi(phenotype)
            kinmatrix = self._mk_centered_kinship(plink)
            gemma, gemma_output = self.run_gemma_assoc_multi(kinmatrix)

        stats = self.process_gemma_out(gemma_output)

        return gemma, stats
