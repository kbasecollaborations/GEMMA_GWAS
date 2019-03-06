import subprocess
import os
import operator
import shutil
from pprint import pprint as pp

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
                        f.write(fids[k] + " " + fids[k] + " " + str(phenotypevals[k]) + "\n")
                        """
                        if iids[fids[k]]:
                            f.write(fids[k] + " " + iids[fids[k]] + " " + str(phenotypevals[k]) + "\n")
                        else:
                            f.write(fids[k] + " NA " + str(phenotypevals[k]) + "\n")
                        """
                    f.close()
                    phenotype_files[x]['file'] = single_pheno_file_path

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

            phenotype_files['multi']['file'] = os.path.join(self.scratch, 'multi.fam')
            phenotypevals = trait_matrix_obj['data']['data']['values']

            with open(phenotype_files['multi']['file'], 'w') as f:
                for k in range(0, len(fids)):
                    famline = fids[k] + " " + fids[k] + " 0 0 0 "
                    for j in range(0, len(phenotypevals)):
                        if str(phenotypevals[j][k]).upper() is 'NONE':
                            famline += " -9"
                        else:
                            famline += " " + str(phenotypevals[j][k])
                    f.write(famline+"\n")
                f.close()

            return phenotype_files
        else:
            raise ValueError('Cannot write data to VCF; invalid WS type (' + ws_type +
                             ').  Supported types is KBaseMatrices.TraitMatrix')

    def _mk_centered_kinship(self, plink_prefixes):
        if 'multi' in plink_prefixes:
            # multivariate analysis
            kin_cmd = ['gemma', '-bfile', os.path.join(self.scratch, plink_prefixes['multi']['plink']), '-gk', '1', '-o',
                       'kinship-multi']
            try:
                proc = subprocess.Popen(kin_cmd, cwd=self.scratch)
                proc.wait()
            except Exception as e:
                exit(e)

            plink_prefixes['multi']['kinship'] = os.path.join(self.scratch, 'output', 'kinship-multi.cXX.txt')

            if not os.path.exists(plink_prefixes['multi']['kinship']):
                exit("Kinship file does not exist: " + plink_prefixes[x]['kinship'])
        else:
            # univariate analysis
            kinship_base_prefix = 'kinship'

            for x in range(0, len(plink_prefixes)):
                kin_cmd = ['gemma', '-bfile', os.path.join(self.scratch, plink_prefixes[x]['plink']), '-gk', '1', '-o', kinship_base_prefix+str(x)]

                try:
                    proc = subprocess.Popen(kin_cmd, cwd=self.scratch)
                    proc.wait()
                except Exception as e:
                    exit(e)

                plink_prefixes[x]['kinship'] = os.path.join(self.scratch, 'output', kinship_base_prefix+str(x)+'.cXX.txt')

                if not os.path.exists(plink_prefixes[x]['kinship']):
                    exit("Kinship file does not exist: "+plink_prefixes[x]['kinship'])

        return plink_prefixes

    def _mk_standardized_kinship(self, assoc_info):
        kinship_cmd = ['gemma', '-bfile', os.path.join(self.scratch, assoc_info['plink']), '-gk', '2', '-o', str(assoc_info['id'])]
        kinship_file = os.path.join(self.scratch, 'output', str(assoc_info['id'])+'.sXX.txt')

        try:
            proc = subprocess.Popen(kinship_cmd, cwd=self.scratch)
            proc.wait()
        except Exception as e:
            exit(e)

        if not os.path.exists(kinship_file):
            raise ValueError('Unable to create Standardized Kinship Matrix. Contact Administrator')

        return kinship_file

    def _mk_plink_bin_uni(self, phenotypes):
        plink_base_prefix = 'plink_variation'
        """
            plink flags for .ped .map and phenotype tsv into plink
            plinkvars = ['--make-bed','--ped',self.local_ped_file,'--map',self.local_map_file,'--pheno',
                    self.local_pheno_file,'--allow-no-sex','--chr','1','--out',self.local_plink_prefix]
        """
        plink_prefixes = []
        for x in range(0, len(phenotypes)):
            plinkvars = ['--make-bed', '--vcf', self.varfile, '--pheno',
                         phenotypes[x]['file'], '--allow-no-sex','--allow-extra-chr', '--out', plink_base_prefix+str(x)]
            plinkcmd = ['plink']

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
                raise IOError('Plink bed doesn\'t exist')
            else:
                plink_prefixes.append(plink_base_prefix+str(x))
                print("----pheno type: "+os.path.basename(phenotypes[x]['file'])[:-4]+"----")
                print("----plink binary generated----")
                print("----" + plink_bed + "----")

            if not os.path.exists(plink_bim):
                raise IOError('Plink bim doesn\'t exist')
            else:
                print("----" + plink_bim + "----")

            if not os.path.exists(plink_fam):
                raise IOError('Plink fam doesn\'t exist')
            else:
                print("----" + plink_fam + "----\n")

        return phenotypes

    def _mk_plink_bin_multi(self, phenotypes):
        plinkvars = ['--make-bed', '--vcf', self.varfile, '--pheno',
                     phenotypes['multi']['file'], '--all-pheno', '--allow-no-sex', '--allow-extra-chr', '--out',
                     'plink_multi']
        plinkcmd = ['plink']

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
        plink_fam = os.path.join(self.scratch, 'plink_multi.fam')

        if not os.path.exists(plink_bed):
            raise IOError('Plink bed doesn\'t exist')
        else:
            print("----phenotypes: " + phenotypes['multi']['id'] + " ----")
            print("----plink binary generated----")
            print("----" + plink_bed + "----")

        if not os.path.exists(plink_bim):
            raise IOError('Plink bim doesn\'t exist')
        else:
            print("----" + plink_bim + "----")

        if not os.path.exists(plink_fam):
            raise IOError('Plink fam doesn\'t exist')
        else:
            os.remove(plink_fam)
            if os.path.exists(phenotypes['multi']['file']):
                shutil.copyfile(phenotypes['multi']['file'], plink_fam)
                phenotypes['multi']['file'] = plink_fam
                print("----" + plink_fam + "----\n")
            else:
                raise ValueError('Constructed multi-phenotype fam file was not created previously.')

        return phenotypes

    def run_assoc_exp(self, params):
        if params['model'] is 0:
            # univariate analysis
            phenotype = self._mk_phenos_from_trait_matrix_uni(params['trait_matrix'])
            plink_prefixes = self._mk_plink_bin_uni(phenotype)
        else:
            # mutlivariate analysis
            phenotype = self._mk_phenos_from_trait_matrix_multi(params['trait_matrix'])
            plink_prefixes = self._mk_plink_bin_multi(phenotype)

        kinmatrix = self._mk_centered_kinship(plink_prefixes)

        for x in range(0, len(kinmatrix)):
            assoc_base_file_prefix = 'gemma_assoc'
            assoc_args = ['-bfile', kinmatrix[x]['plink'], '-k', kinmatrix[x]['kinship'], '-lmm', '4', '-o',
                          assoc_base_file_prefix + kinmatrix[x]['id'] + ".assoc.txt"]
            assoc_cmd = ['gemma']

            for arg in assoc_args:
                assoc_cmd.append(arg)

            assoc_results = kinmatrix

            try:
                proc = subprocess.Popen(assoc_cmd, cwd=self.scratch)
                proc.wait()

                if proc.returncode is -2:
                    # brent error
                    newkinship = self._mk_standardized_kinship(kinmatrix[x])
                    new_assoc_cmd = ['gemma','-bfile', kinmatrix[x]['plink'], '-k', newkinship, '-lmm', '4', '-o',
                                     assoc_base_file_prefix + kinmatrix[x]['id'] + ".assoc.txt"]

                    try:
                        newproc = subprocess.Popen(new_assoc_cmd, cwd=self.scratch)
                        newproc.wait()
                    except Exception as e:
                        exit(e)

                    if not newproc.returncode is -2:
                        # brent error
                        assoc_results[x]['gemma'] = os.path.join(self.scratch, 'output', assoc_base_file_prefix +
                                                                 kinmatrix[x]['id'] + ".assoc.txt")
                    else:
                        assoc_results[x]['gemma'] = ''
                else:
                    assoc_results[x]['gemma'] = os.path.join(self.scratch, 'output',assoc_base_file_prefix +
                                                             kinmatrix[x]['id'] + ".assoc.txt")
            except Exception as e:
                exit(e)

            failed_phenos = []

            if not os.path.exists(assoc_results[x]['gemma']):
                #exit("GEMMA results file does not exist: "+assoc_results[x]['gemma'])
                failed_phenos.append(assoc_results[x]['id'])
            else:
                print("--- gemma results generated: " + assoc_results[x]['gemma'] + "--- \n")

        return assoc_results
