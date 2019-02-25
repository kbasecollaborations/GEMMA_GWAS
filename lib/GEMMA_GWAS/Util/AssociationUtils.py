import subprocess
import os
import operator
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

    def _mk_phenotypes_from_trait_matrix(self, trait_matrix_ref):
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

    def _mk_kinship(self, plink_prefixes):
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

    def _mk_plink_bin(self, phenotypes):
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

    def run_assoc_exp(self, params):
        phenotype = self._mk_phenotypes_from_trait_matrix(params['trait_matrix'])
        plink_prefixes = self._mk_plink_bin(phenotype)
        kinmatrix = self._mk_kinship(plink_prefixes)

        for x in range(0, len(kinmatrix)):
            assoc_base_file_prefix = 'gemma_assoc'
            assoc_args = ['-bfile', kinmatrix[x]['plink'], '-k', kinmatrix[x]['kinship'], '-lmm', '4', '-o', assoc_base_file_prefix+str(x)]
            assoc_cmd = ['gemma']

            for arg in assoc_args:
                assoc_cmd.append(arg)

            try:
                proc = subprocess.Popen(assoc_cmd, cwd=self.scratch)
                proc.wait()
            except Exception:
                exit(e)

            assoc_results = kinmatrix

            assoc_results[x]['gemma'] = os.path.join(self.scratch, 'output', assoc_base_file_prefix+str(x)+".assoc.txt")

            if not os.path.exists(assoc_results[x]['gemma']):
                exit("GEMMA results file does not exist: "+assoc_results[x]['gemma'])
            else:
                print("--- gemma results generated: " + assoc_results[x]['gemma'] + "--- \n")

        return assoc_results
