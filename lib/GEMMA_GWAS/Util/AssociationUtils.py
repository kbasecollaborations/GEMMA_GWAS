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

    def _mk_phenotype_from_trait_matrix(self, trait_matrix_ref):
        trait_matrix_obj = self.dfu.get_objects({'object_refs': [trait_matrix_ref]})['data'][0]

        ws_type = trait_matrix_obj['info'][2]
        if 'KBaseMatrices.TraitMatrix' in ws_type:
            fids = trait_matrix_obj['data']['data']['col_ids']
            iids = self._get_iids(trait_matrix_obj['data']['col_attributemapping_ref'])
            phenotypevals = trait_matrix_obj['data']['data']['values'][0]
            phenotypeids = trait_matrix_obj['data']['data']['row_ids']
            # obj_name = trait_matrix_obj['info'][1]
            phenotype_file_path = os.path.join(self.scratch, 'phenotype.txt')

            if not len(fids) == len(iids) or not len(iids) == len(phenotypevals):
                raise ValueError('Retrieved family-ids, within-family-ids, and phenotype \
                    length of values mismatch')

            pheno_tsv_headers = "FID IID"

            # TODO: Account for multiple phenotypes
            for id in phenotypeids:
                pheno_tsv_headers += " \""+id.upper()+"\""

            pheno_tsv_headers += '\n'

            with open(phenotype_file_path, 'w') as f:
                f.write(pheno_tsv_headers)

                for x in range(0, len(fids)):
                    # TODO: Account for multiple phenotypes
                    if iids[fids[x]]:
                        f.write(fids[x] + " " + iids[fids[x]] + " " + str(phenotypevals[x]) + "\n")
                    else:
                        f.write(fids[x] + " " + fids[x] + " " + str(phenotypevals[x]) + "\n")

                f.close()
            return phenotype_file_path
        else:
            raise ValueError('Cannot write data to VCF; invalid WS type (' + ws_type +
                             ').  Supported types is KBaseMatrices.TraitMatrix')

    def _mk_kinship(self, plink_prefix):
        kinship_prefix = 'kinship'
        kin_cmd = ['gemma', '-bfile', os.path.join(self.scratch, plink_prefix), '-gk', '1', '-o', kinship_prefix]

        try:
            proc = subprocess.Popen(kin_cmd, cwd=self.scratch)
            proc.wait()
        except Exception:
            exit(e)

        kinship_file = os.path.join(self.scratch, 'output', kinship_prefix+'.cXX.txt')

        if not os.path.exists(kinship_file):
            exit("Kinship file does not exist: "+kinship_file)

        return kinship_file

    def _mk_plink_bin(self, phenotype):
        plink_prefix = 'plink_variation'
        """
            plink flags for .ped .map and phenotype tsv into plink
            plinkvars = ['--make-bed','--ped',self.local_ped_file,'--map',self.local_map_file,'--pheno',
                    self.local_pheno_file,'--allow-no-sex','--chr','1','--out',self.local_plink_prefix]
        """
        plinkvars = ['--make-bed', '--vcf', self.varfile, '--pheno',
                     phenotype, '--allow-no-sex','--allow-extra-chr', '--out', plink_prefix]
        plinkcmd = ['plink']

        for arg in plinkvars:
            plinkcmd.append(arg)

        try:
            proc = subprocess.Popen(plinkcmd, cwd=self.scratch)
            proc.wait()
        except Exception:
            exit(e)

        plink_bed = os.path.join(self.scratch, plink_prefix+'.bed')
        plink_bim = os.path.join(self.scratch, plink_prefix+'.bim')
        plink_fam = os.path.join(self.scratch, plink_prefix+'.fam')

        if not os.path.exists(plink_bed):
            raise IOError('Plink bed doesn\'t exist')
        else:
            print("----plink binary generated----")
            print("----" + plink_bed + "----\n")


        if not os.path.exists(plink_bim):
            raise IOError('Plink bim doesn\'t exist')
        else:
            print("----" + plink_bim + "----\n")

        if not os.path.exists(plink_fam):
            raise IOError('Plink fam doesn\'t exist')
        else:
            print("----" + plink_fam + "----\n")

        return plink_prefix

    def run_assoc_exp(self, params):
        phenotype = self._mk_phenotype_from_trait_matrix(params['trait_matrix'])
        plink_prefix = self._mk_plink_bin(phenotype)
        kinmatrix = self._mk_kinship(plink_prefix)

        assoc_file_prefix = 'gemma_assoc'
        assoc_args = ['-bfile', plink_prefix, '-k', kinmatrix, '-lmm', '4', '-o', assoc_file_prefix]
        assoc_cmd = ['gemma']

        for arg in assoc_args:
            assoc_cmd.append(arg)

        try:
            proc = subprocess.Popen(assoc_cmd, cwd=self.scratch)
            proc.wait()
        except Exception:
            exit(e)

        assoc_results_file = os.path.join(self.scratch, 'output', assoc_file_prefix+".assoc.txt")

        if not os.path.exists(assoc_results_file):
            exit("GEMMA results file does not exist: "+assoc_results_file)
        else:
            print("--- gemma results generated: " + assoc_results_file + "--- \n")

        return assoc_results_file
