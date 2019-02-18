import subprocess
import os
from pprint import pprint as pp

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace

class AssociationUtils:
    def __init__(self, config, varfiles):
        self.dfu = DataFileUtil(config['SDK_CALLBACK_URL'])
        self.wsc = Workspace(config['SDK_CALLBACK_URL'])
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
            if col_attr_map_ids[id_dict][1]:
                iids[col_attr_map_ids[id_dict][0]] = col_attr_map_ids[id_dict][1]
            else:
                iids[col_attr_map_ids[id_dict][0]] = "NA"

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

            pheno_tsv_headers = "FID\tIID"

            # TODO: Account for multiple phenotypes
            for id in phenotypeids:
                pheno_tsv_headers += "\t"+id.upper()

            pheno_tsv_headers += '\n'

            with open(phenotype_file_path, 'w') as f:
                f.write(pheno_tsv_headers)

                for x in range(0, len(fids)):
                    # TODO: Account for multiple phenotypes
                    f.write(fids[x]+"\t"+iids[fids[x]]+"\t"+str(phenotypevals[x])+"\n")
                f.close()
            return phenotype_file_path
        else:
            raise ValueError('Cannot write data to VCF; invalid WS type (' + ws_type +
                             ').  Supported types is KBaseMatrices.TraitMatrix')

    def _mk_kinship(self):
        self.local_kinship_prefix = 'test_kinship'
        kin_cmd = ['gemma', '-bfile', self.local_plink_prefix, '-gk', '1', '-o', self.local_kinship_prefix]

        try:
            proc = subprocess.check_output(kin_cmd, cwd=self.local_data_dir)
        except OSError as e:
            exit(e)
        except ValueError as e:
            exit(e)

        self.local_kinship_file = os.path.join(self.local_data_dir,'output',self.local_kinship_prefix+'.cXX.txt')

        if not os.path.exists(self.local_kinship_file):
            exit("Kinship file does not exist: "+self.local_kinship_file)

        return self.local_kinship_file

    def _mk_plink_bin(self, phenotype):
        plink_prefix = 'plink_variation'
        """
            plink flags for .ped .map and phenotype tsv into plink
            plinkvars = ['--make-bed','--ped',self.local_ped_file,'--map',self.local_map_file,'--pheno',self.local_pheno_file,'--allow-no-sex','--chr','1','--out',self.local_plink_prefix]
        """
        plinkvars = ['--make-bed', '--vcf', self.varfile, '--pheno',
                     phenotype, '--allow-no-sex', '--out', plink_prefix]
        plinkcmd = ['plink']

        for arg in plinkvars:
            plinkcmd.append(arg)

        try:
            proc = subprocess.Popen(plinkcmd, cwd=self.scratch)
            proc.wait()
        except:
            exit(e)

        plink_bed = os.path.join(self.scratch, plink_prefix+'.bed')
        plink_bam = os.path.join(self.scratch, plink_prefix+'.bam')
        plink_fam = os.path.join(self.scratch, plink_prefix+'.fam')

        if not os.path.exists(plink_bed) or not os.path.exists(plink_bam) or not os.path.exists(plink_fam):
            #print("----plink binary generated----")
            #print("----"+self.local_plink_bed+"----\n\n")
            #print("----" + self.local_plink_bam + "----\n\n")
            #print("----" + self.local_plink_fam + "----\n\n")
            raise IOError("Plink files do not exist!")

        return [plink_bed, plink_bam, plink_fam]

    def run_assoc_exp(self, trait_matrix_ref):
        phenotype = self._mk_phenotype_from_trait_matrix(trait_matrix_ref)
        plink_prefix = self._mk_plink_bin(phenotype)
        kinmatrix = self._mk_kinship()

        self.local_assoc_results_file_prefix = 'gemma_assoc'
        assoc_args = ['-bfile', plink_prefix, '-k', kinmatrix, '-lmm', '4', '-o', self.local_assoc_results_file_prefix]
        assoc_cmd = ['gemma']

        for arg in assoc_args:
            assoc_cmd.append(arg)

        try:
            proc = subprocess.check_output(assoc_cmd, cwd=self.local_data_dir)
        except (OSError, ValueError) as e:
            exit(e)

        self.local_assoc_results_file = os.path.join(self.local_data_dir, 'output', self.local_assoc_results_file_prefix+".assoc.txt")

        if not os.path.exists(self.local_assoc_results_file):
            exit("GEMMA results file does not exist: "+self.local_assoc_results_file)

        return self.local_assoc_results_file