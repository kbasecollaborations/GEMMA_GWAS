import subprocess
import os

class AssociationUtils:
    def __init__(self, config, files):
        self.scratch = config["scratch"]
        self.local_ped_file = files[0]
        self.local_map_file = files[1]
        self.local_pheno_file = files[2]

        """
            Might need these later down the road:

            self.ws_client = workspaceService(config["workspace-url"])
            self.gfu = GenomeFileUtil(self.callback_url)
            self.au = AssemblyUtil(self.callback_url)
            self.dfu = DataFileUtil(self.callback_url)
            self.genome_api = GenomeAnnotationAPI(self.callback_url)
        """

    def _local_mk_kinship(self):
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

    def _local_get_var_file_list(self):
        return [self.local_ped_file, self.local_map_file, self.local_pheno_file]

    def _local_mk_plink_bin(self):
        self.local_plink_prefix = 'test_plink_bin'
        plinkvars = ['--make-bed','--ped',self.local_ped_file,'--map',self.local_map_file,'--pheno',self.local_pheno_file,'--allow-no-sex','--chr','1','--out',self.local_plink_prefix]
        plinkcmd = ['plink']

        for arg in plinkvars:
            plinkcmd.append(arg)

        try:
            proc = subprocess.check_output(plinkcmd, cwd=self.local_data_dir)
        except OSError as e:
            exit(e)
        except ValueError as e:
            exit(e)

        self.local_plink_bed = os.path.join(self.local_data_dir, self.local_plink_prefix+'.bed')
        self.local_plink_bam = os.path.join(self.local_data_dir, self.local_plink_prefix+'.bed')
        self.local_plink_fam = os.path.join(self.local_data_dir, self.local_plink_prefix+'.bed')

        if not os.path.exists(self.local_plink_bed) or not os.path.exists(self.local_plink_bam) or not os.path.exists(self.local_plink_fam):
            #print("----plink binary generated----")
            #print("----"+self.local_plink_bed+"----\n\n")
            #print("----" + self.local_plink_bam + "----\n\n")
            #print("----" + self.local_plink_fam + "----\n\n")
            exit("Plink files do not exist!")

        return self.local_plink_prefix

    def local_run_assoc(self):
        self.local_data_dir = '/kb/deps/testdata/'
        var_files = self._local_get_var_file_list()
        plink_prefix = self._local_mk_plink_bin()
        kinmatrix = self._local_mk_kinship()

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