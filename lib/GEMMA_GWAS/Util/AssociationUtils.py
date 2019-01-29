import subprocess
import os

class AssociationUtils:
    def __init__(self, config, files):
        self.scratch = config["scratch"]
        self.ctx = config['ctx']
        self.callback_url = config["SDK_CALLBACK_URL"]
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

    def _local_mk_return_kinship(self):
        kinship = 'kinship return'
        return kinship

    def _local_get_var_file_list(self):
        return [self.local_ped_file, self.local_map_file, self.local_pheno_file]

    def _local_mk_plink_bin(self):
        self.local_plink_file = os.path.join(self.local_data_dir, 'test_bin')
        plinkvars = ['--ped '+self.local_ped_file, '--map '+self.local_map_file, '--pheno '+self.local_pheno_file, '--output '+self.local_plink_file]
        plinkcmd = ['plink']

        for arg in plinkvars:
            plinkcmd.append(arg)

        try:
            proc = subprocess.Popen(plinkcmd)
        except OSError as e:
            exit(e)
        except ValueError as e:
            exit(e)

        if os.path.exists(self.local_plink_file):
            print("----plink binary generated----")
            print("----"+self.local_plink_file+"----\n\n")

        return self.local_plink_file

    def local_run_assoc(self):
        self.local_data_dir = '/kb/deps/testdata/'
        var_files = self._local_get_var_file_list()
        plink_bin = self._local_mk_plink_bin()
        kinmatrix = self._local_mk_return_kinship()

        assoc = 'assoc returns'
        return assoc