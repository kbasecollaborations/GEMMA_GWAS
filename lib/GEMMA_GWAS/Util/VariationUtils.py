import os
import vcf

"""
    VariationUtils:

    This will be the class that is associated with handling, parsing, and manipulating the Variant Call Format (VCF)
    file.

    Once the variation data has been input we can then pass to association based utils for plink binary generation,
    kinship generation, and ultimate association data calculations.

    Ideally this would be a good file to use KBase Parallel and Caching Service with as the VCF file can be rather
    large.

    Some ideas that have been dicussed include chopping up the VCF into multiple files, i.e. for each 10k SNPs, and
    processing in parallel
"""

class VariationUtils:
    def __init__(self, config):
        self.scratch = config["scratch"]
        self.ctx = config['ctx'];
        self.callback_url = config["SDK_CALLBACK_URL"]

        """
            Might need these later down the road:

            self.ws_client = workspaceService(config["workspace-url"])
            self.gfu = GenomeFileUtil(self.callback_url)
            self.au = AssemblyUtil(self.callback_url)
            self.dfu = DataFileUtil(self.callback_url)
            self.genome_api = GenomeAnnotationAPI(self.callback_url)
            self.vcf_file_ref = variations['variation_file_reference'] # Shock ID? 
        """

        def _get_local_vcf():
            """
                For local testing we will be using a test vcf file

                :param key: Key to search in Params
                :returns: string filepath
                :raises ValueError: raises an exception if file path doesn't exist
            """

            local_vcf_file = '/kb/testdata/test.vcf'

            try:
                if not os.path.exists(local_vcf_file):
                    raise ValueError('Local vcf file does not exist')
            except ValueError as e:
                exit(str(e))

            return local_vcf_file

        def get_vcf_data(self):
            """
                Get VCF file and return parsed object

                :return: vcf.Reader object: https://pyvcf.readthedocs.io/en/latest/API.html#vcf-reader
            """

            return self.parse_vcf_data(_get_local_vcf())

        def parse_vcf_data(self, file_path):
            """
            Get value of key after checking for its existence
                :param params: Params dictionary haystack
                :param key: Key to search in Params
                :return: Parameter Value
                :raises ValueError: raises an exception if the key doesn"t exist
            """
            # TODO: Validate parsing, i.e. vcf.Reader obj
            return vcf.Reader(open(file_path, 'r'))
