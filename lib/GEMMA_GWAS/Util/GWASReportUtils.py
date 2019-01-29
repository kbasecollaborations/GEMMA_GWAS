import os
import json

"""
    GWASReportUtils:
    
    
"""

class GWASReportUtils:
    def __init__(self, config):
        self.scratch = config["scratch"]
        self.ctx = config['ctx'];
        self.callback_url = config["SDK_CALLBACK_URL"]
        self.htmldir = os.path.join('/kb/module/', '/lib/GEMMA_GWAS/Util/Report/mhplot/')

        print("\n\nHTMLDIR:\n\n"+str(self.htmldir)+"\n\n")
        """
            Might need these later down the road:

            self.ws_client = workspaceService(config["workspace-url"])
            self.gfu = GenomeFileUtil(self.callback_url)
            self.au = AssemblyUtil(self.callback_url)
            self.dfu = DataFileUtil(self.callback_url)
            self.genome_api = GenomeAnnotationAPI(self.callback_url)
            self.vcf_file_ref = variations['variation_file_reference'] # Shock ID? 
        """