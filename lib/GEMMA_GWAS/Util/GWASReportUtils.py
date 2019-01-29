import os
import json
import csv

"""
    GWASReportUtils:
    
    
"""

class GWASReportUtils:
    def __init__(self, config):
        self.scratch = config["scratch"]
        self.ctx = config['ctx'];
        self.callback_url = config["SDK_CALLBACK_URL"]
        self.htmldir = os.path.join('/kb/module', 'lib/GEMMA_GWAS/Util/Report/mhplot/')

        """
            Might need these later down the road:

            self.ws_client = workspaceService(config["workspace-url"])
            self.gfu = GenomeFileUtil(self.callback_url)
            self.au = AssemblyUtil(self.callback_url)
            self.dfu = DataFileUtil(self.callback_url)
            self.genome_api = GenomeAnnotationAPI(self.callback_url)
            self.vcf_file_ref = variations['variation_file_reference'] # Shock ID? 
        """

    def _filter_local_assoc_results(self, assoc_file):
        tsv_unfilered = csv.reader(open(assoc_file, 'r'), delimiter='\t')

        tsv_sorted = sorted(tsv_unfilered, key=lambda row: row[8])

        tsv_filtered_headers = ['SNP', 'CHR', 'BP', 'P']
        self.local_filtered_tsv_file = os.path.join(self.htmldir,'snpdata.tsv')
        tsv_filtered = csv.writer(open(self.local_filtered_tsv_file,'w'),delimiter='\t')
        tsv_filtered.writerow(tsv_filtered_headers)

        assoc_entry_limit = 5000
        k = 0

        for snp in tsv_unfilered:
            if k < assoc_entry_limit:
                tsv_filtered.writerow([row[1],row[0],row[2],row[8]])
            k+=1

        return self.local_filtered_tsv_file

    def mk_html_report(self, assoc_file):
        self.local_assoc_results_file = assoc_file
        self._filter_local_assoc_results(assoc_file)

        return self.htmldir