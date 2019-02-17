import os
import json
import csv
import shutil
import logging

"""
    GWASReportUtils:
    
    
"""

class GWASReportUtils:
    def __init__(self, config):
        self.config = config
        self.scratch = config["scratch"]
        self.ctx = config['ctx'];
        self.callback_url = config["SDK_CALLBACK_URL"]
        #self.htmldir = os.path.join('/kb/module', 'lib/GEMMA_GWAS/Util/Report/mhplot')
        self.htmldir = os.path.join(self.scratch,'mhplot')

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
        gemma_results = open(assoc_file, 'r', newline='')
        tsv_unfiltered = csv.reader(gemma_results, delimiter='\t')
        next(tsv_unfiltered, None) # skip old csv headers

        tsv_sorted = sorted(tsv_unfiltered, key=lambda row: float(row[13]))

        tsv_filtered_headers = "SNP\tCHR\tBP\tP\n"
        self.local_filtered_tsv_file = os.path.join('/kb/module', 'lib/GEMMA_GWAS/Util/Report/mhplot','snpdata.tsv')
        assoc_entry_limit = 5000

        with open(self.local_filtered_tsv_file,'w') as tsv_filtered:
            tsv_filtered.write(tsv_filtered_headers)

            k = 0
            for snp in tsv_sorted:
                if k < assoc_entry_limit:
                    tsv_filtered.write(snp[1]+"\t"+snp[0]+"\t"+snp[2]+"\t"+snp[13]+"\n")
                    k += 1

            tsv_filtered.close()

        return self.local_filtered_tsv_file

    def _copy_html_to_scratch(self):
        dst = self.htmldir
        src = os.path.join('/kb/module', 'lib/GEMMA_GWAS/Util/Report/mhplot')

        try:
            shutil.copytree(src, dst)
        except OSError as why:
            exit('src: '+src+'\ndst: '+dst+'\nerror: '+str(why))

    def mk_html_report(self, assoc_file):
        self.local_assoc_results_file = assoc_file
        self._filter_local_assoc_results(assoc_file)
        self._copy_html_to_scratch()

        logging.info("\n\n\nfiltered:\n")
        os.system("wc -l "+os.path.join(self.htmldir, 'snpdata.tsv'))
        logging.info("\n\n\nunfiltered:\n")
        os.system("wc -l " + assoc_file)
        logging.info("\n\n")

        html_return = {
            'path': self.htmldir,
            'name': "index.html",
            'description': 'Manhattan plot of GEMMA GWAS association tests'
        }

        return html_return