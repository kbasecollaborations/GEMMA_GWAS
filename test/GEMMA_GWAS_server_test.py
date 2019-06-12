# -*- coding: utf-8 -*-
import os
import time
import unittest
import subprocess
import shutil
import json
from configparser import ConfigParser
from pprint import pprint as pp

from GEMMA_GWAS.GEMMA_GWASImpl import GEMMA_GWAS
from GEMMA_GWAS.GEMMA_GWASServer import MethodContext
from GEMMA_GWAS.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.snp2geneClient import snp2gene

class GEMMA_GWASTest(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        token = os.environ.get('KB_AUTH_TOKEN', None)
        config_file = os.environ.get('KB_DEPLOYMENT_CONFIG', None)
        cls.cfg = {}
        config = ConfigParser()
        config.read(config_file)
        for nameval in config.items('GEMMA_GWAS'):
            cls.cfg[nameval[0]] = nameval[1]
        # Getting username from Auth profile for token
        authServiceUrl = cls.cfg['auth-service-url']
        auth_client = _KBaseAuth(authServiceUrl)
        user_id = auth_client.get_user(token)
        # WARNING: don't call any logging methods on the context object,
        # it'll result in a NoneType error
        cls.ctx = MethodContext(None)
        cls.ctx.update({'token': token,
                        'user_id': user_id,
                        'provenance': [
                            {'service': 'GEMMA_GWAS',
                             'method': 'please_never_use_it_in_production',
                             'method_params': []
                             }],
                        'authenticated': 1})
        cls.wsURL = cls.cfg['workspace-url']
        cls.wsClient = Workspace(cls.wsURL)
        cls.serviceImpl = GEMMA_GWAS(cls.cfg)
        cls.scratch = cls.cfg['scratch']
        cls.callback_url = os.environ['SDK_CALLBACK_URL']
        suffix = int(time.time() * 1000)
        cls.wsName = "test_ContigFilter_" + str(suffix)
        cls.snp2gene = snp2gene(cls.callback_url)

        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getContext(self):
        return self.__class__.ctx

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
    """
    def test_plink_install(self):
        FNULL = open(os.devnull, 'w')
        plink = subprocess.check_call(['plink', '--help'],stdout=FNULL,stderr=subprocess.STDOUT)
        # if exit code 0, then binary returned output and did not fault
        # otherwise CalledProcessError is raised
        self.assertEqual(plink, 0)
        FNULL.close()

    def test_gemma_install(self):
        FNULL = open(os.devnull, 'w')
        gemma = subprocess.check_call(['gemma', '-h'], stdout=FNULL, stderr=subprocess.STDOUT)
        # if exit code 0, then binary returned output and did not fault
        # otherwise CalledProcessError is raised
        self.assertEqual(gemma, 0)
        FNULL.close()
    
    def test_SNP2GENE(self):
        shutil.copy('/kb/module/test/sample_data/snpdata-example.tsv', os.path.join(self.scratch, 'snpdata-example.tsv'))
        gwas_results = os.path.join(self.scratch, 'snpdata-example.tsv')
        params = {
            'genome_obj': '25404/2/1',
            'gwas_result_file': gwas_results
        }

        self.snp2gene.annotate_gwas_results(params)

    def test_baselengths(self):
        assembly_obj = self.dfu.get_objects({'object_refs': ['26587/11/1']})['data'][0]
        contigs = assembly_obj['data']['contigs']
        contig_ids = list(contigs.keys())
        contig_ids.sort()

        contig_baselengths = {}
        for id in contig_ids:
            contig_baselengths[id] = contigs[id]['length']

        prev_len = 0

        pp(contig_baselengths)
    """
    def test_GEMMA_GWAS_single_uni(self):
        # Prepare test objects in workspace if needed using
        # self.getWsClient().save_objects({'workspace': self.getWsName(),
        #                                  'objects': []})
        #
        # Run your method by
        # ret = self.getImpl().your_method(self.getContext(), parameters...)
        #
        # Check returned data with
        # self.assertEqual(ret[...], ...) or other unittest methods

        ret = self.serviceImpl.run_gemma_association(self.getContext(), {
            # Single trait testing narrative
            # Univariate analysis
            #'workspace_name': 'rmr:narrative_1552612680986',
            #'variation': '26322/20/1',
            #'trait_matrix': '26322/28/1',
            #'model': 0

            'workspace_name': 'rmr:narrative_1552612680986',
            'variation': '/kb/module/test/sample_data/vcf/LFC_arabidopsis.vcf',
            'trait_matrix': '26322/28/1',
            'model': 0
        })

    """
    def test_GEMMA_GWAS_two_uni(self):
        ret = self.serviceImpl.run_gemma_association(self.getContext(), {
            # Two trait testing narrative
            'workspace_name': 'rmr:narrative_1550791707920',
            'variation': '25488/26/1',
            'trait_matrix': '25488/6/1',
            'model': 0
        })

    def test_GEMMA_GWAS_two_multi(self):
        ret = self.serviceImpl.run_gemma_association(self.getContext(), {
            # Two trait testing narrative
            'workspace_name': 'rmr:narrative_1550791707920',
            'variation': '25488/26/1',
            'trait_matrix': '25488/6/1',
            'model': 1
        })

    def test_GEMMA_GWAS_twenty_uni(self):
        ret = self.serviceImpl.run_gemma_association(self.getContext(), {
            # Twenty trait testing narrative
            'workspace_name': 'rmr:narrative_1550791707920',
            # Small variation
            # 'variation': '25488/26/1',
            # Large variation
            'variation': '25488/28/1',
            'trait_matrix': '25488/36/1',
            'model': 0
        })
    """
