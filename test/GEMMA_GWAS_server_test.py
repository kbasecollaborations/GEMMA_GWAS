# -*- coding: utf-8 -*-
import os
import time
import unittest
import subprocess
from configparser import ConfigParser

from GEMMA_GWAS.GEMMA_GWASImpl import GEMMA_GWAS
from GEMMA_GWAS.GEMMA_GWASServer import MethodContext
from GEMMA_GWAS.authclient import KBaseAuth as _KBaseAuth

from installed_clients.WorkspaceClient import Workspace

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
        ret = cls.wsClient.create_workspace({'workspace': cls.wsName})  # noqa

    @classmethod
    def tearDownClass(cls):
        if hasattr(cls, 'wsName'):
            cls.wsClient.delete_workspace({'workspace': cls.wsName})
            print('Test workspace was deleted')

    def getContext(self):
        return self.__class__.ctx

    # NOTE: According to Python unittest naming rules test method names should start from 'test'. # noqa
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

    def test_dummy_data(self):
        """
            Test if A. thaliana GWAS test data exists and is readable

        """
        test_data_dir = '/kb/deps/testdata/'
        geno_test_data_zip = 'AtPolyDB.zip'
        pheno_test_data = 'LFC_phenotype'

        self.assertTrue(os.path.isdir(test_data_dir))
        self.assertTrue(os.path.isfile(os.path.join(test_data_dir, geno_test_data_zip)))
        self.assertTrue(os.path.isfile(os.path.join(test_data_dir, pheno_test_data)))

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
            'workspace_name': 'pranjan77:narrative_1549050842078',
            'variation': '24237/115/1',
            'trait_matrix': '24237/125/1',
            'model': 0
        })
        
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
    """
    def test_GEMMA_GWAS_twenty_uni(self):
        ret = self.serviceImpl.run_gemma_association(self.getContext(), {
            # Twenty trait testing narrative
            'workspace_name': 'rmr:narrative_1550791707920',
            'variation': '25488/26/1',
            'trait_matrix': '25488/36/1',
            'model': 0
        })
