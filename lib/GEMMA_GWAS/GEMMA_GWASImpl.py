# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os

from installed_clients.KBaseReportClient import KBaseReport
from GEMMA_GWAS.Util.VariationUtils import VariationUtils
from GEMMA_GWAS.Util.AssociationUtils import AssociationUtils

#END_HEADER


class GEMMA_GWAS:
    '''
    Module Name:
    GEMMA_GWAS

    Module Description:
    A KBase module: GEMMA_GWAS
    '''

    ######## WARNING FOR GEVENT USERS ####### noqa
    # Since asynchronous IO can lead to methods - even the same method -
    # interrupting each other, you must be *very* careful when using global
    # state. A method could easily clobber the state set by another while
    # the latter method is running.
    ######################################### noqa
    VERSION = "0.0.1"
    GIT_URL = ""
    GIT_COMMIT_HASH = ""

    #BEGIN_CLASS_HEADER
    #END_CLASS_HEADER

    # config contains contents of config file in a hash or None if it couldn't
    # be found
    def __init__(self, config):
        #BEGIN_CONSTRUCTOR
        self.config = config
        self.config['SDK_CALLBACK_URL'] = os.environ['SDK_CALLBACK_URL']
        self.config['KB_AUTH_TOKEN'] = os.environ['KB_AUTH_TOKEN']
        self.config['test_data_dir'] = os.path.abspath('/kb/testdata')
        # self.config['scratch'] is the tmp directory

        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

        #END_CONSTRUCTOR
        pass


    def run_GEMMA_GWAS(self, ctx, params):
        """
        :param params: instance of type "GemmaGwasInput" -> structure:
           parameter "output_ws" of String, parameter "Variation" of type
           "obj_ref" (An X/Y/Z style reference), parameter "Associations" of
           type "obj_ref" (An X/Y/Z style reference)
        :returns: instance of type "GwasResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_GEMMA_GWAS
        self.config['ctx'] = ctx

        variations = VariationUtils(self.config)
        associations = AssociationUtils(self.config)

        report = KBaseReport(self.config['SDK_CALLBACK_URL'])
        report_msg = "The variation object: "+params['Variation']+"\nThe association object:"+params['Associations']+"\n"
        report_info = report.create({'report': {'objects_created':[],
                                                'text_message': report_msg},
                                                'workspace_name': params['output_ws']})
        output = {
            'report_name': report_info['name'],
            'report_ref': report_info['ref'],
            'ws': params['output_ws']
        }
        #END run_GEMMA_GWAS

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_GEMMA_GWAS return value ' +
                             'output is not type dict as required.')
        # return the results
        return [output]
    def status(self, ctx):
        #BEGIN_STATUS
        returnVal = {'state': "OK",
                     'message': "",
                     'version': self.VERSION,
                     'git_url': self.GIT_URL,
                     'git_commit_hash': self.GIT_COMMIT_HASH}
        #END_STATUS
        return [returnVal]
