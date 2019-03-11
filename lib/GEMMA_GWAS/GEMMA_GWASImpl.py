# -*- coding: utf-8 -*-
#BEGIN_HEADER
import logging
import os
import uuid

from installed_clients.KBaseReportClient import KBaseReport
from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.VariationUtilClient import VariationUtil
from GEMMA_GWAS.Util.InputUtils import InputUtils
from GEMMA_GWAS.Util.AssociationUtils import AssociationUtils
from GEMMA_GWAS.Util.GWASReportUtils import GWASReportUtils

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
    VERSION = "0.0.2"
    GIT_URL = "https://github.com/rroutsong/GEMMA_GWAS.git"
    GIT_COMMIT_HASH = "8050c6b68bc8248777501cc97ba4cde31a85aa74"

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
        self.dfu = DataFileUtil(self.config['SDK_CALLBACK_URL'])
        # self.config['scratch'] is the tmp directory

        logging.basicConfig(format='%(created)s %(levelname)s: %(message)s',
                            level=logging.INFO)

        #END_CONSTRUCTOR
        pass

    def run_gemma_association(self, ctx, params):
        """
        :param params: instance of type "GemmaGwasInput" -> structure:
           parameter "workspace_name" of String, parameter "trait_matrix" of
           type "trait_ref" (KBase style object reference X/Y/Z to a
           KBaseMatrices.TraitMatrix structure @id ws
           KBaseMatrices.TraitMatrix), parameter "variation" of type
           "var_ref" (KBase style object reference X/Y/Z to a @id ws
           KBaseGwasData.Variations)
        :returns: instance of type "GwasResults" -> structure: parameter
           "report_name" of String, parameter "report_ref" of String
        """
        # ctx is the context object
        # return variables are: output
        #BEGIN run_gemma_association

        if 'variation' not in params:
            raise ValueError('Variation KBase reference not set.')
        if 'model' not in params:
            raise ValueError('GEMMA linear mixed model not selected.')
        if 'trait_matrix' not in params:
            raise ValueError('Trait matrix KBase reference not set.')

        # Validate linear mixed model selection
        InputUtils(self.config).validate(params)

        # Get Variation file
        variations = VariationUtil(self.config['SDK_CALLBACK_URL'])
        variation_info = variations.get_variation_as_vcf({
            'variation_ref' : params['variation'],
            'filename': os.path.join(self.config['scratch'], 'variation.vcf')
        })

        # For faster tests, let's use a local file vcf
        # variation_info = {'path': '/kb/module/work/variation.vcf'}

        # Run association tests
        associations = AssociationUtils(self.config, variation_info['path'])
        assoc_file = associations.run_assoc_exp(params)

        # Generate association report
        assoc_report = GWASReportUtils(self.config)
        report_obj = assoc_report.mk_output(params, assoc_file)

        # Send report
        report_client = KBaseReport(self.config['SDK_CALLBACK_URL'])
        report = report_client.create_extended_report(report_obj)

        output = {
            'report_name': report['name'],
            'report_ref': report['ref'],
            'ws': params['workspace_name']
        }

        #END run_gemma_association

        # At some point might do deeper type checking...
        if not isinstance(output, dict):
            raise ValueError('Method run_gemma_association return value ' +
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
