import logging

from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace


class InputUtils:
    def __init__(self, config):
        self.config = config
        self.scratch = config["scratch"]
        self.callback_url = config["SDK_CALLBACK_URL"]
        self.dfu = DataFileUtil(self.callback_url)

    def validate(self, params):
        if params['model'] is 0:
            # univariate linear mixed model
            #TODO This is the place where trait data format validation should be done
            return True

        """
            
            trait_matrix = self.dfu.get_objects({'object_refs': [params['trait_matrix']]})['data'][0]['data']['data']
            phenos = trait_matrix["row_ids"]

            if len(phenos) <= 1:
                logging.info('Running a univariate analysis on the '+phenos[0]+' phenotype.')
                return True
            else:
                logmsg = 'Running a univariate analysis on '+str(len(phenos))+' phenotypes:\n'
                for pheno in phenos:
                    logmsg += pheno + "\n"
                logmsg += "\n"

                logging.info(logmsg)

                return True
        elif params['model'] is 1:
            # multivariate linear mixed model
            trait_matrix = self.dfu.get_objects({'object_refs': [params['trait_matrix']]})['data'][0]['data']['data']
            phenos = trait_matrix["row_ids"]

            if len(phenos) <= 1:
                raise ValueError('Cannot run multivariate analysis on a single/empty phenotype.')
            else:
                logmsg = 'Running a multivariate analysis on '+str(len(phenos))+' phenotypes:'
                for pheno in phenos:
                    logmsg += pheno + '\n'
                logmsg += '\n'

                logging.info(logmsg)

                return True
        else:
            logging.error("Model is " + params['model'] + ", of type: " + str(type(params['model'])))
            raise ValueError('Model selected is not either multivariate or univariate. Contact Module maintainer. ')
            
        """
