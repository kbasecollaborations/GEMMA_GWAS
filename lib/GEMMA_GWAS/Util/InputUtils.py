from installed_clients.DataFileUtilClient import DataFileUtil
from installed_clients.WorkspaceClient import Workspace

class InputUtils:
    def __init__(self, config):
        self.config = config
        self.scratch = config["scratch"]
        self.callback_url = config["SDK_CALLBACK_URL"]
        self.dfu = DataFileUtil(self.callback_url)

    def validate(self, params):
        exit(params)
