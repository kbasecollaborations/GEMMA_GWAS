class AssociationUtils:
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
        """