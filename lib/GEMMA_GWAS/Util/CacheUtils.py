import requests
import json

from pprint import pprint as pp

class CacheUtils:
    def __init__(self, config, cache_id, service_token):
        self.callback = config['SDK_CALLBACK_URL']
        if not self.callback.endswith('/'):
            self.cacheurl = self.callback + '/v1/cache/'
        else:
            self.cacheurl = self.callback + 'v1/cache/'

        self.service_token = service_token
        self.cache_id = cache_id

    def downld_cache(self, destination):
        headers = {'Content-type': 'application/json', 'Authorization': self.service_token}
        endpoint = self.cacheurl + '/' + self.cache_id
        req_call = requests.get(endpoint, headers=headers, stream=True)




    def delete_cache(cls, cache_id, service_token):
        cls.deletecachecurl = 'curl -X DELETE'+ \
                             '-H "Authorization: ' + service_token + \
                             'https://' + cls.cacheurl + 'v1/cache/' + cache_id
        return cls.deletecachecurl

