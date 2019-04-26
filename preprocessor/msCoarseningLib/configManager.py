#import pdb
import configparser as cp
import yaml

class initMoabVariable(object):
    def __init__(self, file_path="variable_settings.yml"):

    pass

class partitionManager(object):
    def __init__(self, file_path = "test_file3.card"):
        pass
        with open(file_path) as file: # Use file to refer to the file object
           data = file.read()
           m = yaml.safe_load(data)
        self.create_tree(m)

        # partition tree
        self.tree = m
        # default config file for each layer
        self.default_config_layer = []


    def create_tree(self,m):
        pass
        print(m)


    def create_root_set(self, config_file):
        self.tree['root'] = {}
        self.tree['root']['config'] = config_file
        self.default_config_layer[0] = config_file


    def add_layer(self):
        pass

    def add_coarse_volume(self):
        pass

        #return f(num) # Prove that function definition has completed
    #
    # def __call__(self,num):
    #     print("inside my_decorator.__call__()")
    #     return self.fun(num)
    #
    # def __getitem__(self, item):
    #     print("inside my_decorator.__call__()")
    #     return item


def readConfig(configInput = "msCoarse.ini"):
    configInput = 'lololita'
    configFile = cp.ConfigParser()
    try:
        configFile.read(configInput)
    except:
        print("Não foi possivel ler o arquivo: "+configInput)
    return configFile


    #pdb.set_trace()
