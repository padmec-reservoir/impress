import yaml
import numpy as np


class variableInit(object):
    def __init__(self, empty=False, file_path=None):
        if file_path is None:
            file_path = 'input_cards/variable_input.yml'
        if empty is False:
            self.read(file_path)
        else:
            self.tree = np.array([[], [], [], [], [], []]).T

    def get_var(self, level):
        ref = self.tree[:, 0].astype('int') == level
        level_entries = self.tree[ref]
        self.tree = np.delete(self.tree, np.where(ref)[0], 0)

    def create_command(self, var_list):
        command_list = []
        string_pattern = '{0} = MoabVariable(self.core, data_size={1}, var_type=\'{2}\', data_format=\'{3}\', name_tag=\'{0}\', data_density=\'{4}\')'
        density = ['int', 'sparse']
        for var in var_list:
            sparse = density[bool(var[-1])]
            command_list.append(string_pattern.format(var[1], var[3], var[2], var[-2], sparse))
        return command_list

    def add_entry(self, name, type, data_size, data_type, sparse=False,
                  level=0):
        entry = [level, name, type, data_size, data_type, sparse]
        self.tree = np.vstack((self.tree, np.array(entry)))

    def read(self, file_path):
        with open(file_path) as file:
            data = file.read()
            tree = yaml.safe_load(data)
        self.tree = self.list_variable(tree)

    def list_variable(self, tree):
        all_var = []
        for var in list(tree):
            all_var.append([int(tree[var]['level']), var,
                            tree[var]['type'],
                            int(tree[var]['data size']),
                            tree[var]['data type'], bool(tree[var]['sparse'])])
        return np.array(all_var)
