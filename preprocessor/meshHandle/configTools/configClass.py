import yaml
import numpy as np


class coarseningInit(object):
    def __init__(self, empty=False, file_path=None):
        if file_path is None:
            file_path = 'input_cards/coarsening.yml'
        if empty is False:
            self.read(file_path)
        else:
            self.tree = {'Scheme': '', 'Simple': '', 'Smart': ''}

    def scheme(self, scheme='smart'):
        self.tree['Scheme'] = scheme

    def simple(self, nx=4, ny=4, nz=4):
        self.tree['Scheme'] = 'simple'
        self.tree['Simple'] = {'nx': nx, 'ny': ny, 'nz': nz}

    def smart(self, path='default', file='mesh.h5m'):
        self.tree['Scheme'] = 'smart'
        self.tree['Smart'] = {'path': path, 'file': file}

    def read(self, file_path):
        with open(file_path) as file:  # Use file to refer to the file object
            data = file.read()
            self.tree = yaml.safe_load(data)

    # def write(self, file_path):
    #     with open(file_path, 'w') as f:
    #         data = yaml.dump(self.tree, f)


class variableInit(object):
    def __init__(self, empty=False, file_path=None):
        self.level_tree = {}
        if file_path is None:
            file_path = 'input_cards/variable_input.yml'
        if empty is False:
            self.read(file_path)
        else:
            self.tree = np.array([[], [], [], [], [], [], []]).T

    def get_var(self, level=0, coarse_num=0):
        ref = self.tree[:, 0].astype('int') == level
        level_entries = self.tree[ref]
        # self.tree = np.delete(self.tree, np.where(ref)[0], 0)
        return self.create_command(level_entries, level, coarse_num)

    def create_command(self, var_list, level, coarse_num):
        command_list = []
        cmd_pattern = 'self.{0} = MoabVariable(self.core, data_size={1}, var_type=\'{2}\', data_format=\'{3}\', name_tag=\'{0}\', data_density=\'{4}\')'
        cmd_pattern_level = 'self.{0} = MoabVariableMS(self.core, data_size={1}, var_type=\'{2}\',level={5}, coarse_num={6}, data_format=\'{3}\', name_tag=\'{0}\', data_density=\'{4}\')'
        cmd_pattern_list = 'self.var_handle_list.append(self.{0})'
        density = ['dense', 'sparse']
        if level == 0:
            for var in var_list:
                sparse = density[bool(var[-2])]
                command_list.append(cmd_pattern.format(var[1], var[3], var[2], var[-3], sparse))
                command_list.append(cmd_pattern_list.format(var[1]))
        else:
            for var in var_list:
                sparse = density[bool(var[-2])]
                command_list.append(cmd_pattern_level.format(var[1], var[3], var[2], var[4], sparse, level, coarse_num))
                command_list.append(cmd_pattern_list.format(var[1]))
        return command_list

    def add_entry(self, name, type, data_size, data_type, sparse=False,
                  level=0):
        entry = [level, name, type, data_size, data_type, sparse, []]
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
                            tree[var]['data type'], bool(tree[var]['sparse']),
                            list()])
        return np.array(all_var)
