import pickle
from pymoab import core, types, rng, topo_util, skinner, tag

class IMPRESSPickler(pickle.Pickler):
    def __init__(self, file):
        super().__init__(file)
        self.file_name = file.name
        self.core_written = False
    
    def persistent_id(self, obj):
        if isinstance(obj, core.Core):
            if not self.core_written:
                obj.write_file(f'{self.file_name}.h5m')
                self.core_written = True
            return ("core",)
        elif isinstance(obj, topo_util.MeshTopoUtil):
            return ("mtu",)
        elif isinstance(obj, skinner.Skinner):
            return ("skinner",)
        elif isinstance(obj, rng.Range):
            array = obj.get_array()
            return ("range", array)
        elif isinstance(obj, tag.Tag):
            name = obj.get_name()
            return ("tag", name)
        return None

class IMPRESSUnpickler(pickle.Unpickler):
    def __init__(self, file):
        super().__init__(file)
        self.file_name = file.name
        self.mb = None
        self.mtu = None
        self.skinner = None
    
    def persistent_load(self, pid):
        type_tag = pid[0]
        if type_tag == "core":
            if self.mb is None:
                mb = core.Core()
                mb.load_file(f'{self.file_name}.h5m')
                self.mb = mb
            return self.mb
        elif type_tag == "mtu":
            if self.mb is None:
                raise pickle.UnpicklingError("core is missing")
            if self.mtu is None:
                self.mtu = topo_util.MeshTopoUtil(self.mb)
            return self.mtu
        elif type_tag == "skinner":
            if self.mb is None:
                raise pickle.UnpicklingError("core is missing")
            if self.skinner is None:
                self.skinner = skinner.Skinner(self.mb)
            return self.skinner
        elif type_tag == "range":
            r = rng.Range(pid[1])
            return r
        elif type_tag == "tag":
            if self.mb is None:
                raise pickle.UnpicklingError("core is missing")
            return self.mb.tag_get_handle(pid[1])
        
        raise pickle.UnpicklingError(f'unsupported persistent object {pid}')
