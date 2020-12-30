class BaseOrdering(object):
    """
    Base class that defines the basic methods and
    attributes for element sorting. It should NOT
    be instantiated.
    """
    
    def __init__(self, mesh_entities, target_dim, interface_dim):
        # A IMPRESS MeshEntities object.
        self.mesh_entities = mesh_entities

        # The dimension of the elements to be ordered.
        self.target_dim = target_dim

        # The dimension of the interface to check for neighbors.
        self.interface_dim = interface_dim
    
    def sort_elements(self):
        raise NotImplementedError
