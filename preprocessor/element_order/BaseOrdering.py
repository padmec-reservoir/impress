class BaseOrdering(object):
    """
    Base class that defines the basic methods and
    attributes for element sorting. It should NOT
    be instantiated.
    """
    
    def __init__(self, mesh_entities, interface_dim):
        """
        Base constructor.

        Parameters
        ----------
        mesh_entities: MeshEntities
            A MeshEntities instance from IMPRESS containing
            the dimension of the elements to be ordered.
        interface_dim: string
            The dimension of the interface shared between elements.
            Admissible values: "nodes", "edges", "faces"
        """

        # An IMPRESS MeshEntities object.
        self.mesh_entities = mesh_entities

        # The dimension of the elements to be ordered.
        self.target_dim = self.mesh_entities.entity_type

        # The dimension of the interface to check for neighbors.
        self.interface_dim = interface_dim
    
    def sort_elements(self):
        raise NotImplementedError
