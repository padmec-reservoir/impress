from .BaseOrdering import BaseOrdering 
import numpy as np

class MPFAD2DOrdering(BaseOrdering):
    """
    Implementation of the 2D MPFA-D entity ordering.
    """

    def __init__(self, mesh_entities, interface_dim):
        """Constructor"""
        super().__init__(mesh_entities, interface_dim)

        if self.mesh_entities not in ("nodes", "edges", "faces"):
            raise ValueError("Entities must be nodes, edges or faces")
    
    
    def sort_elements(self, elements):
        """
        Sort elements according to the MPFA-D ordering.
        """
        visited_entities = [elements[0]]
        elements_set = set(elements)
        num_elements = elements.shape[0]
        while len(visited_entities) < num_elements:
            curr_entity = visited_entities[-1]
            interface_neighbors = self.mesh_entities.bridge_adjacencies(curr_entity, 
                                                                        self.interface_dim, self.target_dim)
            unvisited_entities = elements_set - set(visited_entities)
            next_entities = unvisited_entities & set(interface_neighbors)
            visited_entities.append(next_entities.pop())
        
        ordered_elements = np.array(visited_entities)
        return ordered_elements