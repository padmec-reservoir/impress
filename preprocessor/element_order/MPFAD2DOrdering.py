from .BaseOrdering import BaseOrdering 
from functools import cmp_to_key
import numpy as np

class MPFAD2DOrdering(BaseOrdering):
    """
    Implementation of the 2D MPFA-D entity ordering.
    """

    def __init__(self, mesh_entities, interface_dim):
        """Constructor"""
        super().__init__(mesh_entities, interface_dim)

        if self.mesh_entities.entity_type not in ("nodes", "edges", "faces"):
            raise ValueError("Entities must be nodes, edges or faces")
    
    
    def sort_elements(self, elements, center):
        """
        Sort elements according to the MPFA-D ordering.

        Parameters
        ----------
        elements: Numpy array
            Array containing the indices of the elements to be
            ordered.
        from_element: Numpy array
            Array containing the centroid of the elements distribution.
        Returns
        -------
        Numpy array of ordered elements.
        """
        visited_entities = [elements[0]]
        elements_set = set(elements)
        num_elements = elements.shape[0]
        while len(visited_entities) < num_elements:
            curr_entity = visited_entities[-1]
            interface_neighbors = self.mesh_entities.bridge_adjacencies(curr_entity, 
                                                                        self.interface_dim, 
                                                                        self.target_dim)
            unvisited_entities = elements_set - set(visited_entities)
            next_entities = unvisited_entities & set(interface_neighbors)
            visited_entities.append(next_entities.pop())
        
        ordered_elements = np.array(visited_entities)

        elements_centers = self.mesh_entities.center[ordered_elements][:, 0:2]
        if num_elements > 2:
            # Check the orientation of the elements.
            A, B, C = elements_centers[0], elements_centers[1], elements_centers[2]

            # If clockwise, reverse so it is counterclockwise.
            if (B[0] - A[0])*(C[1] - A[1]) - (C[0] - A[0])*(B[1] - A[1]) < 0:
                ordered_elements = np.flip(ordered_elements)
        elif num_elements == 2:
            A, B = elements_centers[0], elements_centers[1]
            if np.cross(center[0, 0:2] - A, center[0, 0:2] - B) < 0:
                ordered_elements = np.flip(ordered_elements)

        return ordered_elements
