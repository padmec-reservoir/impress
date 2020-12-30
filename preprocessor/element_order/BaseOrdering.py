class BaseOrdering(object):
  """
  Base class that defines the basic methods and
  attributes for element sorting. It should NOT
  be instantiated.
  """
  
  def __init__(self, mesh, target_dim, interface_dim):
    self.mesh = mesh
    self.target_dim = target_dim
    self.interface_dim = interface_dim
  
  def sort_elements(self):
    raise NotImplementedError
