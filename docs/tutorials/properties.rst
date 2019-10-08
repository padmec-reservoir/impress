How to Access Mesh Entities' Properties
=======================================

To perform any kind of simulation, it will be necessary to access informations about the mesh entities such as coordinates from a node or from the center of an edge, for example. Probably, the user will also need informations about the neighbour entities (adjacencies) or even the internal or boundary elements as well. We call these are automatically generated informations of mesh entities' properties.

IMPRESS was developed to systematically provide access to these properties for any kind of mesh entity (nodes, edges, faces or volumes, in a 3D mesh). Every time IMPRESS is executed, several objects are created to represent the mesh entities and inherited by a main class (that represents the mesh itself) from which the user can access these objects.

SCHEME IMAGE


Nodes
-----
-> Coordinates
-> Boundary
-> Internal
-> Adjacencies
-> Flags
-> Flagged Elements
-> Global ID

Edges
-----
-> Center
-> Boundary
-> Internal
-> Adjacencies
-> Connectivities
-> Flags
-> Flagged Elements
-> Global ID

Faces
-----

-> Center
-> Boundary
-> Internal
-> Adjacencies
-> Connectivities
-> Flags
-> Flagged Elements
-> Global ID

Volumes
-------

-> Center
-> Boundary
-> Internal
-> Adjacencies
-> Connectivities
-> Flags
-> Flagged Elements
-> Global ID


Example:

.. code:: python
