How to Access Mesh Entities Properties
======================================

To perform any kind of simulation, it will be necessary to access informations about the mesh entities such as coordinates from a node or from the center of an edge, for example. Probably, the user will also need informations about the adjacents elements or even the internal or boundary elements as well. We call these are automatically generated informations of **properties**.

IMPRESS was developed to systematically provide access to these properties for any kind of mesh entity (nodes, edges, faces or volumes, in a 3D mesh). Every time IMPRESS is executed, several objects are created to represent the mesh entities and inherited by a main class (that represents the mesh itself), FineScaleMesh. A ilustrative scheme of the inheritance scheme follows below:

SCHEME IMAGE

The IMPRESS' execution script instatiates this class creating an object called **M** through which it's possible to access all mesh entities. Furthermore, the user can obtain the mesh entity properties through its objects (which are literally called `nodes`, `edges`, `faces` and `volumes`).

Properties
----------
The **properties** that IMPRESS generates and the objects that represents these properties are described below:

* **Coordinates**: returns the coordinates of an array of elements;

* **Boundary Elements**: returns the global id from all elements located in the mesh boundaries;

* **Internal Elements**: returns the global id from all elements that do not belong to the mesh boundaries;

* **Adjacencies**: returns global id from all elements of the same dimension that are immediately connected to the element;

* **Connectivities**:

* **Flags**:

* **Flagged Elements**:

* **Global ID**:

Consulting Properties
---------------------

.. code:: python

  In [4]: M.nodes.coords[:]
  Out[4]:
  array([[ 0.,  0.,  0.],
        [ 1.,  0.,  0.],
        [ 2.,  0.,  0.],
        ...,
        [18., 20., 20.],
        [19., 20., 20.],
        [20., 20., 20.]])

Coarse Scale Mesh Entities Properties
-------------------------------------

Coming soon!
