Getting Started
===============
Mesh file
---------
To use the preprocessor, open the preprocessor.py file in the root directory and insert the name of the mesh file to be preprocessed and its dimension. For example, with the line below, IMPRESS would look for '20.h5m' file and assume that this mesh is tridimensional:

.. code:: python

   M = msh('20.h5m', dim = 3)

Coarsening Settings
-------------------
If you are interested in using a coarse mesh, it will be necessary to inform the partitioner scheme and coarsening ratio desired in the file msCoarse.yml. So far, IMPRESS handles only the Geometric Simple Partitioner Cube Based.


Variable Settings
-----------------

Run IMPRESS
-----------
Now, save preprocessor.py file and run it!
