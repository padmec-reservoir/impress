Via Repo
--------

The only way to use IMPRESS is to explicitly download it from our `repo <https://github.com/padmec-reservoir/impress>`_ or simply `fork <https://help.github.com/en/articles/fork-a-repo>`_ it.

Dependencies
------------

The simplest way to run IMPRESS is to use `Docker <https://www.docker.com/>`_ as a container software with a image developed by our team that provides all the dependecies of IMPRESS. The image file is stored in a directory called ``docker-image``. In that case, the procedure to use IMPRESS is:

* `Install <https://www.docker.com/get-started>`_ Docker correctly;
* Build the image stored in the directory called ``docker-image``;

  1. Open your terminal and navigate to ``docker-image`` directory
  2. Build the ``Dockerfile``

.. code:: python

   docker build -t image_name .

* Run Docker using the image you've just built through your terminal with the command:

.. code:: python

   sudo docker run -t -it -v path/to/repo/directory:/root image_name bash -c "cd /root; ipython"

If all steps were completed properly, you're running `iPython <https://ipython.org/>`_ inside your recently created container and all files from IMPRESS repo must be available to be accessed.
