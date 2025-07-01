Installation
============

``autoDCR`` runs on Python3. As of version 0.3.0, it's made using `poetry <https://python-poetry.org/docs/>`_ and `typer <https://typer.tiangolo.com/>`_. It can be installed from PyPI:

``pip install autoDCR``

(Note that the `last non-package version (0.2.7) is available as a release from the Github repo <https://github.com/JamieHeather/autoDCR/releases/tag/v0.2.7>`_, which is the version used in the `2022 stitchr paper <https://doi.org/10.1093/nar/gkac190>`_.)

In order to download the necessary TCR germline reference data, `IMGTgeneDL <https://github.com/JamieHeather/IMGTgeneDL>`_ is also required. If it's not automatically installed alongside ``stitchr``, it can be installed with:

``pip install IMGTgeneDL``

General usage
-------------

* AutoDCR is run in the command line using a number of different subcommands, e.g. ``autoDCR subcommand -x -y -z parameters``.
* The ``--help / -h`` flag can be used to list either:
    * The different subcommands available (``autoDCR -h``)...
    * ... Or the parameters associated with a subcommand (``autoDCR subcommand -h``).
* This help page gives details as to the types and default values for the different editable parameters.

