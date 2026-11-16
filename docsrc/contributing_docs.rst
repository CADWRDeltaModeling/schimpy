Building the Documentation
==========================

Documentation sources are in ``docsrc/``. The normal build generates the API
reference from the current Python sources, then builds HTML in ``docs/html/``.
API generation is not a separate step after adding modules.

Environment
-----------

Open a Conda-initialized terminal at the repository root. Use
``environment.yml`` for the runtime/test dependencies and
``environment-doc.yml`` as a documentation overlay.

For a new environment, the convenience helper performs creation, the overlay
update, and editable installation. On Windows PowerShell:

.. code-block:: console

   .\build_doc_env.bat

On Linux/macOS, with Conda available in the terminal:

.. code-block:: console

   bash build_doc_env.sh

The helper does not delete or replace an existing ``schimpy_doc`` environment,
does not build the HTML, and does not leave your calling terminal activated.
After successful setup, run ``conda activate schimpy_doc``.
Use ``--help`` to see what the helper does.

Alternatively, perform the same steps explicitly:

.. code-block:: console

   conda env create -n schimpy_doc -f environment.yml
   conda env update -n schimpy_doc -f environment-doc.yml
   conda activate schimpy_doc
   python -m pip install --no-deps -e .

Do not use ``--prune`` with the overlay: it does not describe the full
environment. For an existing compatible documentation environment, use the
update step rather than deleting the environment.

Recent Conda releases supporting ``--format cep-24`` can combine the two files
in one creation command:

.. code-block:: console

   conda create -n schimpy_doc --format cep-24 --file environment.yml --file environment-doc.yml

Then activate the environment and install the editable checkout as above.
The overlay supplies Sphinx 8.2 or later within the 8.x series, required for
the built-in API generator. An older Sphinx environment cannot run this build.

Build and Preview
-----------------

With ``schimpy_doc`` active:

.. code-block:: console

   cd docsrc
   make html

On Windows PowerShell, use ``.\make.bat html``. On Linux/macOS,
``make html`` uses the Makefile; a ``make`` executable is required.
``doall.bat`` remains a compatibility alias for the same build.
For a clean rebuild, run ``make clean`` followed by ``make html``.

Open ``docs/html/index.html`` from the repository root in a browser. Check the
build messages for failed imports, inspect the API pages and the pages you
edited, and verify that the expected functions and CLI help are present.
A successful build alone does not establish complete API documentation.

Generated API sources live in ``docsrc/_api/`` and are ignored by git.
``schimpy.rst`` includes the generated package page, preserving the existing
API page location. The default Sphinx layout is used.
Do not edit generated API sources or commit generated HTML.

Add handwritten pages to the appropriate toctree and update Python docstrings
for API documentation. CI uses the same runtime YAML, documentation overlay,
editable checkout installation, and ``make html`` build, then publishes
``docs/html`` through GitHub Pages actions.

GitHub Pages
------------

In repository **Settings -> Pages -> Build and deployment**, select
**GitHub Actions** as the source. There is no branch or ``/docs`` source folder
to select: ``docs/html`` is the workflow's artifact path, not a folder served
from the repository.

The **Sphinx: Render docs** workflow runs on pushes to ``main``/``master``,
matching release tags, or manual dispatch from the **Actions** tab. The build
job uploads the HTML artifact; only a successful build can start the separate
deployment job. Check both jobs: a successful build does not establish that
deployment succeeded.

The ``github-pages`` environment must allow the publishing branch or tag.
If deployment is rejected, check **Settings -> Environments -> github-pages**
and any required reviewers or deployment-branch restrictions.

The canonical site is ``https://cadwrdeltamodeling.github.io/schimpy/``.
After the first successful deployment, check the home page, API reference and
CLI help. A stale ``gh-pages`` branch is not used by this workflow and need not
be deleted to make Actions deployment work.
Before retiring that branch, update any downstream intersphinx configuration
that reads its ``objects.inv`` to use
``https://cadwrdeltamodeling.github.io/schimpy/objects.inv`` instead.

Do not commit the generated site or add ``.nojekyll`` to support branch-based
publishing. Actions deployment does not use that mechanism. Pull-request checks,
branch protection and warnings-as-errors can be introduced separately.
