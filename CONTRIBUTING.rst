.. highlight:: shell

============
Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

You can contribute in many ways:

Types of Contributions
----------------------

Report Bugs
~~~~~~~~~~~

Report bugs at https://github.com/CADWRDeltaModeling/schimpy/issues.

If you are reporting a bug, please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
~~~~~~~~

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
~~~~~~~~~~~~~~~~~~

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
~~~~~~~~~~~~~~~~~~~

schimpy could always use more documentation, whether as part of the
official schimpy docs, in docstrings, or even on the web in blog posts,
articles, and such.

Submit Feedback
~~~~~~~~~~~~~~~

The best way to send feedback is to file an issue at https://github.com/CADWRDeltaModeling/schimpy/issues.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)

Get Started!
------------

Ready to contribute? Here's how to set up `schimpy` for local development.

1. Fork the `schimpy` repo on GitHub.
2. Clone your fork locally::

    $ git clone git@github.com:your_name_here/schimpy.git

3. Install your local copy into a conda enviornment. Assuming you have conda installed, 
this is how you set up your fork for local development::

    $ # run commands in setup_devenv.bat on command line
    # you now have a dev_schimpy conda environment
    $ cd schimpy/
    $ conda activate dev_schimpy
    $ pip install -e . # Do this last as you don't want pip pulling in any dependencies

4. Create a branch for local development::

    $ git checkout -b name-of-your-bugfix-or-feature

   Now you can make your changes locally.

5. When you're done making changes, check that your changes pass the
   tests::

    $ cd tests
    $ pytest


6. Commit your changes and push your branch to GitHub::

    $ git add .
    $ git commit -m "Your detailed description of your changes."
    $ git push origin name-of-your-bugfix-or-feature

7. Submit a pull request through the GitHub website.

Pull Request Guidelines
-----------------------

Before you submit a pull request, check that it meets these guidelines:

1. The pull request should include tests.
2. If the pull request adds functionality, the docs should be updated. Put
   your new functionality into a function with a docstring, and add the
   feature to the list in README.rst.
3. The pull request should work for Python 3.5, 3.6, 3.7 and 3.8, and for conda build. 
    Do a local conda build 

    $ conda build conda.recipe
    
   and make sure that it builds a conda package

Tips
----

To run a subset of tests::

$ pytest tests.test_schimpy


Deploying
---------

A reminder for the maintainers on how to deploy.
Make sure all your changes are committed (including an entry in HISTORY.rst).
Then run::

$ git tag -a major.minor.patch
$ git push
$ git push --tags

Jenkins build will deploy to cadwr-dms channel within the hour
