# The Biostar Engine

Software for better science.

The **Biostar Engine** is a [Python][python] and [Django][django] is a collection of web applications that support science education.

Each web application may be deployed individually or in combination with all the others.

Specifically, the engine contains applications that support:

- Sharing and demonstrating data analysis scripts, see [Bioinformatics Recipes][recipes] (completed and operatinal)
- Community and social interactions via a Q&A (Question and Answer forum), see: [Biostars: Bioinformatics Explained][biostars] (in beta testing)
- Publishing documentation, see [Biostar Handbook][handbook] (in alpha stage)

[python]: https://www.python.org/
[django]: https://www.djangoproject.com/
[biostars]: https://www.biostars.org
[recipes]: https://www.bioinformatics.recipes
[handbook]: https://www.biostarhandbook.com
[conda]: https://conda.io/docs/

## Installation

The **Biostar Engine**  requires [Python 3.6][python] or above.

Our installation instructions rely on [conda][conda] though other alternatives are equally viable.  Users may use `virtualenv`, `pipenv`, `homebrew`, `apt-get` etc, or they may opt to not use environment management tools at all.

    # Create a virtual environment.
    conda create -y --name engine python=3.6
    source activate engine

    # Clone the source server code and the recipe code.
    git clone https://github.com/biostars/biostar-engine.git

    # Switch to the biostar-engine directory.
    cd biostar-engine

    # Install server dependencies.
    pip install -r conf/pip_requirements.txt

The installation is now complete. All server management commands run through `make`. To initialize and run the test site use:

      make reset serve

Visit `localhost:8000` to view the site. Use the command

    python manage.py

to see all the additional commands that may be executed.

To run bioinformatics oriented recipes further configuration of your environment may be necessary:

    conda config --add channels r
    conda config --add channels conda-forge
    conda config --add channels bioconda

    # Install the conda requirements.
    conda install --file conf/conda_requirements.txt

More information on [Deploying the Biostar Engine](docs/engine-deploy.md)
