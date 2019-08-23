## Biostar Central

### Software for better science.

**Biostar Central** is a [Python][python] and [Django][django] based collection of web applications that support scientific practice and education. Each web application may be deployed individually or in combination with all the others.

The goal of the project is to produce software with straigtforward installation, and minimal dependencies that works on any computing platform that supports Python.

Applications that are currently feature complete:

- `recipes` a web app that runs scripts via a web interface, see: [Bioinformatics Recipes][recipes]
- `forum` a web app that runs a Q&A forum inspired by StackOverflow, see: [Biostars Q&A][biostars]

Note (July 25, 2019): The public site runs the code from the `biostar2016` branch. The version in the master branch is the development version that we plan to release later.

In addition the `biostar-central` code repository includes applications that provide generic utility.

- `accounts` an app that manages user accounts
- `emailer` an app that can send emails

Future plans include adding applications such as:

- `publisher` a web app that allows publishing static webpages with various group level access rights, see [Biostar Handbook][handbook]
- `courses` a web app that allows for publishing training materials, courses, quizzes over the web, see the courses in the [Biostar Handbook][handbook]
- `archive` a web application that allows publishing scientific works.

In each case the driving philopsphy is that of decentralization. We write code that allows others to recreate the same services that we run.

[python]: https://www.python.org/
[django]: https://www.djangoproject.com/
[biostars]: https://www.biostars.org
[recipes]: https://www.bioinformatics.recipes
[handbook]: https://www.biostarhandbook.com
[conda]: https://conda.io/docs/


## Installation

The code in **Biostar Central**  requires [Python 3.6][python] or above.

Our installation instructions rely on [conda][conda] though other alternatives for managing python environments are equally viable.


    # Create a virtual environment.
    conda create -y --name engine python=3.6
    
    # Activate the python environment.
    conda activate engine

    # Clone the source server code and the recipe code.
    git clone https://github.com/ialbert/biostar-central.git

    # Switch to the biostar-engine directory.
    cd biostar-central

    # Install server dependencies.
    pip install -r conf/pip_requirements.txt

The installation is now complete. To test the code run:

    make forum test
    make recipes test

## Run

All server management commands run through `make`. We provide different initializations depending on each app. For example, to run a demonstration version of the `forum` app execute:

    make forum_demo

To run the demonstration version of the `recipes` app execute:

    make recipes_demo

In each case visit <http://localhost:8000> to view the site.

The Makefile has several tasks that demonstrate the commands that may be chained together:

    make forum init serve
    make recipes init serve


## Documentation

Additional documentation in [docs/index.md](docs/README.md)


