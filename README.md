## Biostar Central

### Software for better science.

**Biostar Central** is a [Python][python] and [Django][django] based collection of web applications that support science education. Each web application may be deployed individually or in combination with all the others.

Applications that are currently feature complete:

- `recipes` is an app that runs scripts via a web interface. See: [Bioinformatics Recipes][recipes]
- `forum` is an app that runs a Q&A forum inspired by StackOverflow. See: [Biostars Q&A][biostars]

In addition `biostar-central` includes applications that provide generic utility.

- `accounts` is an app that manages user accounts
- `emailer` is an app that can send emails

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
    make recipe test

## Run

All server management commands run through `make`. We provide different initializations depending on each app. For example, to run a demonstration version of the `forum` app execute:

    make forum_demo

To run the demonstration version of the `recipes` app execute:

    make recipes_demo

In each case visit `http://localhost:8000` to view the site.

The Makefile has several tasks that demonstrate the main be chained together:

    make forum init serve
    make recipes init serve


## Documentation

---

> **Note**
>
> *The new version of the [Biostar Q&A][biostars] is in beta testing!*
> *The public site currently runs the code in the `biostar2016` branch. We will switch to the new version once testing completes.*

---

Additional documentation in [docs/index.md](docs/index.md)


