FROM continuumio/miniconda3:4.11.0

COPY . /biostar-central/
WORKDIR /biostar-central

RUN conda env create -f environment.yml

RUN apt-get update && \
    apt-get -y install build-essential

SHELL ["conda", "run", "-n", "engine", "/bin/bash", "-c"]

RUN pip install -r /biostar-central/conf/requirements.txt

RUN conda config --add channels r
RUN conda config --add channels conda-forge
RUN conda config --add channels bioconda
RUN conda install --file /biostar-central/conf/conda-packages.txt

RUN python manage.py migrate --settings biostar.forum.settings
RUN python manage.py collectstatic --noinput -v 0 --settings biostar.forum.settings

RUN make forum init 
RUN make forum test

ENTRYPOINT ["make", "forum", "serve"]
