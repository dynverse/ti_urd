FROM dynverse/dynwrap:bioc

RUN apt-get install -y libudunits2-dev

RUN R -e 'devtools::install_github("farrellja/URD")'

LABEL version 0.1.2

ADD . /code

ENTRYPOINT Rscript /code/run.R
