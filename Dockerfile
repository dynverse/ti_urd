FROM dynverse/dynwrap_latest:v0.1.0

ARG GITHUB_PAT

RUN apt-get update && apt-get install -y libudunits2-dev

RUN R -e 'devtools::install_github("farrellja/URD")'

COPY definition.yml run.R example.sh /code/

ENTRYPOINT ["/code/run.R"]
