ARG itag=latest
ARG r_version=3.5.2
ARG base_name=r-${r_version}-mkl
FROM uwgac/$base_name:$itag

RUN git clone https://github.com/broadinstitute/TOPMed_AFib_pipeline.git

RUN R CMD BATCH /TOPMed_AFib_pipeline/install.R
