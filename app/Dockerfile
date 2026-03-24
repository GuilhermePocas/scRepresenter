
FROM rocker/shiny:4.5 AS builder

RUN apt-get update && apt-get install -y \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libhdf5-dev \
    libpng-dev \
    libjpeg-dev \
    libcairo2-dev \
    libxt-dev \
    python3 python3-venv python3-pip \
    && rm -rf /var/lib/apt/lists/*

RUN python3 -m venv /opt/pyzell \
    && /opt/pyzell/bin/pip install --upgrade pip \
    && /opt/pyzell/bin/pip install numpy pandas anndata h5py scipy

RUN R -e "install.packages(c('shiny','plotly','DT','ggplot2','dplyr','tidyr','networkD3','RColorBrewer','htmlwidgets'), repos='https://cloud.r-project.org/')"
RUN R -e "if (!requireNamespace('BiocManager')) install.packages('BiocManager'); BiocManager::install(c('SingleCellExperiment','anndataR','Seqinfo','SummarizedExperiment','GenomicRanges','S4Vectors', 'rhdf5'), ask=FALSE)"

COPY app.R /srv/shiny-server/

RUN echo "options(shiny.maxRequestSize = 10 * 1024^3)" >> /usr/local/lib/R/etc/Rprofile.site

EXPOSE 3838

ENV RETICULATE_PYTHON=/opt/pyzell/bin/python3

CMD ["/usr/bin/shiny-server"]