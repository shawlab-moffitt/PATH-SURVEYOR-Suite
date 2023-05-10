FROM rocker/r-base:latest
LABEL maintainer="Tim Shaw <timothy.shaw@moffitt.org>"
RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*
RUN install.r shiny
RUN R -e 'install.packages(c("pbkrtest", "lme4", "car", "FactoMineR", "ggpubr", "shiny","shinydashboard","shinythemes","shinyjqui","pheatmap","RColorBrewer","ggdendro","factoextra","dplyr","DT","viridis","readr","shinycssloaders","stringr","tools","plotly","reshape2", "BiocManager"))'
RUN R -e 'BiocManager::install("clusterProfiler", ask = F)'
RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" > /usr/lib/R/etc/Rprofile.site
RUN addgroup --system app \
    && adduser --system --ingroup app app
#WORKDIR "/Users/4472414/DockerShiny/Jaccard_Pathway_Connectivity"
COPY ./app/ /Users/4472414/DockerShiny/Jaccard_Pathway_Connectivity
RUN chown app:app -R /Users/4472414/DockerShiny/Jaccard_Pathway_Connectivity
RUN chown app:app -R /usr/local/lib/R/site-library
USER app
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/Users/4472414/DockerShiny/Jaccard_Pathway_Connectivity', launch.browser = FALSE)"]
