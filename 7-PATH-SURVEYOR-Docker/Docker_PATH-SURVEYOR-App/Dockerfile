FROM rocker/r-base:latest
LABEL maintainer="Tim Shaw <timothy.shaw@moffitt.org>"
RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    libcurl4-gnutls-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    libssh2-1-dev \
    libv8-dev \
    libnode-dev \
    cmake \
    && rm -rf /var/lib/apt/lists/*
RUN install.r shiny
RUN R -e 'install.packages(c("pbkrtest", "lme4", "car", "FactoMineR", "ggpubr", "shiny","shinydashboard","shinythemes","shinyjqui","pheatmap","RColorBrewer","ggdendro","factoextra","dplyr","DT","viridis","readr","shinycssloaders","stringr","tools","plotly","reshape2", "BiocManager","ggplot2","survival","readr","shinycssloaders","survminer","gridExtra","ggrepel"))'
RUN R -e 'BiocManager::install("clusterProfiler", ask = F)'
RUN apt-get update && apt-get install -y --no-install-recommends \
    libxml2 \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

RUN R -e 'install.packages(c("remotes", "juicyjuice", "xml2"))'
RUN R -e 'remotes::install_github("rstudio/gt")'
RUN R -e 'install.packages(c("gtsummary"))'
RUN R -e 'install.packages(c("gridtext", "ggtext", "survminer"))'
RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" > /usr/lib/R/etc/Rprofile.site
RUN addgroup --system app \
    && adduser --system --ingroup app app
COPY ./app/ /Users/4472414/DockerShiny/PATH-SURVEYOR-UserInput
RUN chown app:app -R /Users/4472414/DockerShiny/PATH-SURVEYOR-UserInput
RUN chown app:app -R /usr/local/lib/R/site-library
RUN R -e 'BiocManager::install("GSVA", ask = F)'
USER app
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/Users/4472414/DockerShiny/PATH-SURVEYOR-UserInput', launch.browser = FALSE)"]
#CMD ["R", "-e", "BiocManager::install('GSVA', ask = F)"]
#CMD ["R", "-e", "install.packages('https://cran.r-project.org/src/contrib/survminer_0.4.9.tar.gz', repos=NULL)"]
#CMD ["R", "-e", "install.packages('https://cran.r-project.org/src/contrib/ggtext_0.1.2.tar.gz', repos=NULL)"]
#CMD ["R", "-e", "install.packages('https://cran.r-project.org/src/contrib/gridtext_0.1.5.tar.gz', repos=NULL)"]
#CMD ["R", "-e", "remotes::install_github('rstudio/gt')"]
#CMD ["R", "-e", "install.packages('https://cran.r-project.org/src/contrib/xml2_1.3.4.tar.gz', repos=NULL)"]
#CMD ["R", "-e", "install.packages('https://cran.r-project.org/src/contrib/juicyjuice_0.1.0.tar.gz', repos=NULL)"]
#CMD ["R", "-e", "install.packages('', repos=NULL)"]
#CMD ["R", "-e", "install.packages('https://cran.r-project.org/src/contrib/V8_4.3.0.tar.gz', repos=NULL)"]
