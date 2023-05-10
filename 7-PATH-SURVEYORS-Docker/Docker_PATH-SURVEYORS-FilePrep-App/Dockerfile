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
RUN R -e 'install.packages(c("pbkrtest", "lme4", "car", "FactoMineR", "ggpubr", "shiny", "shinydashboard", "dplyr", "DT", "readr", "tools"))'
RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" > /usr/lib/R/etc/Rprofile.site
RUN addgroup --system app \
    && adduser --system --ingroup app app
COPY ./app/ /Users/alyss/OneDrive/Documents/R/Docker_PATH_SURVEYORS_FilePrep
RUN chown app:app -R /Users/alyss/OneDrive/Documents/R/Docker_PATH_SURVEYORS_FilePrep
RUN chown app:app -R /usr/local/lib/R/site-library
USER app
EXPOSE 3838
CMD ["R", "-e", "shiny::runApp('/Users/alyss/OneDrive/Documents/R/Docker_PATH_SURVEYORS_FilePrep', launch.browser = FALSE)"]
