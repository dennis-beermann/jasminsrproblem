FROM ghcr.io/rocker-org/devcontainer/r-ver:4.2

# Install R's X11 runtime dependencies
RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libx11-6 \
    libxss1 \
    libxt6 \
    libxext6 \
    libsm6 \
    libice6 \
    xdg-utils \
  && rm -rf /var/lib/apt/lists/*

# Install specific devtools version from source since we need it to install
# all other packages
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/usethis/usethis_2.1.6.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/devtools/devtools_2.4.4.tar.gz', repos=NULL, type='source')"
RUN R -e "install.packages('https://cran.r-project.org/src/contrib/Archive/cpp11/cpp11_0.4.4.tar.gz', repos=NULL, type='source')"

# # Install other R packages using devtools
RUN R -e "install.packages('rmarkdown')"
RUN R -e "install.packages('gghighlight')"
RUN R -e "install.packages('GGally')"
RUN R -e "remotes::install_github('vqv/ggbiplot')"
RUN R -e "devtools::install_github('Vitek-Lab/MSstatsConvert')"
RUN R -e "devtools::install_github('Vitek-Lab/MSstats')"
RUN R -e "devtools::install_github('Vitek-Lab/MSstatsTMT')"
RUN R -e "devtools::install_github('Vitek-Lab/MSstatsPTM')"
RUN R -e "devtools::install_github('Vitek-Lab/MSstatsLiP', build_vignettes = TRUE)"