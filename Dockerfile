FROM ubuntu:24.04
ARG DEBIAN_FRONTEND=noninteractive

# 1. Install R and System Dependencies
# r-base-dev is CRITICAL for compiling the new Rcpp C++ code
RUN apt-get update && apt-get install -y \
    r-base \
    r-base-dev \
    libcurl4-gnutls-dev \
    libxml2-dev \
    libssl-dev \
    libcairo2-dev \
    libxt-dev \
    libgomp1 \
    make \
    curl \
    git \
    bedtools \
    pigz \
    tabix \
    libdeflate-tools \
    && rm -rf /var/lib/apt/lists/*

# 2. OPTIMIZATION: Configure Posit Binary Repository for Ubuntu Noble
# This allows installing R packages as pre-compiled binaries (10x faster)
RUN mkdir -p /usr/lib/R/etc && \
    echo 'options(repos = c(CRAN = "https://packagemanager.posit.co/cran/__linux__/noble/latest"))' >> /usr/lib/R/etc/Rprofile.site && \
    echo 'options(HTTPUserAgent = sprintf("R/%s R (%s)", getRversion(), paste(getRversion(), R.version$platform, R.version$arch, R.version$os)))' >> /usr/lib/R/etc/Rprofile.site

# 3. Install pak (Fast dependency manager)
RUN Rscript -e "install.packages('pak', repos = 'https://r-lib.github.io/p/pak/stable')"

WORKDIR /opt/dpclust

# 4. OPTIMIZATION: Cache dependency installation layer
# Copy DESCRIPTION first so that changes to code don't invalidate the dependency cache.
COPY DESCRIPTION .
# Install dependencies defined in DESCRIPTION (Rcpp, ggplot2, etc.)
RUN Rscript -e "pak::local_install_dev_deps()"

# 5. Copy the rest of the code and install the package
COPY . .
# Remove any local artifacts to ensure clean build
RUN rm -rf src/*.o src/*.so
# Install DPClust (compiles the C++ code here)
RUN Rscript -e "pak::local_install()"

USER ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]