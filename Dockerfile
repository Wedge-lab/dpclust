FROM ubuntu:16.04

USER root

RUN apt-get update && apt-get -y install r-base

RUN R -q -e 'source("http://bioconductor.org/biocLite.R"); biocLite(c("optparse","KernSmooth","ks","lattice","ggplot2","gridExtra"))'

RUN mkdir -p /opt/dpclust
COPY . /opt/dpclust/
RUN R -q -e 'install.packages("/opt/dpclust", repos=NULL, type="source")'

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
