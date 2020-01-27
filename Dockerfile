FROM bioconductor/release_base2

ENV DEBIAN_FRONTEND=noninteractive
ENV R_LIBS_USER=$HOME/R
ENV R_LIBS_SITE=/opt/R
RUN apt-get update 
RUN apt-get -y install build-essential vim man
RUN apt-get -y install tk libglu1-mesa-dev mesa-common-dev mesa-utils libsdl2-2.0-0 libsdl2-dev glew-utils libglew-dev libsdl2-net-2.0-0 libsdl2-net-dev libsdl2-mixer-2.0-0 libsdl2-mixer-dev
RUN apt-get -y install libc6-dbg gdb valgrind libfreetype6-dev
RUN mkdir /opt/R
RUN R -e "install.packages(c('ggplot2','grid', 'gridExtra','pheatmap','R.utils', 'Rcpp', 'RcppArmadillo', 'rgl', 'roxygen2'), lib = '/opt/R/')" 
RUN R -e "BiocManager::install(c('multtest','SingleCellExperiment', 'scmap','M3Drop', 'monocle', 'DESeq2', 'DropletUtils'))" 
RUN R -e "install.packages(c('Seurat', 'gProfileR','Rtsne'), repos='http://cran.rstudio.com/', lib='/opt/R/')"
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda.sh
RUN bash /opt/miniconda.sh -b -p /opt/miniconda 
RUN rm /opt/miniconda.sh
WORKDIR /opt/display
COPY src /opt/display/src
COPY Images /opt/display/Images
COPY Makefile /opt/display/Makefile
RUN make
ENV PATH="/opt/display/:${PATH}"
CMD ["helloworld"]
