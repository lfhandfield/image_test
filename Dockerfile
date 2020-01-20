FROM bioconductor/release_base2

ENV DEBIAN_FRONTEND=noninteractive
ENV R_LIBS_USER=$HOME/R
ENV R_LIBS_SITE=/opt/R
RUN apt-get update 
RUN apt-get -y install build-essential vim man
RUN apt-get -y install libglu1-mesa-dev mesa-common-dev mesa-utils libsdl2-2.0-0 libsdl2-dev glew-utils libglew-dev libsdl2-net-2.0-0 libsdl2-net-dev 
#RUN apt-get -y install libc6-dbg gdb valgrind libfreetype6-dev
RUN R -e "install.packages(c('ggplot2','grid', 'gridExtra','pheatmap','R.utils', 'Rcpp', 'rgl', 'roxygen2'))" 
#RUN R -e "install.packages(c('Seurat', 'gProfileR','Rtsne'), repos='http://cran.rstudio.com/')"
#RUN R -e "BiocManager::install(c('SingleCellExperiment', 'scmap','M3Drop', 'monocle', 'DEseq2', 'DropletUtils'))" 
WORKDIR /opt/display
COPY src /opt/display/src
COPY Images /opt/display/Images
COPY Makefile /opt/display/Makefile
RUN make
ENV PATH="/opt/display/:${PATH}"
CMD ["helloworld"]
