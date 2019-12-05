FROM bioconductor/release_base2

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update 
RUN apt-get -y install build-essential
RUN apt-get -y install libglu1-mesa-dev mesa-common-dev mesa-utils libsdl2-2.0-0 libsdl2-dev glew-utils libglew-dev libsdl2-net-2.0-0 libsdl2-net-dev 
RUN apt-get -y install libc6-dbg gdb valgrind
WORKDIR /opt/display
COPY src /opt/display/src
COPY Images /opt/display/Images
COPY Makefile /opt/display/Makefile
RUN ls
RUN make
RUN R -e "install.packages(c('ggplot2','grid', 'gridExtra','pheatmap'),dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages(c('Seurat', 'gProfileR','Rtsne'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "BiocManager::install('SingleCellExperiment', 'scmap','M3Drop', 'monocle', 'DEseq2', 'DropletUtils')" 
RUN R -e "install.packages('rgl',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages('R.utils',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
ENV PATH="/opt/display/:${PATH}"
CMD ["helloworld"]
