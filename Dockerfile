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
RUN R -e "install.packages('ggplot2',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages('gridExtra',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages('SingleCellExperiment',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages('pheatmap',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages('Seurat',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('DEseq2',dependencies=TRUE, repos='http://cran.rstudio.com/')"
RUN R -e "install.packages('scmap',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages('M3Drop',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages('rgl',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages('gProfileR',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages('monocle',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages('grid',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages('Rtsne',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages('R.utils',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
RUN R -e "install.packages('DropletUtils',dependencies=TRUE, repos='http://cran.rstudio.com/')" 
ENV PATH="/opt/display/:${PATH}"
CMD ["helloworld"]
