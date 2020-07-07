FROM bioconductor/release_base2

ENV DEBIAN_FRONTEND=noninteractive
ENV R_LIBS_USER=$HOME/R
ENV R_LIBS_SITE=/opt/R
RUN apt-get update 
RUN apt-get -y install build-essential vim man
RUN apt-get -y install tk libglu1-mesa-dev mesa-common-dev mesa-utils libsdl2-2.0-0 libsdl2-dev glew-utils libglew-dev libsdl2-net-2.0-0 libsdl2-net-dev libsdl2-mixer-2.0-0 libsdl2-mixer-dev
RUN apt-get -y install libarchive-dev libc6-dbg gdb valgrind libfreetype6-dev squashfs-tools xclip
RUN mkdir /opt/R; R -e "install.packages(c('ggplot2','grid', 'gridExtra','pheatmap','R.utils', 'Rcpp', 'RcppArmadillo', 'rgl', 'roxygen2'), lib = '/opt/R/')" ; R -e "BiocManager::install(c('multtest','SingleCellExperiment', 'scmap','M3Drop', 'monocle', 'DESeq2', 'DropletUtils'), lib = '/opt/R/')" ; R -e "install.packages(c('Seurat', 'gProfileR','Rtsne'), repos='http://cran.rstudio.com/', lib='/opt/R/')"
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda.sh ; bash /opt/miniconda.sh -b -p /opt/miniconda ; /opt/miniconda/condabin/conda install -c bioconda -c conda-forge snakemake ; rm /opt/miniconda.sh
RUN apt-get -y install libarchive-dev
RUN apt-get update && apt-get install -y chromium-browser wget
# RUN wget https://github.com/singularityware/singularity/releases/download/2.5.2/singularity-2.5.2.tar.gz -O /opt/singularity.tar.gz ; tar xvf /opt/singularity.tar.gz ; rm /opt/singularity.tar.gz ; cd singularity-2.5.2; ./configure --prefix=/opt/singularity ; make ; sudo make install; cd .. ; rm -r -f singularity-2.5.2/
ENV PATH="/opt/display/:/opt/miniconda/bin/:${PATH}"
CMD ["helloworld"]
