FROM bioconductor/release_base2

ENV DEBIAN_FRONTEND=noninteractive
ENV R_LIBS_USER=$HOME/R
ENV R_LIBS_SITE=/opt/R
RUN apt-get -y update 
RUN apt-get -y install build-essential vim man
RUN apt-get -y install tk libglu1-mesa-dev mesa-common-dev mesa-utils libsdl2-2.0-0 libsdl2-dev glew-utils libglew-dev libsdl2-net-2.0-0 libsdl2-net-dev libsdl2-mixer-2.0-0 libsdl2-mixer-dev
RUN apt-get -y install libarchive-dev libc6-dbg gdb valgrind libfreetype6-dev squashfs-tools xclip
RUN mkdir /opt/R; R -e "install.packages(c('ggplot2','grid', 'gridExtra','pheatmap','R.utils', 'Rcpp', 'RcppArmadillo', 'rgl', 'roxygen2'), lib = '/opt/R/')" ; R -e "BiocManager::install(c('multtest','SingleCellExperiment', 'scmap','M3Drop', 'monocle', 'slingshot', 'DESeq2', 'DropletUtils'), lib = '/opt/R/')" ; R -e "install.packages(c('Seurat', 'gProfileR','Rtsne'), repos='http://cran.rstudio.com/', lib='/opt/R/')"
RUN wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /opt/miniconda.sh ; bash /opt/miniconda.sh -b -p /opt/miniconda ; /opt/miniconda/condabin/conda install -c bioconda -c conda-forge snakemake ; rm /opt/miniconda.sh
RUN apt-get -y install libarchive-dev python3 python3-pip python3-setuptools
RUN python3 -m pip install -U numpy; python3 -m pip install -U pandas; python3 -m pip install -U scipy; python3 -m pip install -U scanpy;  python3 -m pip install -U matplotlib ;  python3 -m pip install -U anndata
RUN python3 -m pip install -U scvi ;  python3 -m pip install -U scirpy;  python3 -m pip install -U ipywidgets
RUN python3 -m pip install -U vireoSNP ; python3 -m pip install -U scrublet; python3 -m pip install -U rbcde ; python3 -m pip install cellSNP
RUN apt-get -y install libcurl4-gnutls-dev
RUN R -e "install.packages(c('devtools', 'SoupX'), lib='/opt/R/'); library(devtools) ; with_libpaths(new = '/opt/R/', install_github('campbio/celda'))"

CMD ["helloworld"]

