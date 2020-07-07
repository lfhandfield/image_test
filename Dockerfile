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
# install manually all the missing libraries
#RUN wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb
#RUN dpkg -i google-chrome-stable_current_amd64.deb --fix-missing; apt-get -fy install

#RUN apt-get install -y --fix-missing gconf-service libasound2 libatk1.0-0 libcairo2 libcups2 libfontconfig1 libgdk-pixbuf2.0-0 libgtk-3-0 libnspr4 libpango-1.0-0 libxss1 fonts-liberation libappindicator1 libnss3 lsb-release xdg-utils
# install chrome
#RUN wget https://dl.google.com/linux/direct/google-chrome-stable_current_amd64.deb
#RUN dpkg -i google-chrome-stable_current_amd64.deb; apt-get -fy install
# RUN wget https://github.com/singularityware/singularity/releases/download/2.5.2/singularity-2.5.2.tar.gz -O /opt/singularity.tar.gz ; tar xvf /opt/singularity.tar.gz ; rm /opt/singularity.tar.gz ; cd singularity-2.5.2; ./configure --prefix=/opt/singularity ; make ; sudo make install; cd .. ; rm -r -f singularity-2.5.2/
#ENV PATH="/opt/display/:/opt/miniconda/bin/:${PATH}"
#RUN add-apt-repository universe ; apt update ; apt upgrade
RUN apt-get install -y software-properties-common
RUN apt-add-repository 'deb http://ftp.uk.debian.org/debian stretch main contrib non-free'
RUN cat /etc/apt/sources.list
RUN apt-add-repository -r 'deb http://ftp.uk.debian.org/debian stretch main contrib non-free'
RUN apt-get update && apt upgrade -y
RUN apt install chromium-browser
CMD ["helloworld"]
