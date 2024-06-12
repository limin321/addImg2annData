FROM satijalab/seurat:4.3.0

RUN mkdir /home/anndata
COPY policy.xml addImage.sh h5ad2rds.R rds2h5adsc.R rds2annImg.py /usr/local/bin

RUN apt-get update && \
    apt-get install -y imagemagick bc python3.9 python3-pip && \
    rm /etc/ImageMagick-6/policy.xml && \
    cp /usr/local/bin/policy.xml /etc/ImageMagick-6/policy.xml && \
    pip install pandas matplotlib scanpy==1.9.6 && \
    pip install --upgrade numba numpy


RUN R -e "install.packages('remotes'); \
  remotes::install_version('data.table', '1.14.8')"
RUN R -e "install.packages('remotes'); \
  remotes::install_version('findpython', '1.0.8')"
RUN R -e "install.packages('remotes'); \
  remotes::install_version('R6', '2.5.1')"
RUN R -e "install.packages('remotes'); \
  remotes::install_version('jsonlite', '1.8.4')"
RUN R -e "install.packages('remotes'); \
  remotes::install_version('argparser', '0.7.2')"
RUN R -e "install.packages('remotes'); \
  remotes::install_version('dplyr', '1.1.2')"


# Set working dir to /home/anndata
WORKDIR /home/anndata

# Run your R script when the container starts
CMD ["bash", "addImage.sh"]



