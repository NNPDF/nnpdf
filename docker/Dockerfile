# Download base image ubuntu latest
FROM ubuntu:latest

# LABEL about the custom image
LABEL maintainer="NNPDF"
LABEL description="Image with NNPDF tagged code."

ENV CONDA_URL https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
ENV CONDA_FILE Miniconda3-latest-Linux-x86_64.sh

# Disable Prompt During Packages Installation
ARG DEBIAN_FRONTEND=noninteractive

# # Update Ubuntu Software repository
RUN apt update

# Install wget
RUN apt install -y wget

# Istall conda
RUN wget "$CONDA_URL"      && \
    chmod +x "$CONDA_FILE" && \
    ./"$CONDA_FILE" -b     && \
    rm "$CONDA_FILE"

ENV PATH /root/miniconda3/bin:$PATH

# Setup conda channels
RUN conda config --append channels conda-forge                                    && \
    conda config --prepend channels https://packages.nnpdf.science/public/        && \
    conda config --set show_channel_urls true                                     && \
    conda init bash

# Copy conda package to image
COPY conda-bld /root/conda-bld

# install the nnpdf framework
RUN conda create -n nnpdfenv -c file:///root/conda-bld/linux-64 nnpdf && \
    echo "conda activate nnpdfenv" >> /root/.bashrc

# delete package
RUN rm -r /root/conda-bld

# run bash as entrypoint
ENTRYPOINT bash
