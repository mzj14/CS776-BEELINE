FROM ubuntu:18.04

LABEL maintainer = "Zijun Ma <mazijun2017@outlook.com>"

USER root

WORKDIR /

COPY runCorto.R /

RUN mkdir data/

RUN apt-get update && apt-get install -y tzdata software-properties-common

RUN apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9

RUN add-apt-repository "deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/"

RUN apt-get update && apt-get install -y r-base time

RUN R -e "install.packages('corto')"

RUN R -e "install.packages('optparse')"
