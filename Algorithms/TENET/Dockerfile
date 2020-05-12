FROM ubuntu:18.04

LABEL maintainer = "Zijun Ma <mazijun2017@outlook.com>"

USER root

RUN mkdir runTENET/

COPY runTENET/ /runTENET/

WORKDIR runTENET/

RUN mkdir data/

RUN apt-get update && apt-get install -y python3 python3-pip time libblacs-mpi-dev openssh-server default-jre

RUN pip3 install JPype1 numpy


