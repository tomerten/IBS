FROM continuumio/miniconda3

RUN conda install cmake make 
RUN conda install gcc_linux-64

RUN git clone https://github.com/tomerten/IBS.git

RUN apt-get update && apt-get -y install build-essential
RUN cd IBS && bash docker_build_all.sh

