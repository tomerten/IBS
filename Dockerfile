FROM continuumio/miniconda3

RUN conda install cmake make 
RUN conda install gcc_linux-64

RUN git clone https://github.com/tomerten/IBS.git

RUN apt-get update && apt-get -y install build-essential
RUN cd IBS && bash docker_build_all.sh

# setup for using on Binder
# see https://mybinder.readthedocs.io/en/latest/tutorials/dockerfile.htmla
#
# Note:
# Need to add user for running container -> --user 1000
ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

RUN cp -r /IBS ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}

