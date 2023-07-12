FROM ubuntu:jammy

# Change the default shell from sh to bash
#SHELL ["/bin/bash", "-c"]

# Prevent user interactions required by r-base
ENV DEBIAN_FRONTEND=noninteractive 
# Update package list and install required packages
USER root
RUN apt-get update && \
    apt-get install -y \
    zip \
    libnetcdf-dev \
    openssh-client \
    r-base \
    gdal-bin \
    libgdal-dev \
    build-essential

# Install R packages
RUN R -e "install.packages('chron')"
RUN R -e "install.packages('fields')"
RUN R -e "install.packages('gplots')"
RUN R -e "install.packages('ncdf4')"
RUN R -e "install.packages('parallel')"
#RUN R -e "install.packages('rgdal')"
RUN R -e "install.packages('raster')"
RUN R -e "install.packages('sp')"
RUN R -e "install.packages('zoo')"
RUN R -e "install.packages('apcluster')"
RUN R -e "install.packages('compiler')"
RUN R -e "install.packages('RColorBrewer')"
RUN R -e "install.packages('colorspace')"
RUN R -e "install.packages('maps')"

# Create a user
ARG USER=cardamom
RUN adduser --disabled-password --gecos '' $USER
RUN adduser $USER sudo; echo '%sudo ALL=(ALL) NOPASSWD:ALL' >> /etc/sudoers
RUN chown -R $USER:$USER /home/$USER
USER $USER
ENV HOME /home/$USER
ENV USER $USER
WORKDIR $HOME

# Switch to the new user
USER $USER

# Install CARDAMOM
RUN mkdir CARDAMOM
ADD --chown=$USER:$USER . CARDAMOM/

