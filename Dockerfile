FROM ubuntu:22.04

# https://serverfault.com/a/797318
ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y
RUN apt-get upgrade -y

RUN apt-get install -y git make cmake g++ libcurl4-openssl-dev


#
# Create a user
#

RUN useradd -m user

RUN apt-get install -y sudo
RUN echo "user ALL=(ALL) NOPASSWD:ALL" > /etc/sudoers

USER user
WORKDIR /home/user
ENV PATH="${PATH}:/home/user/.local/bin"
