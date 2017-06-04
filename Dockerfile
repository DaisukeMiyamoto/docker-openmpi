FROM ubuntu:16.04

MAINTAINER Daisuke Miyamoto <miyamoto@brain.imi.i.u-tokyo.ac.jp>
ENV HOME /home/mpi_user

RUN apt-get update \
    && apt-get install -y \
        locales \
        wget \
        gcc \
        g++ \
        build-essential \
        openmpi-bin \
        openmpi-common \
        libopenmpi-dev \
        ssh \
    && localedef -i en_US -c -f UTF-8 -A /usr/share/locale/locale.alias en_US.UTF-8 \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get autoclean

RUN useradd mpi_user \
    && mkdir -p ${HOME} \
    && mkdir -p {HOME}/.ssh \
    && mkdir /var/run/sshd

COPY ssh-keys/id_rsa.mpi ${HOME}/.ssh/id_rsa
COPY ssh-keys/id_rsa.mpi.pub ${HOME}/.ssh/authorized_keys

WORKDIR ${HOME}
RUN  chown -R mpi_user:mpi_user ${HOME} \
    && chmod 700 .ssh \
    && chmod 600 .ssh/*

EXPOSE 22
#USER mpi_user
CMD ["/usr/sbin/sshd", "-D"]
#RUN chmod 600 .ssh/*

