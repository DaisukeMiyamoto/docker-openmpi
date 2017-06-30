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
    && mkdir -p ${HOME}/.ssh \
    && mkdir /var/run/sshd \
    && mkdir ${HOME}/mpi_examples \
    && chsh -s /bin/bash mpi_user

COPY ssh-keys/id_rsa.mpi ${HOME}/.ssh/id_rsa
COPY ssh-keys/id_rsa.mpi.pub ${HOME}/.ssh/authorized_keys
COPY ssh-keys/config ${HOME}/.ssh/config
COPY compose-start.sh ${HOME}
COPY make-hostfile.sh ${HOME}
COPY mpi_examples/ ${HOME}/mpi_examples

WORKDIR ${HOME}
RUN  chown -R mpi_user:mpi_user ${HOME} \
    && chmod 700 .ssh \
    && chmod 600 .ssh/* \
    && chmod +x *.sh \
    && touch /etc/openmpi/openmpi-default-hostfile \
    && chmod 666 /etc/openmpi/openmpi-default-hostfile
RUN cd mpi_examples \
    && make

EXPOSE 22
CMD ["/usr/sbin/sshd", "-D"]
