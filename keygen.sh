#!/bin/bash

yes | ssh-keygen -N "" -f docker-image/ssh-keys/id_rsa.mpi
chmod 700 docker-image/ssh-keys
chmod 600 docker-image/ssh-keys/*
