#!/bin/bash

ssh -p 2222 -i docker-image/ssh-keys/id_rsa.mpi -o StrictHostKeyChecking=no mpi_user@localhost
