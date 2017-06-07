#!/bin/bash

ls -1 /hosts > /etc/openmpi/openmpi-default-hostfile
mpirun hostname
