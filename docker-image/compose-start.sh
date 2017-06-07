#!/bin/bash

touch /hosts/$(hostname)
/usr/sbin/sshd -D
