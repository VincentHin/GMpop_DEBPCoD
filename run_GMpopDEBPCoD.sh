#!/bin/bash

make allclean all

parallel < arglist.txt
