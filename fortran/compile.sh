#!/bin/bash


gfortran generator.f90 -o generator
gfortran generator_elipsoid.f90 -o generator_elipsoid
