################################
# Makefile
#
# author: He Zhang
# edited by: 11/2019
################################

CC=g++
DEPS=src/LinearSampling.h src/backtrace.cpp src/Utils/energy_parameter.h src/Utils/feature_weight.h src/Utils/intl11.h src/Utils/intl21.h src/Utils/intl22.h src/Utils/utility_v.h src/Utils/utility.h src/Utils/logspace.h

CFLAGS=-std=c++11 -O3
.PHONY : clean linearsampling
objects=bin/linearsampling_v

linearsampling: src/LinearSampling.cpp $(DEPS) 
		chmod +x linearsampling
		mkdir -p bin
		$(CC) src/LinearSampling.cpp $(CFLAGS) -o bin/linearsampling_v 

clean:
	-rm $(objects)