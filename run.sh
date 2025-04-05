#!/bin/bash
g++ -fPIC -shared -o libWignerUtils.so CWignerUtils.cpp CWignerSource.cpp `root-config --cflags --libs`