#!/bin/bash

for file in ./*
do
    filename="${file%.*}.pin"
	cat $file | ../../../build/ARR >> $filename
done