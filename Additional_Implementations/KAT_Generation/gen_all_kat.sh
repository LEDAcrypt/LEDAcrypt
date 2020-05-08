#!/bin/bash

for i in build/bin/KEM-CPA/*
do
	./$i
done
for i in build/bin/KEM/*
do
	./$i
done
for i in build/bin/PKC/*
do
	./$i
done
