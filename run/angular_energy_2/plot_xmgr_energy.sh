#!/bin/sh
sed 's/INFILE/files\/'''${1}'''0.dat/' < xmgr_batch.com > tmp1.com
sed 's/OUTFILE/print_files\/'''${1}'''0.pdf/' < tmp1.com > tmp2.com
gracebat -nosafe -batch tmp2.com
sed 's/INFILE/files\/'''${1}'''1.dat/' < xmgr_batch.com > tmp1.com
sed 's/OUTFILE/print_files\/'''${1}'''1.pdf/' < tmp1.com > tmp2.com
gracebat -nosafe -batch tmp2.com
sed 's/INFILE/files\/'''${1}'''2.dat/' < xmgr_batch.com > tmp1.com
sed 's/OUTFILE/print_files\/'''${1}'''2.pdf/' < tmp1.com > tmp2.com
gracebat -nosafe -batch tmp2.com
rm tmp1.com tmp2.com
