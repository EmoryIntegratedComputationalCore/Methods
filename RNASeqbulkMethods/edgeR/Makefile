# -*- coding: utf-8 -*- 
#
#Author: JR
#Maintainer: JR
#License:  2019, EICC, GPL v2 or later
# 
#Emory/EICC/Methods/RNASeqbulkMethods/edgeR/Makefile
#----------------------------------------------------
#

HERE := $(shell Rscript --vanilla -e 'write(here::here(), stdout())')

.PHONY:	all clean

all: output/edgeR.pdf

clean:
	
					-rm -r output/*
		
output/edgeR.pdf: \
					src/edgeR.Rmd
			Rscript --vanilla -e "rmarkdown::render('$<')"
			mv src/edgeR.pdf $@
# done.
