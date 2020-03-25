# -*- coding: utf-8 -*- 
#
#Author: JR
#Maintainer: JR
#License:  2019, EICC, GPL v2 or later
# 
#Emory/EICC/Methods/RNASeqbulkMethods/DESeq2/Makefile
#----------------------------------------------------
#
HERE := $(shell Rscript --vanilla -e 'write(here::here(), stdout())')

.PHONY:	all clean

all: output/DESeq2.pdf

clean:
	
					-rm -r output/*
		
output/DESeq2.pdf: \
					src/DESeq2.Rmd
			Rscript --vanilla -e "rmarkdown::render('$<')"
			mv src/DESeq2.pdf $@
# done.
