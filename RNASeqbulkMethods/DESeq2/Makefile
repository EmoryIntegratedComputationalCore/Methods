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

all: output/DEGAnalysiswithDESeq2.pdf

clean:
	
					-rm -r output/*
		
output/DEGAnalysiswithDESeq2.pdf: \
					src/DEGAnalysiswithDESeq2.Rmd
			-rm src/DEGAnalysiswithDESeq2.pdf
			Rscript --vanilla -e "rmarkdown::render('$<')"
			mv src/DEGAnalysiswithDESeq2.pdf $@
# done.
