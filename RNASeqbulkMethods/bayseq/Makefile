# -*- coding: utf-8 -*- 
#
#Author: JR
#Maintainer: JR
#License:  2019, EICC, GPL v2 or later
# 
#Emory/EICC/Methods/RNASeqbulkMethods/bayseq/Makefile
#----------------------------------------------------
#

HERE := $(shell Rscript --vanilla -e 'write(here::here(), stdout())')

.PHONY:	all clean

all: output/bayseq.pdf

clean:
	
					-rm -r output/*
		
output/bayseq.pdf: \
					src/bayseq.Rmd
			-rm src/bayseq.pdf
			Rscript --vanilla -e "rmarkdown::render('$<')"
			mv src/bayseq.pdf $@
# done.
