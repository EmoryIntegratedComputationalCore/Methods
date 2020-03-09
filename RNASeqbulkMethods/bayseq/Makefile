# -*- coding: utf-8 -*- 
#
#Author: JR
#Maintainer: JR
#License:  2019, EICC, GPL v2 or later
# 
#Emory/EICC/Methods/RNASeqbulkMethods/baySeq/Makefile
#----------------------------------------------------
#

HERE := $(shell Rscript --vanilla -e 'write(here::here(), stdout())')

.PHONY:	all clean

all:	output/DEGAnalysiswithbayseq_underconstruction.pdf

clean:
					-rm -r output/*
	
	
output/DEGAnalysiswithbayseq_underconstruction.pdf: \
					src/DEGAnalysiswithbayseq_underconstruction.Rmd
			-rm src/DEGAnalysiswithbayseq_underconstruction.pdf
			Rscript --vanilla -e "rmarkdown::render('$<')"
			mv src/DEGAnalysiswithbayseq_underconstruction.pdf $@ 
# done.
