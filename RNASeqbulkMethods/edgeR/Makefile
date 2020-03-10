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

all: output/DEGAnalysiswithedgeR.pdf

clean:
	
					-rm -r output/*
		
output/DEGAnalysiswithedgeR.pdf: \
					src/DEGAnalysiswithedgeR.Rmd
			-rm src/DEGAnalysiswithedgeR.pdf
			Rscript --vanilla -e "rmarkdown::render('$<')"
			mv src/DEGAnalysiswithedgeR.pdf $@
# done.
