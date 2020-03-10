# -*- coding: utf-8 -*- 
#
#Author: JR
#Maintainer: JR
#License:  2019, EICC, GPL v2 or later
# 
#Emory/EICC/Methods/MicrobiomeMethods/16s/Makefile
#----------------------------------------------------
#

HERE := $(shell Rscript --vanilla -e 'write(here::here(), stdout())')

.PHONY:	all clean

all: output/16sAnalysiswithphyloseq_LDM.pdf

clean:
	
					-rm -r output/*
		
output/16sAnalysiswithphyloseq_LDM.pdf: \
					src/16sAnalysiswithphyloseq_LDM.Rmd
			-rm src/16sAnalysiswithphyloseq_LDM.pdf
			Rscript --vanilla -e "rmarkdown::render('$<')"
			mv src/16sAnalysiswithphyloseq_LDM.pdf $@
# done.
