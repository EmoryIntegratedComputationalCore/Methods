# -*- coding: utf-8 -*- 
#
#Author: JR
#Maintainer: JR
#License:  2019, EICC, GPL v2 or later
# 
#Emory/EICC/Methods/Miscellaneous/qvalue/Makefile
#----------------------------------------------------
#

HERE := $(shell Rscript --vanilla -e 'write(here::here(), stdout())')

.PHONY:	all clean

all: output/qvalue.pdf

clean:
	
					-rm -r output/*
		
output/qvalue.pdf: \
					src/qvalue.Rmd
			-rm output/qvalue.pdf
			Rscript --vanilla -e "rmarkdown::render('$<')"
			mv src/qvalue.pdf $@
# done.
