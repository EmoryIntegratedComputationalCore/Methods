# -*- coding: utf-8 -*- 
#
#Author: JR
#Maintainer: JR
#License:  2019, EICC, GPL v2 or later
# 
#Emory/EICC/Methods/Miscellaneous/PROPER/Makefile
#----------------------------------------------------
#

HERE := $(shell Rscript --vanilla -e 'write(here::here(), stdout())')

.PHONY:	all clean

all: output/PROPER.pdf

clean:
	
					-rm -r output/*
		
output/PROPER.pdf: \
					src/PROPER.Rmd
			-rm src/PROPER.pdf
			Rscript --vanilla -e "rmarkdown::render('$<')"
			mv src/PROPER.pdf $@
# done.
