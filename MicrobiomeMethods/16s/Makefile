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

all: output/16s.pdf

clean:
	-rm -r output/*
		
output/16s.pdf: \
					src/16s.Rmd
					-rm src/16s.pdf
			Rscript --vanilla -e "rmarkdown::render('$<')"
			mv src/16s.pdf $@
# done.
