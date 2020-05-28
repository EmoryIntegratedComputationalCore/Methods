# -*- coding: utf-8 -*- 
#
#Author: JR
#Maintainer: JR
#License:  2019, EICC, GPL v2 or later
# 
#Emory/EICC/Methods/ProteomicsMethods/Proteus/Makefile
#----------------------------------------------------
#

HERE := $(shell Rscript --vanilla -e 'write(here::here(), stdout())')

.PHONY:	all clean

all: output/Proteus.pdf

clean:
					-rm -r output/*
					
output/Proteus.pdf: \
					src/Proteus.Rmd
			Rscript --vanilla -e "rmarkdown::render('$<')"
			mv src/Proteus.pdf $@
# done.
