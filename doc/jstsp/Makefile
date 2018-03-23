TEXOPT = -interaction nonstopmode -file-line-error-style
#all: $(TEX_SRC).pdf
TEX_SRC = $(wildcard *.tex)
TARGETS := $(patsubst %.tex,%.pdf,$(TEX_SRC))
TEXIN := .:./inc:
TEXFO := ./inc/tfm/::

all: $(TARGETS)


%.pdf: %.tex
	TEXFONTS=$(TEXFO) TEXINPUTS=$(TEXIN) pdflatex $(TEXOPT) -draftmode $<
	TEXFONTS=$(TEXFO) TEXINPUTS=$(TEXIN) pdflatex $(TEXOPT) -draftmode -interaction batchmode $<
	TEXFONTS=$(TEXFO) TEXINPUTS=$(TEXIN) pdflatex $(TEXOPT) -interaction batchmode $<


clean:
	rm -f *.aux *.log *.out *.ps *.toc *.nav *.snm *.dvi *.bbl *.blg *.idx *.glo *.ist *~


purge: clean
	-rm -f *.pdf

