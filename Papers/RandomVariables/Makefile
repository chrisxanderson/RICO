SHELL = /bin/sh

PROJECT  = RandomVariables
TARGET   = $(PROJECT).pdf

IMAGES   = ./Images
DEPS     = *.tex $(patsubst %.svg,%.eps,$(wildcard $(IMAGES)/*.svg))
VPATH    = $(IMAGES)

LATEX    = latex
BIBTEX   = bibtex
DVIPDF   = dvipdf
INKSCAPE = inkscape 
RM       = rm

all: $(TARGET)

$(TARGET) : $(PROJECT).dvi
	$(DVIPDF) $(PROJECT).dvi

$(PROJECT).dvi : $(PROJECT).bbl $(DEPS)
	$(LATEX) $(PROJECT)
	$(LATEX) $(PROJECT)

$(PROJECT).bbl : $(PROJECT).bib $(DEPS)
	@echo "in the .bbl rule"
	@echo '$$(DEPS) = $(DEPS)'
	$(LATEX) $(PROJECT)
	$(BIBTEX) $(PROJECT)

$(IMAGES)/%.eps : $(IMAGES)/%.svg
	@echo "in the .eps rule '$@' and '$<'"
	$(INKSCAPE) --export-eps=$@ $<

clean :
	$(RM) -rf *.aux *.dvi *.pdf *.log 
	$(RM) -rf *.blg *.bbl
	$(RM) -rf $(IMAGES)/*.eps


