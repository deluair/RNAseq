# Makefile for compiling LaTeX manuscript

# Main file name (without .tex extension)
MAIN = soybean_drought_rnaseq

# LaTeX compiler
LATEX = pdflatex
BIBTEX = bibtex

# Default target
all: $(MAIN).pdf

# Compile PDF
$(MAIN).pdf: $(MAIN).tex references.bib
	$(LATEX) $(MAIN).tex
	$(BIBTEX) $(MAIN)
	$(LATEX) $(MAIN).tex
	$(LATEX) $(MAIN).tex

# Clean auxiliary files
clean:
	rm -f *.aux *.bbl *.blg *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz

# Clean all generated files including PDF
distclean: clean
	rm -f $(MAIN).pdf

# View PDF (macOS)
view: $(MAIN).pdf
	open $(MAIN).pdf

# Continuous compilation with latexmk
watch:
	latexmk -pdf -pvc $(MAIN).tex

.PHONY: all clean distclean view watch 