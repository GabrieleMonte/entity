# Build the DPDM reproduction document.
#   module load texlive/2023 && make
MAIN = main
all: $(MAIN).pdf
$(MAIN).pdf: $(MAIN).tex
	latexmk -pdf -interaction=nonstopmode -halt-on-error $(MAIN).tex
clean:
	latexmk -C
	rm -f *.aux *.log *.out *.toc *.fls *.fdb_latexmk
.PHONY: all clean
