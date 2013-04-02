DOCFILE=report
all: $(DOCFILE).pdf

$(DOCFILE).pdf: $(DOCFILE).tex
	pdflatex $(DOCFILE).tex
	bibtex   $(DOCFILE).aux
	pdflatex $(DOCFILE).tex
	pdflatex $(DOCFILE).tex

clean:
	rm -f *.aux *.out *.log *.pdf *.toc *.bbl *.blg
