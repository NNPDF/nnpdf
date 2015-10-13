all: electrodynamics.tex
	latex electrodynamics.tex
	bibtex electrodynamics.aux
	latex electrodynamics.tex
	dvipdf electrodynamics.dvi
	rm -rf *.out *.toc *.log *.aux *.bbl *.blg

clean:
	rm -rf *.out *.toc *.log *.aux *.bbl *.blg *.pdf
