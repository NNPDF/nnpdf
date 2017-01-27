%.html: %.md
	pandoc -f markdown-markdown_in_html_blocks+raw_html --to html5 -s -S --toc $< --template report.html5 --css report.css  -o $@


%.png: %.dot
	dot -Tpng $< -o  $@

.PHONY: all

DOTS = $(wildcard *.dot)
MDS = $(wildcard *.md)
all: $(patsubst %.md, %.html, $(MDS)) $(patsubst %.dot, %.png, $(DOTS))


.PHONY: rsync
rsync:
	rsync -avz ./ nnpdf:WEB/validphys-docs/
