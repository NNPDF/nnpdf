all: tutorial.md
	pandoc -f markdown-markdown_in_html_blocks+raw_html --to html5 -s -S --toc tutorial.md --css report.css  -o guide.html
	pandoc -f markdown-markdown_in_html_blocks+raw_html --to html5 -s -S --toc plotting_format.md --css report.css  -o plotting_format.html

