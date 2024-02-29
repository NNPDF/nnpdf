# Minimal makefile for Sphinx documentation
#

# You can set these variables from the command line.
SPHINXOPTS    =
SPHINXBUILD   = sphinx-build
SOURCEDIR     = source
BUILDDIR      = build

# Put it first so that "make" without argument is like "make help".
help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: help Makefile

# Catch-all target: route all unknown targets to Sphinx using the new
# "make mode" option.  $(O) is meant as a shortcut for $(SPHINXOPTS).
%: Makefile MakeTheoryCSV
	@if test $@ != "clean"; then \
		sphinx-apidoc -f -o ./$(SOURCEDIR)/modules/validphys ../../validphys2/src/validphys/ ; \
		sphinx-apidoc -f -o ./$(SOURCEDIR)/modules/n3fit-code ../../n3fit/src/n3fit/ ; \
	fi
	@$(SPHINXBUILD) -M $@ "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

.PHONY: rsync
rsync: html
	rsync -avzr ./$(BUILD)/html/ nnpdf@vp.nnpdf.science:WEB/sphinx-docs/

view: html
	$(BROWSER) $(BUILDDIR)/html/index.html

MakeTheoryCSV:
	sqlite3 -header -csv ../../validphys2/src/validphys/datafiles/theory.db "select * from TheoryIndex;" > source/theories.csv
