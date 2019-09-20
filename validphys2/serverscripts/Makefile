.PHONY: all
all:
	git ls-files | rsync -avz --exclude 'Makefile' --exclude '.gitignore'  --files-from=- . nnpdf:
