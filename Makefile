SHELL := /bin/bash
VCFS := $(wildcard VCF/*.vcf)
FASTAS := $(patsubst VCF/%.vcf,fasta/%.fasta,$(VCFS))

all: $(FASTAS)

fasta/%.fasta: VCF/%.vcf
	python3 V2F.py $< $@

help:
		@echo "Available commands:"
		@echo "  make env     - Setup the python virtual environment."
		@echo "  make install - Install dependencies."
		@echo "  make data    - Run data preprocessing script."
		@echo "  make clear   - Remove VCF directory."
		@echo "  make run     - Run the main script."
		@echo "  make clean   - Remove virtual environment and other cleanup."

env:
		python3 -m venv --system-site-packages venv
		# Activate the virtual environment with:
		# source venv/bin/activate
		# And then run setenv.sh manually with:
		# . setenv.sh

install:
		python3 -m pip install --upgrade pip
		pip install -r requirements.txt

data:
		python3 data.py

clear:
		rm -fr VCF

run:
		python3 main.py

clean:
		rm -fr venv
		rm -fr fasta/*.fasta


