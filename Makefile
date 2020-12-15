# Makefile for "CNCP Playground"
#
# Copyright (C) 2020 Simon Dobson
#
# This file is part of cncp-playground, an interactive view of network science 
#
# cncp-playground is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# cncp-playground is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with cncp-playground. If not, see <http://www.gnu.org/licenses/gpl.html>.

PROJECT = cncp-playground


# ----- Sources -----

# Content
INDEX = index.md
TEXT =
NOTEBOOKS =

# Image files
RAW_IMAGES = \
	cc-by-nc-sa.png \
	sd.png
SVG_IMAGES =
GENERATED_IMAGES = $(SVG_IMAGES:.svg=.png)

# Bibliography
BIBLIOGRAPHY = bibliography.bib

# License
LICENSE = LICENSE

# All content
CONTENT = \
	$(INDEX) \
	$(TEXT) \
	$(NOTEBOOKS) \
	$(RAW_IMAGES) $(SVG_IMAGES) $(GENERATED_IMAGES) \
	$(BIBLIOGRAPHY)


# ----- Tools -----

# Root directory
ROOT = $(shell pwd)

# Base commands
PYTHON = python3
IPYTHON = ipython
JUPYTER = jupyter
RSYNC = rsync
PIP = pip
VIRTUALENV = $(PYTHON) -m venv
ACTIVATE = . $(VENV)/bin/activate
INKSCAPE = inkscape
CONVERT = convert
TR = tr
CAT = cat
SED = sed
RM = rm -fr
CP = cp
CHDIR = cd
MKDIR = mkdir -p
ZIP = zip -r
UNZIP = unzip
WGET = wget
ECHO = echo

# Datestamp
DATE = `date`

# Requirements and venv
VENV = venv3
REQUIREMENTS = requirements.txt

# Constructed commands
RUN_SERVER = PYTHONPATH=. $(JUPYTER) lab


# ----- Top-level targets -----

# Default prints a help message
help:
	@make usage

# Run the notebook server
live: env
	$(ACTIVATE) && $(RUN_SERVER)

# Build a development venv
.PHONY: env
env: $(VENV)

$(VENV):
	$(VIRTUALENV) $(VENV)
	$(ACTIVATE) && $(PIP) install wheel && $(PIP) install -r requirements.txt


# Clean up the build
clean:
	$(RM) $(BOOK_DIR) $(GENERATED_IMAGES)

# Clean up everything, including the venv (which is quite expensive to rebuild)
reallyclean: clean
	$(RM) $(VENV)


# ----- Construction rules -----

.SUFFIXES: .svg .pdf .png .ipynb .html .md .rst

.svg.pdf:
	$(INKSCAPE) $*.svg --export-pdf=$*.pdf

.svg.png:
	$(CONVERT) $*.svg $*.png


# ----- Usage -----

define HELP_MESSAGE
Editing:
   make live         run the notebook server

Maintenance:
   make env          create a virtual environment
   make clean        clean-up the build
   make reallyclean  delete the venv as well

endef
export HELP_MESSAGE

usage:
	@echo "$$HELP_MESSAGE"


