DOCUMENT=ipred-examples

TEXT= \
	$(DOCUMENT).Rnw \
	$(DONE)

BIBS= \
	ipred.bib \
	$(DONE) 

all	: pdf
	make clean	

ps	: $(DOCUMENT).ps

pdf	: $(DOCUMENT).pdf

$(DOCUMENT).ps	:  $(DOCUMENT).dvi ; dvips -Z -Pcmz -Pamz $(DOCUMENT) -o $(DOCUMENT).ps

$(DOCUMENT).tex	: $(DOCUMENT).Rnw
	echo "library(tools); Sweave(\"$(DOCUMENT).Rnw\");" | R --no-save --no-restore

$(DOCUMENT).dvi	: $(TEXT) $(BIBS)
	latex $(DOCUMENT) && \
	bibtex $(DOCUMENT) && \
	latex $(DOCUMENT) && \
	latex $(DOCUMENT) && \
	latex $(DOCUMENT)

$(DOCUMENT).pdf	: $(DOCUMENT).dvi
	pdflatex  $(DOCUMENT)

clean   :
	-rm -f `ls  $(DOCUMENT).* | sed -e /\.Rnw/d | sed -e /\.pdf/d | sed -e /\.tex/d`
	-rm -f *~
