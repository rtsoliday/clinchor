
all:: clinchor.ps clinchor.pdf html

clinchor.ps: clinchor.tex
	rm -f clinchor.aux clinchor.dvi clinchor.log clinchor.toc
	latex clinchor
	latex clinchor
	dvips clinchor -o clinchor.ps

clinchor.pdf: clinchor.ps
	ps2pdf13 clinchor.ps clinchor.pdf

html:
	mkdir -p clinchor
	htlatex clinchor "html,3" "" "-dclinchor/"
	rm -f *.html *.png *.css *.4ct *.tmp *.dvi *.idv *.4tc *.lg *.xref *.log *.aux 

clean::
	rm -f clinchor.aux clinchor.dvi clinchor.log clinchor.toc clinchor.ps clinchor.pdf
	rm -rf clinchor
