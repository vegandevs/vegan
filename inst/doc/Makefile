all: FAQ-vegan.pdf decision-vegan.pdf intro-vegan.pdf diversity-vegan.pdf NEWS.html
FAQ-vegan.pdf: FAQ-vegan.texi
	texi2dvi --pdf --clean FAQ-vegan.texi
decision-vegan.pdf: decision-vegan.tex
	texi2dvi --pdf --clean decision-vegan.tex
	-rm -f decision-vegan-0*.* Rplots.*
intro-vegan.pdf: intro-vegan.tex
	texi2dvi --pdf --clean intro-vegan.tex
	-rm -f intro-vegan-0*.* Rplots.*
diversity-vegan.pdf: diversity-vegan.tex
	texi2dvi --pdf --clean diversity-vegan.tex
	-rm -f diversity-vegan-0*.* Rplots.*
NEWS.html: ../NEWS.Rd
	"$(R_HOME)/bin/R" CMD Rd2txt -t html ../NEWS.Rd -o NEWS.html

