userGuide: userGuide.Rnw
	echo "Stangle(\"userGuide.Rnw\")" | R --no-save --no-restore;	
	R --vanilla --slave < userGuide.R
	echo "Sweave(\"userGuide.Rnw\")" | R --no-save --no-restore
	# texi2dvi --pdf --clean --quiet userGuide.tex
	pdflatex userGuide.tex
