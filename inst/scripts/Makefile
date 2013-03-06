userGuide: userGuide.Rnw
	echo "Stangle(\"userGuide.Rnw\")" | R --no-save --no-restore;	
	R --vanilla --slave < userGuide.R
	echo "Sweave(\"userGuide.Rnw\")" | R --no-save --no-restore
	pdflatex userGuide.tex
figure3Guide: figure3Guide.Rnw
	echo "Sweave(\"figure3Guide.Rnw\")" | R --no-save --no-restore
	pdflatex figure3Guide.tex
