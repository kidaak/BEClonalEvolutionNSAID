TOP=../..
PKG=${shell cd ${TOP};pwd}
SUITE=doRUnit.R
#R=${R_HOME}/bin/R

all: inst test

inst: # Install package

	cd ${TOP}/..;\
	R CMD INSTALL ${PKG} --library=R/i686-pc-linux-gnu-library/2.12

test: # Run unit tests
	export RCMDCHECK=FALSE;\
	export RUNITFILEPATTERN="$(file)";\
	cd ${TOP}/tests;\
	R --vanilla --slave < ${SUITE}

