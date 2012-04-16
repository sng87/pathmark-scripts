THISDIR = $(CURDIR)

init.sh :
	echo \
	export PATH=$(THISDIR)/bin:\$${PATH} > init.sh
	echo \
	export PYTHONPATH=$(THISDIR)/bin:\$${PYTHONPATH} >> init.sh
	echo \
	setenv PATH $(THISDIR)/bin:\$${PATH} > init.csh
	echo \
	setenv PYTHONPATH $(THISDIR)/bin:\$${PYTHONPATH} >> init.csh

clean :
	rm -f init.sh init.csh
