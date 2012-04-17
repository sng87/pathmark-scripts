THISDIR = $(CURDIR)

init.sh :
	echo \
	export PATH=$(THISDIR)/bin:\$${PATH} > init.sh
	echo \
	if [ -n "\$${PYTHONPATH+x}" ] >> init.sh
	echo \
	then >> init.sh
	echo \
	  export PYTHONPATH=$(THISDIR)/bin:\$${PYTHONPATH} >> init.sh
	echo \
	else >> init.sh
	echo \
	  export PYTHONPATH=$(THISDIR)/bin >> init.sh
	echo \
	fi >> init.sh
	echo \
	setenv PATH $(THISDIR)/bin:\$${PATH} > init.csh
	echo \
	if \$$?PYTHONPATH then >> init.csh
	echo \
	  setenv PYTHONPATH $(THISDIR)/bin:\$${PYTHONPATH} >> init.csh
	echo \
	else >> init.csh
	echo \
	  setenv PYTHONPATH $(THISDIR)/bin >> init.csh
	echo \
	endif >> init.csh

clean :
	rm -f init.sh init.csh
