#!/bin/bash

logfile=$HOME/ricar/Scripts/regression/logs/log-`date +%Y-%m-%d`
function run {
    cd ..
    cp Modules/Conf-petites-matrices.pm Modules/Conf.pm  >> $logfile 2>&1
    ./GenNaturalDocs.pl -y                               >> $logfile 2>&1
    cd regression
    leavemealone start "tests de non regression PaStiX"  >> $logfile
    echo "people on vulcain"                             >> $logfile 2>&1
    who                                                  >> $logfile 2>&1
    echo "people on hephaistos"                          >> $logfile 2>&1
    ssh hephaistos "leavemealone start \"tests de non regression PaStiX\"" >> $logfile
    ssh hephaistos who                                   >> $logfile 2>&1
    cp Conf/conf.pastix.py Conf/conf.py                  >> $logfile 2>&1
    python regression.py --export --username=lacoste     >> $logfile 2>&1 
    python sendMail.py                                   >> $logfile 2>&1
    LANG=en_US.UTF-8 svn info svn+ssh://lacoste@scm.gforge.inria.fr/svn/ricar | grep Revision: |cut -c11- > version.txt
    leavemealone stop                                    >> $logfile
    ssh hephaistos "leavemealone stop"                   >> $logfile
}

# Check if the machine isn't locked
date                                                     > $logfile 2>&1
leavemealone status | grep "is locked" > /dev/null 2>&1
if [ $? -eq 0 ]; then
    leavemealone status                                  >> $logfile 2>&1
else
    cd ~/ricar/Scripts/regression                        >> $logfile 2>&1
    svn up ../..                                         >> $logfile 2>&1
    ssh hephaistos "leavemealone status" | grep "is locked" > /dev/null 2>&1
    if [ $? -eq 0 ]; then
	ssh hephaistos "leavemealone status" 
    else
        # check if the version file exists
	if [ -f version.txt ]; then
	    #check if the last version was the same
	    echo "comparaison de la version"                 >> $logfile 2>&1
	    cat version.txt | grep `LANG=en_US.UTF-8  svn info svn+ssh://lacoste@scm.gforge.inria.fr/svn/ricar | grep Revision: |cut -c11-` > /dev/null 2>&1
	    if [ $? -eq 0 ]; then
		echo "Version déjà testée"                   >> $logfile 2>&1
	    else
		run
	    fi
	else
	    run
	fi
    fi
fi