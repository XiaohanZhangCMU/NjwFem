#!/bin/sh

if [ -x "`which svnversion`" ]; then
    rev=`LC_ALL=C svnversion | sed -e 's/ //'`
    if [ "$rev" != 'exporté' ]
    then
	echo $rev;
	exit;
    fi;
fi

if [ -f Revision ]
then 
    cat ./Revision
else
    if [ -f ../../Revision ]
    then
	cat ../../Revision
    else
	echo "NoRevisionFound"
    fi
fi;
