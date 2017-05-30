#!/bin/bash

VISITSRC=${1}
WORKSRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

RELPATH() {
    python -c "import os.path; print os.path.relpath('${1}', '$WORKSRC/src')"
}

COPYFILE() {
    install -Dv ${VISITSRC}/$1 ${WORKSRC}/src-ref/$1
}

READFILES() {
    for i in "$1"/*; do
	path=$(RELPATH $i)
	if [ -d "$i" ];then
            echo "Entering directory: $path"
            READFILES "$i"
	elif [ -f "$i" ]; then
            echo " -- Processing: $path"
	    COPYFILE $path
	fi
    done
}

# adding files with same filename
READFILES $WORKSRC/src

# removed files
echo "Processing removed files"
for i in $(cat removelist.txt); do 
    echo " -- Processing: $i" 
    COPYFILE $i
done
