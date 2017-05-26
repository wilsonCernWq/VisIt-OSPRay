#!/bin/bash

VISITSRC=${1}
WORKSRC=${2}

HELP()
{
    echo "Usage: bash deploy.sh <destination-dir> <source-dir>"
}

if [ "${VISITSRC}" == "" ]; then
    HELP
    exit
fi

if [ "${WORKSRC}" == "" ]; then
    HELP
    exit
fi

if [ "$3" == "source" ]; then 

    rsync -vr \
	--exclude '*.txt'   \
	--exclude '*.cmake' \
	--exclude '*.in'    \
	--exclude '*.sh'    \
	--exclude '.*'      \
	--exclude 'README*' \
	--exclude 'exclude' \
	${WORKSRC}/* ${VISITSRC}

elif [ "$3" == "cmake" ]; then 

    rsync -vr \
	--exclude '*.C'     \
	--exclude '*.h'     \
	--exclude '*.cpp'   \
	--exclude '*.hpp'   \
	--exclude '*.cxx'   \
	--exclude '*.hxx'   \
	--exclude '*.xml'   \
	--exclude '*.code'   \
	--exclude '*.java'   \
	--exclude '*.sh'    \
	--exclude '.*'      \
	--exclude 'README*' \
	--exclude 'exclude' \
	${WORKSRC}/* ${VISITSRC}

fi
