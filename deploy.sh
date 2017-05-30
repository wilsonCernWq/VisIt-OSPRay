#!/bin/bash

VISITSRC=${1}
OPTION=${2}
WORKSRC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

HELP()
{
    echo "Usage: bash deploy.sh <destination-dir> <option>"
}

if [ "${VISITSRC}" == "" ]; then
    HELP
    exit
fi

if [ "$OPTION" == "source" ]; then 

    rsync -vr \
	--exclude '*.txt'   \
	--exclude '*.cmake' \
	--exclude '*.in'    \
	--exclude '.*'      \
	${WORKSRC}/src/* ${VISITSRC}

elif [ "$OPTION" == "cmake" ]; then 

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
	--exclude '.*'      \
	${WORKSRC}/src/* ${VISITSRC}

elif [ "$OPTION" == "clean" ]; then 

    rm -r ${VISITSRC}/avt/Plotter/OSPRay
    rm ${VISITSRC}/avt/Filters/imgMetaData.*
    rm ${VISITSRC}/avt/Filters/avtImgCommunicator.*

fi
