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
	--exclude '.*'      \
	${WORKSRC}/src/* ${VISITSRC}

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
	--exclude '.*'      \
	${WORKSRC}/src/* ${VISITSRC}

elif [ "$3" == "clean" ]; then 

    rm -r ${VISITSRC}/avt/Plotter/OSPRay
    rm ${VISITSRC}/avt/Filters/imgMetaData.*
    rm ${VISITSRC}/avt/Filters/avtImgCommunicator.*

fi
