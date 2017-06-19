#!/bin/bash

convert2sed()
{
    echo $(echo $1 | sed 's_/_\\/_g')
}

declare -a plibs=(
    '/home/sci/qwu/software/Kepler/VisIt/visit/vtk/6.1.0/linux-x86_64_gcc-4.8/lib'
    '/home/sci/qwu/software/Kepler/VisIt/visit/llvm/4.0.0/linux-x86_64_gcc-4.8/lib'
)
declare -a pbins=(
    '/home/sci/qwu/software/Kepler/VisIt/visit/llvm/4.0.0/linux-x86_64_gcc-4.8/bin'
)

declare -a pincs=(
    '/home/sci/qwu/software/Kepler/VisIt/visit/llvm/4.0.0/linux-x86_64_gcc-4.8/include'
)

if [[ "$1" == "activate" ]]; then

    ## now loop through the above array
    for i in "${plibs[@]}"
    do
	export LIBRARY_PATH=$i:${LIBRARY_PATH}
	export LD_LIBRARY_PATH=$i:${LD_LIBRARY_PATH}
    done
    for i in "${pbins[@]}"
    do
	export PATH=$i:${PATH}
    done
    for i in "${pincs[@]}"
    do
	export CPATH=$i:${CPATH}
    done

elif [[ "$1" == "deactivate" ]]; then

    ## now loop through the above array
    for i in "${plibs[@]}"
    do
	ri=$(convert2sed $i)
	export LIBRARY_PATH=$(echo "$LIBRARY_PATH" | sed -e 's/'$ri'://g')
	export LD_LIBRARY_PATH=$(echo "$LD_LIBRARY_PATH" | sed -e 's/'$ri'://g')
    done
    ## now loop through the above array
    for i in "${pbins[@]}"
    do
	ri=$(convert2sed $i)
	export PATH=$(echo "$PATH" | sed -e 's/'$ri'://g')
    done
    ## now loop through the above array
    for i in "${pincs[@]}"
    do
	ri=$(convert2sed $i)
	export CPATH=$(echo "$CPATH" | sed -e 's/'$ri'://g')
    done

fi
