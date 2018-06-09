function bv_ispc_initialize
{
    export DO_ISPC="no"
}

function bv_ispc_enable
{
    DO_ISPC="yes"
}

function bv_ispc_disable
{
    DO_ISPC="no"
}

function bv_ispc_depends_on
{
    echo ""
}

function bv_ispc_info
{
    export ISPC_VERSION=${ISPC_VERSION:-"1.9.2"}
    if [[ "$OPSYS" == "Darwin" ]] ; then
        export ISPC_FILE=${ISPC_FILE:-"ispc-v${ISPC_VERSION}-osx.tar.gz"}
        export ISPC_URL=${ISPC_URL:-"http://sdvis.org/ospray/download/dependencies/osx/"}
    else
        export ISPC_FILE=${ISPC_FILE:-"ispc-v${ISPC_VERSION}-linux.tar.gz"}
        export ISPC_URL=${ISPC_URL:-"http://sdvis.org/ospray/download/dependencies/linux/"}
    fi
    export ISPC_COMPATIBILITY_VERSION=${ISPC_COMPATIBILITY_VERSION:-"${ISPC_VERSION}"}
    export ISPC_BUILD_DIR=${ISPC_BUILD_DIR:-"${ISPC_VERSION}"}
    export ISPC_INSTALL_DIR_NAME=ispc-v$ISPC_VERSION-linux
    export ISPC_MD5_CHECKSUM=""
    export ISPC_SHA256_CHECKSUM=""
}

function bv_ispc_print
{
    printf "%s%s\n" "ISPC_FILE=" "${ISPC_FILE}"
    printf "%s%s\n" "ISPC_VERSION=" "${ISPC_VERSION}"
    printf "%s%s\n" "ISPC_COMPATIBILITY_VERSION=" "${ISPC_COMPATIBILITY_VERSION}"
    printf "%s%s\n" "ISPC_BUILD_DIR=" "${ISPC_BUILD_DIR}"
}

function bv_ispc_host_profile
{
    if [[ "$DO_ISPC" == "yes" ]] ; then
        echo >> $HOSTCONF
        echo "##" >> $HOSTCONF
        echo "## ISPC" >> $HOSTCONF
        echo "##" >> $HOSTCONF
        echo "VISIT_OPTION_DEFAULT(VISIT_ISPC_ROOT \${VISITHOME}/ispc/\${VISITARCH}/$ISPC_INSTALL_DIR_NAME)" \
            >> $HOSTCONF
    fi
}

function bv_ispc_print_usage
{
    #ispc does not have an option, it is only dependent on ispc.
    printf "%-15s %s [%s]\n" "--ispc" "Build ISPC" "$DO_ISPC"
}

function bv_ispc_ensure
{
    if [[ "$DO_ISPC" == "yes" ]] ; then
        ensure_built_or_ready "ispc" $ISPC_VERSION $ISPC_BUILD_DIR $ISPC_FILE $ISPC_URL
        if [[ $? != 0 ]] ; then
            ANY_ERRORS="yes"
            DO_ISPC="no"
            error "Unable to build ISPC.  ${ISPC_FILE} not found."
        fi
    fi
}

function bv_ispc_dry_run
{
    if [[ "$DO_ISPC" == "yes" ]] ; then
        echo "Dry run option not set for ISPC."
    fi
}

# ***************************************************************************
# build_ispc
#
# Modifications:
#
# ***************************************************************************

function build_ispc
{
    # Unzip the ISPC tarball and copy it to the VisIt installation.
    info "Installing prebuilt ISPC"    
    tar zxvf $ISPC_FILE
    mkdir $VISITDIR/ispc
    mkdir $VISITDIR/ispc/$VISITARCH
    cp -R $ISPC_INSTALL_DIR_NAME "$VISITDIR/ispc/$VISITARCH"
    rm -rf $ISPC_INSTALL_DIR_NAME
    if [[ "$DO_GROUP" == "yes" ]] ; then
        chmod -R ug+w,a+rX "$VISITDIR/ispc/$VISITARCH"
        chgrp -R ${GROUP} "$VISITDIR/ispc/$VISITARCH"
    fi
    cd "$START_DIR"
    info "Done with ISPC"
    return 0
}

function bv_ispc_is_enabled
{
    if [[ $DO_ISPC == "yes" ]]; then
        return 1    
    fi
    return 0
}

function bv_ispc_is_installed
{
    check_if_installed "ispc"
    if [[ $? == 0 ]] ; then
        return 1
    fi
    return 0
}

function bv_ispc_build
{
    if [[ "$DO_ISPC" == "yes" ]] ; then
        check_if_installed "ispc"
        if [[ $? == 0 ]] ; then
            info "Skipping build of ISPC"
        else
            build_ispc
            if [[ $? != 0 ]] ; then
                error "Unable to build or install ISPC.  Bailing out."
            fi
            info "Done building ISPC"
        fi
    fi
}

