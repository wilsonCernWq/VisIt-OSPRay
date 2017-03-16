## ======================================================================================= ##
## Copyright 2014-2015 Texas Advanced Computing Center, The University of Texas at Austin  ##
##                                                                                         ##
## Licensed under the BSD 3-Clause License, (the "License"); you may not use this file     ##
## except in compliance with the License.                                                  ##
## A copy of the License is included with this software in the file LICENSE.               ##
## If your copy does not contain the License, you may obtain a copy of the License at:     ##
##                                                                                         ##
##     http://opensource.org/licenses/BSD-3-Clause                                         ##
##                                                                                         ##
## Unless required by applicable law or agreed to in writing, software distributed under   ##
## the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY ##
## KIND, either express or implied.                                                        ##
## See the License for the specific language governing permissions and limitations under   ##
## limitations under the License.                                                          ##
## ======================================================================================= ##
# Find OSPRay
# defines:
# OSPRAY_FOUND
# OSPRAY_INCLUDE_DIRS
# OSPRAY_LIBRARIES

find_package(ospray REQUIRED)
include_directories(${OSPRAY_INCLUDE_DIRS})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(OSPRAY 
  DEFAULT_MSG
  OSPRAY_INCLUDE_DIR
  OSPRAY_LIBRARIES
)

set(OSPRAY_LIBRARY_DIRS ${OSPRAY_INCLUDE_DIRS}/../lib64/)
