Components for point clouds and depth images using pcl

## Add Remotes

    $ conan remote add camposs "https://conan.campar.in.tum.de/api/conan/conan-camposs"
    $ conan remote add ubitrack "https://conan.campar.in.tum.de/api/conan/conan-ubitrack"

## For Users: Use this package

### Basic setup
# PCL needs in windows
conan remote add sight https://conan.ircad.fr/artifactory/api/conan/sight

$ mkdir build
$ cd build
$ conan install ..

# to use the example direcory
$ conan build .. -pf ../example

$ cd example
$ conan install conan/conanfile.py



