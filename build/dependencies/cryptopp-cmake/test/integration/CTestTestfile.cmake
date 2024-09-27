# CMake generated Testfile for 
# Source directory: /Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration
# Build directory: /Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/dependencies/cryptopp-cmake/test/integration
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test([=[cryptopp-int-install-default-configure]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "-GNinja" "-D" "TEST_CMAKE_FILES_DIR=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/cmake" "-D" "INT_TEST_CMAKE_FILES_DIR=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/cmake" "-D" "TEST_EXAMPLE_SOURCES_DIR=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/example-src" "-D" "USE_CCACHE=TRUE" "-D" "CMAKE_USER_MAKE_RULES_OVERRIDE=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/cmake/ResetInitialCompilerOptions.cmake" "-D" "CMAKE_INSTALL_PREFIX=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/install/int-install-default" "-S" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/int-install-default" "-B" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/int-install-default" "-D" "CPM_cryptopp-cmake_SOURCE=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/../.." "-D" "CMAKE_VERBOSE_MAKEFILE=ON" "-D" "CRYPTOPP_MINIMUM_CMAKE_VERSION=3.20" "-D" "CMAKE_BUILD_TYPE=Debug")
set_tests_properties([=[cryptopp-int-install-default-configure]=] PROPERTIES  FIXTURES_SETUP "int-install-default-config" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-default" install_dir "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/install/int-install-default" _BACKTRACE_TRIPLES "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;21;add_test;/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-default-build]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/int-install-default" "--config" "Debug")
set_tests_properties([=[cryptopp-int-install-default-build]=] PROPERTIES  FIXTURES_REQUIRED "int-install-default-config" FIXTURES_SETUP "int-install-default-build" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-default" _BACKTRACE_TRIPLES "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;55;add_test;/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-default-install]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/int-install-default" "--target" "install" "--config" "Debug")
set_tests_properties([=[cryptopp-int-install-default-install]=] PROPERTIES  FIXTURES_REQUIRED "int-install-default-build" FIXTURES_SETUP "int-install-default-install" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-default" _BACKTRACE_TRIPLES "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;65;add_test;/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-default-checks]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/int-install-default" "--target" "do-checks" "--config" "Debug")
set_tests_properties([=[cryptopp-int-install-default-checks]=] PROPERTIES  FIXTURES_REQUIRED "int-install-default-install" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-default" _BACKTRACE_TRIPLES "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;76;add_test;/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-prefix-configure]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "-GNinja" "-D" "TEST_CMAKE_FILES_DIR=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/cmake" "-D" "INT_TEST_CMAKE_FILES_DIR=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/cmake" "-D" "TEST_EXAMPLE_SOURCES_DIR=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/example-src" "-D" "USE_CCACHE=TRUE" "-D" "CMAKE_USER_MAKE_RULES_OVERRIDE=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/cmake/ResetInitialCompilerOptions.cmake" "-D" "CMAKE_INSTALL_PREFIX=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/install/int-install-prefix" "-S" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/int-install-prefix" "-B" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/int-install-prefix" "-D" "CPM_cryptopp-cmake_SOURCE=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/../.." "-D" "CMAKE_VERBOSE_MAKEFILE=ON" "-D" "CRYPTOPP_MINIMUM_CMAKE_VERSION=3.20" "-D" "CMAKE_BUILD_TYPE=Debug")
set_tests_properties([=[cryptopp-int-install-prefix-configure]=] PROPERTIES  FIXTURES_SETUP "int-install-prefix-config" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-prefix" install_dir "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/install/int-install-prefix" _BACKTRACE_TRIPLES "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;21;add_test;/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-prefix-build]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/int-install-prefix" "--config" "Debug")
set_tests_properties([=[cryptopp-int-install-prefix-build]=] PROPERTIES  FIXTURES_REQUIRED "int-install-prefix-config" FIXTURES_SETUP "int-install-prefix-build" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-prefix" _BACKTRACE_TRIPLES "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;55;add_test;/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-prefix-install]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/int-install-prefix" "--target" "install" "--config" "Debug")
set_tests_properties([=[cryptopp-int-install-prefix-install]=] PROPERTIES  FIXTURES_REQUIRED "int-install-prefix-build" FIXTURES_SETUP "int-install-prefix-install" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-prefix" _BACKTRACE_TRIPLES "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;65;add_test;/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-install-prefix-checks]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/int-install-prefix" "--target" "do-checks" "--config" "Debug")
set_tests_properties([=[cryptopp-int-install-prefix-checks]=] PROPERTIES  FIXTURES_REQUIRED "int-install-prefix-install" LABELS "cryptopp;cryptopp-integration-tests;cryptopp-install-prefix" _BACKTRACE_TRIPLES "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;76;add_test;/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-find-package-configure]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "-GNinja" "-D" "TEST_EXAMPLE_SOURCES_DIR=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/example-src" "-D" "CRYPTOPP_CMAKE_INSTALL_ROOT=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/install" "-S" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/int-find-package" "-B" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/int-find-package" "-D" "CMAKE_VERBOSE_MAKEFILE=ON" "-D" "CRYPTOPP_MINIMUM_CMAKE_VERSION=3.20" "-D" "CMAKE_BUILD_TYPE=Debug" "-D" "cryptopp_DIR=/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/install/int-install-default/share/cmake/cryptopp")
set_tests_properties([=[cryptopp-int-find-package-configure]=] PROPERTIES  FIXTURES_REQUIRED "int-install-default-install" FIXTURES_SETUP "int-find-package-config" LABELS "cryptopp;cryptopp-find_package" _BACKTRACE_TRIPLES "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;95;add_test;/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")
add_test([=[cryptopp-int-find-package-build]=] "/opt/homebrew/Cellar/cmake/3.27.7/bin/cmake" "--build" "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/build/test-dirs/int-find-package" "--config" "Debug")
set_tests_properties([=[cryptopp-int-find-package-build]=] PROPERTIES  FIXTURES_REQUIRED "int-find-package-config" LABELS "cryptopp;cryptopp-find_package" _BACKTRACE_TRIPLES "/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;122;add_test;/Users/dulumrae/Desktop/Kübra Hoca veya Ağ Güvenliği/Almancılar/Cecilia/CECILIA-Simple/dependencies/cryptopp-cmake/test/integration/CMakeLists.txt;0;")