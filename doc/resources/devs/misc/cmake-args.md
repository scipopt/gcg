# CMake Arguments {#cmake-args}

|             Argument             | Value Type |             Default Value             |                                     Description                                      |
|:--------------------------------:|:----------:|:-------------------------------------:|:------------------------------------------------------------------------------------:|
|        BLISS_INCLUDE_DIR         |    PATH    |      BLISS_INCLUDE_DIR-NOTFOUND       |                                          -                                           |
|          BLISS_LIBRARY           |  FILEPATH  |        BLISS_LIBRARY-NOTFOUND         |                                          -                                           |
|          BUILD_TESTING           |    BOOL    |                  ON                   |                                          -                                           |
|             CLIQUER              |    BOOL    |                  ON                   |                                          -                                           |
|       CLIQUER_INCLUDE_DIR        |    PATH    |     CLIQUER_INCLUDE_DIR-NOTFOUND      |                                          -                                           |
|         CLIQUER_LIBRARY          |  FILEPATH  |       CLIQUER_LIBRARY-NOTFOUND        |                                          -                                           |
|         CMAKE_BUILD_TYPE         |   STRING   |              no default               |                                          -                                           |
|       CMAKE_INSTALL_PREFIX       |    PATH    |              /usr/local               |                                          -                                           |
| CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS |    BOOL    |                  ON                   |                                          -                                           |
|             COVERAGE             |    BOOL    |                  OFF                  |                                          -                                           |
|       COVERAGE_CTEST_ARGS        |   STRING   |              no default               |                                          -                                           |
|              CPLEX               |    BOOL    |                  OFF                  |                       Enable CPLEX solver for problem solving                        |
|             CXXONLY              |    BOOL    |                  OFF                  |                                          -                                           |
|             DEBUGSOL             |    BOOL    |                  OFF                  |                                          -                                           |
|             EXPRINT              |   STRING   |                 cppad                 |                                          -                                           |
|               GCG                |    BOOL    |                  ON                   |                                          -                                           |
|          GCG_DEV_BUILD           |    BOOL    |                  OFF                  |          compile SCIP and SoPlex located in lib/scip-git and lib/soplex-git          |
|               GMP                |    BOOL    |                  ON                   |                                          -                                           |
|          GMPXX_LIBRARY           |  FILEPATH  | /usr/lib/x86_64-linux-gnu/libgmpxx.so |                                          -                                           |
|         GMP_INCLUDE_DIRS         |    PATH    |     /usr/include/x86_64-linux-gnu     |                                          -                                           |
|           GMP_LIBRARY            |  FILEPATH  |  /usr/lib/x86_64-linux-gnu/libgmp.so  |                                          -                                           |
|               GSL                |    BOOL    |                  ON                   |                                          -                                           |
|              HMETIS              |    BOOL    |                  ON                   |                                          -                                           |
|        HMETIS_EXECUTABLE         |  FILEPATH  |            /opt/bin/hmetis            |                                          -                                           |
|              IPOPT               |    BOOL    |                  ON                   |                                          -                                           |
|            IPOPT_DIR             |    PATH    |                 /usr                  |                                          -                                           |
|              LEGACY              |    BOOL    |                  OFF                  |                                          -                                           |
|               LIBM               |  FILEPATH  |   /usr/lib/x86_64-linux-gnu/libm.so   |                                          -                                           |
|               LPS                |   STRING   |                  spx                  |            What solver to use for LP solving, soplex (spx) or cplex (cpx)            |
|             LPSCHECK             |    BOOL    |                  OFF                  |                                          -                                           |
|                MT                |    BOOL    |                  OFF                  |                                          -                                           |
|           NOBLKBUFMEM            |    BOOL    |                  OFF                  |                                          -                                           |
|             NOBLKMEM             |    BOOL    |                  OFF                  |                                          -                                           |
|             NOBUFMEM             |    BOOL    |                  OFF                  |                                          -                                           |
|              OPENMP              |    BOOL    |                  OFF                  | Use GCG parallelization - combine with PARASCIP (both during compiling SCIP and GCG) |
|             PARASCIP             |    BOOL    |                  OFF                  |                               Use SCIP parallelization                               |
|     POLYSCIP_MAX_NUMBER_OBJS     |   STRING   |                  30                   |                                          -                                           |
|             READLINE             |    BOOL    |                  ON                   |                                          -                                           |
|         SANITIZE_ADDRESS         |    BOOL    |                  OFF                  |                                          -                                           |
|         SANITIZE_MEMORY          |    BOOL    |                  OFF                  |                                          -                                           |
|         SANITIZE_THREAD          |    BOOL    |                  OFF                  |                                          -                                           |
|        SANITIZE_UNDEFINED        |    BOOL    |                  ON                   |                                          -                                           |
|              SHARED              |    BOOL    |                  ON                   |                                          -                                           |
|               SYM                |   STRING   |                 bliss                 |                                          -                                           |
|               TPI                |   STRING   |                 none                  |                                          -                                           |
|              WORHP               |    BOOL    |                  OFF                  |                                          -                                           |
|              ZIMPL               |    BOOL    |                  ON                   |                                          -                                           |
|               ZLIB               |    BOOL    |                  ON                   |                                          -                                           |
|               libm               |  FILEPATH  |   /usr/lib/x86_64-linux-gnu/libm.so   |                                          -                                           |
