execute_process(COMMAND python -c "import ${LIBRARY}.version ; print(${LIBRARY}.version.__file__)" OUTPUT_VARIABLE VERSION_PATH)
STRING(STRIP ${VERSION_PATH} VERSION_PATH)
file(WRITE ${VERSION_PATH} ${GIT_VERSION})
