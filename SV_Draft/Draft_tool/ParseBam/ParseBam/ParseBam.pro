TEMPLATE = app
#CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
INCLUDEPATH += ../../../../ShareLibrary/
LIBS += -lz
CONFIG += thread
CONFIG += c++11

HEADERS += \
    clsparsebam.h \
    ../../../../ShareLibrary/clsbasealgorithm.h \
    ../../../../ShareLibrary/clsbwa.h \
    ../../../../ShareLibrary/clsfastareader.h \
    ../../../../ShareLibrary/clsfastqreader.h \
    ../../../../ShareLibrary/clsreadconfigini.h \
    clsconfig.h \
    ../../../../ShareLibrary/clsvcf1000genome.h \
    ../../../../ShareLibrary/clsmuscle.h \
    clscomparison.h \
    clsdebug.h \
    ../../../../ShareLibrary/clsvelvet.h \
    ../../../../ShareLibrary/LocalAssembly/local_alignment.h \
    ../../../../ShareLibrary/LocalAssembly/stdaln.h \
    ../../../../ShareLibrary/clsblast.h \
    ../../../../ShareLibrary/clskmeralgorithm.h \
    clsdrawimage.h

SOURCES += \
    clsparsebam.cpp \
    main.cpp \
    ../../../../ShareLibrary/clsbasealgorithm.cpp \
    ../../../../ShareLibrary/clsbwa.cpp \
    ../../../../ShareLibrary/clsfastareader.cpp \
    ../../../../ShareLibrary/clsfastqreader.cpp \
    ../../../../ShareLibrary/clsreadconfigini.cpp \
    clsconfig.cpp \
    ../../../../ShareLibrary/clsvcf1000genome.cpp \
    ../../../../ShareLibrary/clsmuscle.cpp \
    clscomparison.cpp \
    clsdebug.cpp \
    ../../../../ShareLibrary/clsvelvet.cpp \
    ../../../../ShareLibrary/LocalAssembly/local_alignment.cpp \
    ../../../../ShareLibrary/clsblast.cpp \
    ../../../../ShareLibrary/clskmeralgorithm.cpp \
    clsdrawimage.cpp

unix:!macx: LIBS += -L$$PWD/../../../../bamtools/install/lib/ -lbamtools

INCLUDEPATH += $$PWD/../../../../bamtools/install/include/bamtools
DEPENDPATH += $$PWD/../../../../bamtools/install/include/bamtools

DISTFILES += \
    Draft_Solution

unix:!macx: LIBS += -L$$PWD/../../../../../Software/opencv/install/lib/ -lopencv_core

unix:!macx: LIBS += -L$$PWD/../../../../../Software/opencv/install/lib/ -lopencv_highgui

unix:!macx: LIBS += -L$$PWD/../../../../../Software/opencv/install/lib/ -lopencv_imgproc

unix:!macx: LIBS += -L$$PWD/../../../../../Software/opencv/install/lib/ -lopencv_imgcodecs

INCLUDEPATH += $$PWD/../../../../../Software/opencv/install/include
DEPENDPATH += $$PWD/../../../../../Software/opencv/install/include
