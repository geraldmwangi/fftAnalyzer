TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt
LIBS +=-lmgl-qt -lmgl -lmgl-wnd -lsndfile  -lfftw3f -lm -lsamplerate
SOURCES += main.cpp \
    waveform.cpp \
    waveformplot.cpp \
    fftanalysis.cpp \
    multiRansac.cpp \
    frequencysieve.cpp \
    autocorrelator.cpp

HEADERS += \
    waveform.h \
    waveformplot.h \
    fftanalysis.h \
    multiRansac.hxx \
    frequencysieve.h \
    autocorrelator.h
