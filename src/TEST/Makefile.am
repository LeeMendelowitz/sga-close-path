bin_PROGRAMS = testDummyThread

testDummyThread_CPPFLAGS = \
    -I$(top_srcdir)/Util \
    -I$(top_srcdir)/Bigraph \
    -I$(top_srcdir)/SuffixTools \
    -I$(top_srcdir)/StringGraph \
    -I$(top_srcdir)/Concurrency \
    -I$(top_srcdir)/Algorithm \
    -I$(top_srcdir)/SQG \
    -I$(top_srcdir)/Scaffold \
    -I$(top_srcdir)/GraphDiff \
    -I$(top_srcdir)/GraphSearch \
    -I$(top_srcdir)/Thirdparty

testDummyThread_LDADD = \
    $(top_builddir)/Scaffold/libscaffold.a \
    $(top_builddir)/GraphDiff/libgraphdiff.a \
    $(top_builddir)/StringGraph/libstringgraph.a \
    $(top_builddir)/Concurrency/libconcurrency.a \
    $(top_builddir)/Algorithm/libalgorithm.a \
    $(top_builddir)/SuffixTools/libsuffixtools.a \
    $(top_builddir)/Bigraph/libbigraph.a \
    $(top_builddir)/Util/libutil.a \
    $(top_builddir)/SQG/libsqg.a \
    $(top_builddir)/GraphSearch/libgraphsearch.a \
    $(top_builddir)/Thirdparty/libthirdparty.a 

testDummyThread_LDFLAGS = -pthread

testDummyThread_SOURCES = ThreadDummyTest.cpp \
                          ThreadDummy.h
