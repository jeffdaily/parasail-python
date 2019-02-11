try:
    import parasail
except ImportError:
    import sys, os
    myPath = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, myPath + '/../')
    import parasail

def print_cigar_attributes(cigar):
    print(cigar)
    print(cigar.seq)
    print(cigar.len)
    print(cigar.beg_query)
    print(cigar.beg_ref)
    print(cigar.decode)

def print_traceback_attributes(traceback):
    print(traceback)
    print(traceback.query)
    print(traceback.comp)
    print(traceback.ref)

def test0():
    result = parasail.sw("asdf","asdf",10,1,parasail.blosum62)
    try:
        cigar = result.cigar
        print_cigar_attributes(cigar)
    except AttributeError:
        pass
    else:
        raise Error("AttributeError expected, Failure")
def test1():
    result = parasail.sw("asdf","asdf",10,1,parasail.blosum62)
    try:
        traceback = result.traceback
        print_traceback_attributes(traceback)
    except AttributeError:
        pass
    else:
        raise Error("AttributeError expected, Failure")

def test2():
    result = parasail.sw_trace("asdf","asdf",10,1,parasail.blosum62)
    cigar = result.cigar
    print_cigar_attributes(cigar)

def test3():
    result = parasail.sw_trace("asdf","asdf",10,1,parasail.blosum62)
    traceback = result.traceback
    print_traceback_attributes(traceback)

if __name__ == '__main__':
    print("running tests")
    test0()
    test1()
    test2()
    test3()
