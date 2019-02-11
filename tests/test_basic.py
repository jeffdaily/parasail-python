try:
    import parasail
except ImportError:
    import sys, os
    myPath = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, myPath + '/../')
    import parasail

def test1():
    result = parasail.sw("asdf","asdf",10,1,parasail.blosum62)
    assert(result.score == 20)
    del result

def test2():
    result = parasail.sw("asdf","asdf",10,1,parasail.pam50)
    assert(result.score == 27)
    del result

def test3():
    matrix = parasail.matrix_create("acgt", 1, -1)
    result = parasail.sw("acgt","acgt",10,1,matrix)
    assert(result.score == 4)
    del result
    del matrix

def test4():
    profile = parasail.profile_create_8("asdf", parasail.blosum62)
    result = parasail.sw_striped_profile_8(profile,"asdf",10,1)
    assert(result.score == 20)
    del result
    del profile

if __name__ == '__main__':
    print("running tests")
    test1()
    test2()
    test3()
    test4()
