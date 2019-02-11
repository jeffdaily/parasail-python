try:
    import parasail
except ImportError:
    import sys, os
    myPath = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, myPath + '/../')
    import parasail

def print_traceback_attributes(traceback):
    print(traceback)
    print(traceback.query)
    print(traceback.comp)
    print(traceback.ref)

def test1():
    matrix = parasail.matrix_create("ACGT", 2, 1)
    result = parasail.sw_trace("ACGT","AcgT",10,1,matrix)
    traceback = result.traceback
    print_traceback_attributes(traceback)

def test21():
    matrix = parasail.matrix_create("ACGTacgt", 2, 1, True)
    result = parasail.sw_trace("ACGT","AcgT",10,1,matrix)
    traceback = result.traceback
    print_traceback_attributes(traceback)

def test22():
    matrix = parasail.matrix_create("ACGTacgt", 2, 1, True)
    result = parasail.sw_trace("ACGT","AcgT",10,1,matrix)
    traceback = result.get_traceback(case_sensitive=True)
    print_traceback_attributes(traceback)

def test3():
    parasail.set_case_sensitive(True)
    matrix = parasail.matrix_create("ACGTacgt", 2, 1)
    result = parasail.sw_trace("ACGT","AcgT",10,1,matrix)
    traceback = result.traceback
    print_traceback_attributes(traceback)

def test4():
    parasail.set_case_sensitive(True)
    matrix = parasail.matrix_create("ACGT", 2, 1)
    result = parasail.sw_trace("ACGT","AcgT",10,1,matrix)
    traceback = result.traceback
    print_traceback_attributes(traceback)

if __name__ == '__main__':
    print("running tests")
    test1()
    test21()
    test22()
    test3()
    test4()
