import parasail
from unittest import TestCase, main


class Tests(TestCase):

    def print_traceback_attributes(self, traceback):
        print(traceback)
        print(traceback.query)
        print(traceback.comp)
        print(traceback.ref)

    def test1(self):
        matrix = parasail.matrix_create("ACGT",  2,  1)
        result = parasail.sw_trace("ACGT", "AcgT", 10, 1, matrix)
        traceback = result.traceback
        self.print_traceback_attributes(traceback)

    def test21(self):
        matrix = parasail.matrix_create("ACGTacgt",  2,  1,  True)
        result = parasail.sw_trace("ACGT", "AcgT", 10, 1, matrix)
        traceback = result.traceback
        self.print_traceback_attributes(traceback)

    def test22(self):
        matrix = parasail.matrix_create("ACGTacgt",  2,  1,  True)
        result = parasail.sw_trace("ACGT", "AcgT", 10, 1, matrix)
        traceback = result.get_traceback(case_sensitive=True)
        self.print_traceback_attributes(traceback)

    def test3(self):
        parasail.set_case_sensitive(True)
        matrix = parasail.matrix_create("ACGTacgt",  2,  1)
        result = parasail.sw_trace("ACGT", "AcgT", 10, 1, matrix)
        traceback = result.traceback
        self.print_traceback_attributes(traceback)

    def test4(self):
        parasail.set_case_sensitive(True)
        matrix = parasail.matrix_create("ACGT",  2,  1)
        result = parasail.sw_trace("ACGT", "AcgT", 10, 1, matrix)
        traceback = result.traceback
        self.print_traceback_attributes(traceback)


if __name__ == '__main__':
    main()
