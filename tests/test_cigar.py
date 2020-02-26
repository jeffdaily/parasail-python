import parasail
from unittest import TestCase, main


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


class Tests(TestCase):

    def test0(self):
        result = parasail.sw("asdf", "asdf", 10, 1, parasail.blosum62)
        with self.assertRaises(AttributeError):
            print_cigar_attributes(result.cigar)

    def test1(self):
        result = parasail.sw("asdf", "asdf", 10, 1, parasail.blosum62)
        with self.assertRaises(AttributeError):
            print_traceback_attributes(result.traceback)

    def test2(self):
        result = parasail.sw_trace("asdf", "asdf", 10, 1, parasail.blosum62)
        print_cigar_attributes(result.cigar)

    def test3(self):
        result = parasail.sw_trace("asdf", "asdf", 10, 1, parasail.blosum62)
        print_traceback_attributes(result.traceback)


if __name__ == '__main__':
    main()
