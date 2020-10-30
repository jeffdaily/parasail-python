import parasail
from unittest import TestCase, main, expectedFailure


class Tests(TestCase):

    def test1(self):
        p = parasail.ssw_init("asdf", parasail.blosum62, 1)
        r = parasail.ssw_profile(p, "asdf", 10, 1)

        print(p.s1)
        print(p.s1Len)
        print(r.cigarLen)
        print(r.cigar[0])

        r = parasail.sw_trace("asdf", "asdf", 10, 1, parasail.blosum62)
        c = r.cigar
        print(c.len)
        print(c.seq[0])
        print(c.decode)

        p = parasail.profile_create_8("asdf", parasail.blosum62)
        r = parasail.sw_trace_striped_profile_8(p, "asdf", 10, 1)
        c = r.cigar
        print(c.len)
        print(c.seq[0])

        r = parasail.sw_trace("asdf", "asdf", 10, 1, parasail.blosum62)
        print(r.query)
        print(r.ref)


if __name__ == '__main__':
    main()
