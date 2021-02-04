import parasail
from unittest import TestCase, main

print("parasail is using C lib version", parasail.version())

class Tests(TestCase):

    def test1(self):
        result = parasail.sw("asdf", "asdf", 10, 1, parasail.blosum62)
        self.assertEqual(result.score, 20)

    def test2(self):
        result = parasail.sw("asdf", "asdf", 10, 1, parasail.pam50)
        self.assertEqual(result.score, 27)

    def test3(self):
        matrix = parasail.matrix_create("acgt", 1, -1)
        result = parasail.sw("acgt", "acgt", 10, 1, matrix)
        self.assertEqual(result.score, 4)

    def test4(self):
        print('profile = parasail.profile_create_8("asdf", parasail.blosum62)')
        profile = parasail.profile_create_8("asdf", parasail.blosum62)
        print('result = parasail.sw_striped_profile_8(profile, "asdf", 10, 1)')
        result = parasail.sw_striped_profile_8(profile, "asdf", 10, 1)
        print('self.assertEqual(result.score, 20)')
        self.assertEqual(result.score, 20)


if __name__ == '__main__':
    main(verbosity=3)
