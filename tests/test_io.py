import os
import parasail
from unittest import TestCase, main
from tempfile import NamedTemporaryFile


def create_input_file(filename):
    with open(filename, 'w') as fp:
        fp.write('''>AF0017_1 COG1250 # Protein_GI_number: 11497638 # Func_class: I Lipid transport and metabolism  # Function: 3-hydroxyacyl-CoA dehydrogenase # Organism: Archaeoglobus fulgidus
MMVLEIRNVAVIGAGSMGHAIAEVVAIHGFNVKLMDVSEDQLKRAMEKIEEGLRKSYERGYISEDPEKVLKRIEATADLIEVAKDADLVIEAIPEIFDLKKKVFSEIEQYCPDHTIFATNTSSLSITKLAEATKRPEKFIGMHFFNPPKILKLLEIVWGEKTSEETIRIVEDFARKIDRIIIHVRKDVPGFIVNRIFVTMSNEASWAVEMGEGTIEEIDSAVKYRLGLPMGLFELHDVLGGGSVDVSYHVLEYYRQTLGESYRPSPLFERLFKAGHYGKKTGKGFYDWSEGKTNEVPLRAGANFDLLRLVAPAVNEAAWLIEKGVASAEEIDLAVLHGLNYPRGLLRMADDFGIDSIVKKLNELYEKYNGEERYKVNPVLQKMVEEGKLGRTTGEGFYKYGD
>AF0017_2_COG1024_#_Protein_GI_number: 11497638 # Func_class: I Lipid transport and metabolism  # Function: Enoyl-CoA hydratase/carnithine racemase # Organism: Archaeoglobus fulgidus
GNYEFVKVEKEGKVGVLKLNRPRRANALNPTFLKEVEDALDLLERDEEVRAIVIAGEGKNFCAGADIAMFASGRPEMVTEFSQLGHKVFRKIEMLSKPVIAAIHGAAVGGGFELAMACDLRVMSERAFLGLPELNLGIIPGWGGTQRLAYYVGVSKLKEVIMLKRNIKPEEAKNLGLVAEVFPQERFWDEVMKLAREVAELPPLAVKYLKKVIALGTMPALETGNLAESEAGAVIALTDDVAEGIQAFNYRRKPNFRGR
''')


class Tests(TestCase):

    def setUp(self):
        self.tmpfile = NamedTemporaryFile()
        create_input_file(self.tmpfile.name)

    def tearDown(self):
        os.unlink(self.tmpfile.name)
        os.path.exists(self.tempfile.name)

    def work(self):
        sequences = parasail.sequences_from_file(self.tmpfile)
        print(len(sequences))
        print(len(sequences[0]))
        print(sequences[0])
        print(len(sequences[-1]))
        print(sequences[-1])
        print(sequences[0][0])
        with self.assertRaises(TypeError):
            print(sequences['asdf'])
        with self.assertRaises(TypeError):
            print(sequences[0]['asdf'])
        with self.assertRaises(IndexError):
            print(sequences[1000000])
        with self.assertRaises(IndexError):
            print(sequences[-1000000])
        with self.assertRaises(IndexError):
            print(sequences[0][1000000])
        with self.assertRaises(IndexError):
            print(sequences[0][-1000000])
        print("name:    '{}'".format(sequences[0].name))
        print("comment: '{}'".format(sequences[0].comment))
        print("seq:     '{}'".format(sequences[0].seq))
        print("qual:    '{}'".format(sequences[0].qual))
        print("characters: {}".format(sequences.characters))
        print("shortest:   {}".format(sequences.shortest))
        print("longest:    {}".format(sequences.longest))
        print("mean:       {}".format(sequences.mean))
        print("stddev:     {}".format(sequences.stddev))
        result = parasail.sw(str(sequences[0]), str(sequences[1]), 10, 1, parasail.blosum62)
        print(result.score)
        result = parasail.sw(sequences[0], sequences[1], 10, 1, parasail.blosum62)
        print(result.score)


if __name__ == '__main__':
    main()
