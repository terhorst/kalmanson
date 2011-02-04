import unittest, formulas, alphabet, consecutive_ones

# The range of split system sizes to test over. 
MIN_N = 6
MAX_N = 7

class TestFormulas(unittest.TestCase):

    def test_formulas(self):
        for n in range(MIN_N, MAX_N + 1):
            d = alphabet.all_F(alphabet.__get_mats(n))
            for k,v in d.iteritems():
                try:
                    f = formulas.getattr("s%s" % "".join(map(str, k)))
                except AttributeError:
                    self.assertEqual(len(v), 0)

                self.assertEqual(len(v), f(n))

    def test_fvector3(self):
        for n in range(MIN_N, MAX_N + 1):
            v = fvector(n,3)
            self.assertEqual
