import unittest, random
import split_system, pq_split_system
from splits import Split
from sage.all import Set
import kalmanson

class TestSplits(unittest.TestCase):
    def setUp(self):
        X = range(1,6)
        self.circular_ss = Set(Split(X,map(int, a))
                for a in [[1,2],[2,3],[3,4],[4,5],[5,1]])
        self.non_circular_ss = self.circular_ss + Set([Split(X,[1,4])])

    def ss_tester(self, SsCls):
        ss = SsCls(self.circular_ss)
        self.assertTrue(ss.is_circular())
        ss = SsCls(self.non_circular_ss)
        self.assertTrue(not ss.is_circular())

    def test_split_system(self):
        sscls = split_system.SplitSystem
        self.ss_tester(sscls)
        self.kalmanson_tester(sscls)

    def test_pq_split_system(self):
        sscls = pq_split_system.PQSplitSystem
        self.ss_tester(sscls)
        self.kalmanson_tester(sscls)

    def kalmanson_tester(self, SsCls):
        n = 6
        kn = kalmanson.kalmanson_fan(6)
        X = range(n)
        cones = random.sample([c for cn in kn.cones()[3:] for c in cn], 500)
        for c in cones:
            ss = SsCls( [ Split(X,s[0]) for s in kalmanson.cone_to_split_system(n,c) ] )
            self.assertTrue(ss.is_circular())

if __name__ == '__main__':
    unittest.main()
