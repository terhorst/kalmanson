import unittest
import split_system, pq_split_system
from splits import Split
from sage.all import Set

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
        return self.ss_tester(split_system.SplitSystem)

    def test_pq_split_system(self):
        return self.ss_tester(pq_split_system.PQSplitSystem)

if __name__ == '__main__':
    unittest.main()
