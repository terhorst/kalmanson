import qualified Data.MemoCombinators as Memo
import qualified Data.IntSet as IntSet
import qualified Data.Set as Set
import List (inits, tails, delete)
import Control.Parallel.Strategies
import Control.Parallel

sf s = Set.fromList s
st s = Set.toList s

isf s = IntSet.fromList s
ist s = IntSet.toList s

type Block = IntSet.IntSet
type Split = Set.Set Block
type SplitSystem = Set.Set Split

permutationsOf [] = [[]]
permutationsOf xs = [x:xs' | x <- xs, xs' <- permutationsOf (delete x xs)]

combinationsOf 0 _ = [[]]
combinationsOf _ [] = []
combinationsOf k (x:xs) = map (x:) (combinationsOf (k-1) xs) ++ combinationsOf k xs

combinations k n = combinationsOf k [1..n]

cartProd (set:sets) = let cp = cartProd sets in [x:xs | x <- set, xs <- cp]
cartProd [] = [[]]

makeSplit :: Int -> [Int] -> Split
makeSplit n block = sf [a, IntSet.difference (isf [1..n]) a] where a = isf block 

wct :: Set.Set Block -> Bool
wct ss = (IntSet.null $ foldr1 (IntSet.intersection) ss') ||
		(any (IntSet.null . foldl1 IntSet.difference) $ permutationsOf ss')
		where ss' = st ss
weaklyCompatibleTriple = (Memo.wrap toTriple fromTriple memoString) wct where
	memoString = Memo.list Memo.char
	toTriple str = (read str) :: Set.Set Block
	fromTriple trip = show trip

ctrp :: [Split] -> Bool
ctrp triple = all (weaklyCompatibleTriple) [ sf lst | lst <- cartProd $ map (st) triple ]
checkTriple = (Memo.wrap toTriple fromTriple memoString) ctrp where
	memoString = Memo.list Memo.char
	toTriple str = (read str) :: [Split]
	fromTriple trip = show trip

isWeaklyCompatible :: SplitSystem -> Bool
isWeaklyCompatible ss = all (checkTriple) (combinationsOf 3 $ st ss)

isCircular :: SplitSystem -> Bool
isCircular ss = isWeaklyCompatible ss' where
	ss' = Set.union ss $ sf [ mergeSplits s1 s2 |
		(s1:s2:_) <- combinationsOf 2 $ st ss,
		(IntSet.size $ IntSet.intersection (splitB s1) (splitB s2)) > 0 ]

splitA :: Split -> Block
splitA sp 	| IntSet.member 1 a = a
			| otherwise = b
	where (a:b:_) = Set.toList sp

splitB :: Split -> Block
splitB sp = a where (a:_) = Set.toList $ Set.delete (splitA sp) sp

joinSplitSystems :: SplitSystem -> SplitSystem -> SplitSystem
joinSplitSystems ss1 ss2 = Set.union ss1 ss2

mergeSplits :: Split -> Split -> Split
mergeSplits s1 s2 = sf [ IntSet.intersection a1 a2, IntSet.union b1 b2 ]
	where 	(a1:a2:b1:b2:_) = [splitA s1, splitA s2, splitB s1, splitB s2]

allSplits :: Int -> Set.Set Split
allSplits n = sf [ sf [sp, IntSet.difference x sp] |
	k <- [2..n-2],
	sp <- map isf $ combinations k n ]
	where x = isf [1 .. n]

css :: Int -> Int -> Set.Set SplitSystem
css n k = sf [ s | (s,circ) <- zip ss isCirc, circ == True ] where
		ss = map (sf) $ combinationsOf k $ st $ allSplits n
		isCirc = (parMap rwhnf) (isCircular) ss
circularSplitSystems = (Memo.memo2 Memo.integral Memo.integral) css

fVector :: Int -> [Int]
fVector n = map (Set.size) [ circularSplitSystems n k | k <- [1 .. n*(n-3) `div` 2] ]

main = print $ fVector 5
