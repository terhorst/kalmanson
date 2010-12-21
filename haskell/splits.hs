import qualified Data.Set as Set
import List (inits, tails, delete)

sf s = Set.fromList s
st s = Set.toList s

type Block = Set.Set Int
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
makeSplit n block = sf [a, Set.difference (sf [1..n]) a] where a = sf block 

weaklyCompatibleTriple :: Set.Set Block -> Bool
weaklyCompatibleTriple ss = (Set.null $ foldr1 (Set.intersection) ss') ||
		(any (Set.null . foldl1 Set.difference) $ permutationsOf ss')
		where ss' = st ss

isWeaklyCompatible :: SplitSystem -> Bool
isWeaklyCompatible ss =
	all (weaklyCompatibleTriple) [ sf lst |
		triple <- combinationsOf 3 $ st ss,
		lst <- cartProd $ map (st) triple ]

isCircular :: SplitSystem -> Bool
isCircular ss = isWeaklyCompatible ss' where
	ss' = Set.union ss $ sf [ mergeSplits s1 s2 |
		(s1:s2:_) <- combinationsOf 2 $ st ss,
		(Set.size $ Set.intersection (splitB s1) (splitB s2)) > 0 ]

splitA :: Split -> Block
splitA sp 	| Set.member 1 a = a
			| otherwise = b
	where (a:b:_) = Set.toList sp

splitB :: Split -> Block
splitB sp = a where (a:_) = Set.toList $ Set.delete (splitA sp) sp

joinSplitSystems :: SplitSystem -> SplitSystem -> SplitSystem
joinSplitSystems ss1 ss2 = Set.union ss1 ss2

mergeSplits :: Split -> Split -> Split
mergeSplits s1 s2 = sf [ Set.intersection a1 a2, Set.union b1 b2 ]
	where 	(a1:a2:b1:b2:_) = [splitA s1, splitA s2, splitB s1, splitB s2]

allSplits :: Int -> Set.Set Split
allSplits n = sf [ sf [sp, Set.difference x sp] |
	k <- [2..n-2],
	sp <- map sf $ combinations k n ]
	where x = sf [1 .. n]

circularSplitSystems :: Int -> Int -> Set.Set SplitSystem
-- circularSplitSystems n k 	| k == 1 = sf [ sf [sp] | sp <- st $ allSplits n ]
--	 							| otherwise = Set.filter (isCircular)
--									(sf $ map (sf) $ combinationsOf k $ st $ allSplits n)
circularSplitSystems n k = Set.filter (isCircular) 
	(sf $ map (sf) $ combinationsOf k $ st $ allSplits n)

fVector :: Int -> [Int]
fVector n = map (Set.size) [ circularSplitSystems n k | k <- [1 .. n*(n-3) `div` 2] ]

main = print $ fVector 6
