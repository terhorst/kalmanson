import Data.Bits
import System.Environment
import qualified Data.MemoCombinators as Memo
import qualified Data.IntSet as IntSet
import qualified Data.Set as Set
import qualified Data.Map
import List (inits, tails, delete, (\\))
import Control.Parallel.Strategies
import Control.Parallel

sf s = IntSet.fromList s
st s = IntSet.toList s

type Split = Int
type SplitSystem = IntSet.IntSet

permutationsOf [] = [[]]
permutationsOf xs = [x:xs' | x <- xs, xs' <- permutationsOf (delete x xs)]

combinationsOf 0 _ = [[]]
combinationsOf _ [] = []
combinationsOf k (x:xs) = map (x:) (combinationsOf (k-1) xs) ++ combinationsOf k xs

combinations k n = combinationsOf k [1..n]

combinationsOfSet k s = 1

cartProd (set:sets) = let cp = cartProd sets in [x:xs | x <- set, xs <- cp]
cartProd [] = [[]]

memoSs = Memo.wrap IntSet.fromAscList IntSet.toAscList (Memo.list Memo.integral)

makeSplit :: Int -> [Int] -> Split
makeSplit n block = let
	a | elem 1 block = block
	  | otherwise = [1..n] \\ block 
	in
	foldr1 (.|.) (map (bit . (flip subtract) n) a)

fromSplit :: Int -> Int -> [Int]
fromSplit n sp = [ i | i <- [1..n], 2^(n-i) .&. sp == 2^(n-i) ]

maskedComplement :: Int -> Int -> Int
maskedComplement n x = (complement x) .&. (2^n - 1)

wct :: Int -> SplitSystem -> Bool
wct n ss = let
			ss' = st ss
			comp = maskedComplement n
			c = [ concat [ss' \\ combo, map (comp) combo] | combo <- combinationsOf 2 ss' ]
			myfoldr = foldr1 (.&.) :: [Int] -> Int
		in
		(myfoldr ss' .&. (2^n - 1) == 0) || (any ((==0) .  myfoldr) c)
weaklyCompatibleTriple = Memo.memo2 Memo.integral memoSs wct

ctrp :: Int -> SplitSystem -> Bool 
ctrp n triple = all ((weaklyCompatibleTriple n). sf) $ 
	cartProd (map (\x -> [x, complement x]) $ st triple)
checkTriple = Memo.memo2 Memo.integral memoSs ctrp

isWeaklyCompatible :: Int -> SplitSystem -> Bool
isWeaklyCompatible n ss = all ((checkTriple n). sf) (combinationsOf 3 $ st ss)

isCircular :: Int -> SplitSystem -> Bool
isCircular n ss = isWeaklyCompatible n ss' where
	comp = maskedComplement n
	ss' = IntSet.union ss $ sf [ mergeSplits a b | 
		(s1:s2:_) <- combinationsOf 2 $ st ss,
		(a:b:_) <- cartProd [ [s1, comp s1], [s2, comp s2] ],
		comp a .&. comp b > 0 ]

joinSplitSystems :: SplitSystem -> SplitSystem -> SplitSystem
joinSplitSystems ss1 ss2 = IntSet.union ss1 ss2

mergeSplits :: Split -> Split -> Split
mergeSplits = (.&.)

allSplits :: Int -> [Split]
allSplits n = unique [ makeSplit n block | k <- [2..n-2], block <- combinations k n ]
-- allSplits n = [2^j + 1 .. 2^(j+1) - 4] where j = n-1

ass :: Int -> Int -> Set.Set SplitSystem
ass n 1 = Set.fromList [ IntSet.singleton sp | sp <- allSplits n ]
ass n k = Set.fromList $ map (foldr1 IntSet.union) (combinationsOf k (Set.toList $ allSplitSystems n 1))
allSplitSystems = (Memo.memo2 Memo.integral Memo.integral) ass

css :: Int -> Int -> Set.Set SplitSystem
css n k = Set.fromList [ s | (s,circ) <- zip ss isCirc, circ == True ] where
		ss = Set.toList $ possCircularSplitSystems n k
		isCirc = (parMap rwhnf) (isCircular n) ss
circularSplitSystems = (Memo.memo2 Memo.integral Memo.integral) css

-- Return all the splits systems in K_n of length k which are invalid
possCircularSplitSystems :: Int -> Int -> Set.Set SplitSystem
possCircularSplitSystems n k = Set.difference (allSplitSystems n k) $	
	Set.fromList [ joinSplitSystems a b | j <- [1 .. (k-1)], 
						a <- Set.toList $ allSplitSystems n (k - j),
						b <- Set.toList $ 
								Set.difference (allSplitSystems n j) (circularSplitSystems n j) ]

-- countLifts :: Int -> Int -> Data.Map.Map Int (Set.Set SplitSystem)
-- countLifts n k = 
	-- let 
		-- sskm1 = circularSplitSystems n (k-1)
		-- ss = circularSplitSystems n k
		-- subs = [ s1 | s1 <- st sskm1, s2 <- st ss, Set.isSubsetOf s1 s2 ]
	-- in
		-- Data.Map.fromListWith (Set.union) [ (occurences subs k, sf [k]) | k <- unique subs ]

unique lst = st $ sf $ lst

occurences [] k = 0
occurences (x:xs) k 
	| x == k = 1 + rec
	| otherwise = 0 + rec
	where rec = occurences xs k

fVector :: Int -> [Int]
fVector n = map (Set.size) [ circularSplitSystems n k | k <- [1 .. n*(n-3) `div` 2] ]

partialfVector :: Int -> Int -> [Int]
partialfVector n k = map (Set.size) [ circularSplitSystems n j | j <- [1 .. k] ]

main = do
	(n:k:_) <- getArgs
	print $ partialfVector (read n) (read k)
