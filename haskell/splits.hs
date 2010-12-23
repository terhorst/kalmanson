import System.Environment
import qualified Data.MemoCombinators as Memo
import qualified Data.IntSet as IntSet
import qualified Data.Set as Set
import qualified Data.Map
import List (inits, tails, delete)
import Control.Parallel.Strategies
import Control.Parallel
import Control.Monad.State.Lazy as State

sf s = Set.fromList s
st s = Set.toList s

isf s = IntSet.fromList s
ist s = IntSet.toList s

type Block = IntSet.IntSet
type Split = Set.Set Block
type SplitSystem = Set.Set Split
type SplitSystemMap = Data.Map.Map SplitSystem Bool
type SplitSystemMapS = State SplitSystemMap Bool

permutationsOf [] = [[]]
permutationsOf xs = [x:xs' | x <- xs, xs' <- permutationsOf (delete x xs)]

combinationsOf 0 _ = [[]]
combinationsOf _ [] = []
combinationsOf k (x:xs) = map (x:) (combinationsOf (k-1) xs) ++ combinationsOf k xs

combinations k n = combinationsOf k [1..n]

cartProd (set:sets) = let cp = cartProd sets in [x:xs | x <- set, xs <- cp]
cartProd [] = [[]]

-- Memoization
type StateMap a b = State (Data.Map.Map a b) b
memoizeM :: (Show a, Show b, Ord a) => 
            ((a -> StateMap a b) -> (a -> StateMap a b)) -> (a -> b)
memoizeM t x = evalState (f x) Data.Map.empty where
  g x = do
    y <- t f x  
    m <- get
    put $ Data.Map.insert x y m
    newM <- get
    return y
  f x = get >>= \m -> maybe (g x) return (Data.Map.lookup x m)

-- Begin code

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

-- checkTripleM :: Monad m => (SplitSystem -> m Bool) -> SplitSystem -> m Bool
-- checkTripleM triple = return $ all (weaklyCompatibleTriple) [ sf lst | lst <- cartProd $ map (st) $ st triple ]
-- checkTriple triple = memoizeM checkTripleM triple
ctrp :: SplitSystem -> Bool 
ctrp triple = all (weaklyCompatibleTriple) [ sf lst | lst <- cartProd $ map (st) $ st triple ]
checkTriple :: SplitSystem -> Bool 
checkTriple = ctrp

isWeaklyCompatible :: SplitSystem -> Bool
isWeaklyCompatible ss = all (checkTriple . sf) (combinationsOf 3 $ st ss)

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

ass :: Int -> Int -> Set.Set SplitSystem
ass n 1 = sf [ sf [ sp ] | sp <- st $ allSplits n ]
ass n k = sf [ Set.union a b | a <- st $ allSplitSystems n (k-1), b <- st $ allSplitSystems n 1 ]
allSplitSystems = (Memo.memo2 Memo.integral Memo.integral) ass

css :: Int -> Int -> Set.Set SplitSystem
css n k = sf [ s | (s,circ) <- zip ss isCirc, circ == True ] where
		ss = map (sf) $ combinationsOf k $ st $ allSplits n
		isCirc = (parMap rwhnf) (isCircular) ss
circularSplitSystems = (Memo.memo2 Memo.integral Memo.integral) css

countLifts :: Int -> Int -> Data.Map.Map Int (Set.Set SplitSystem)
countLifts n k = 
	let 
		sskm1 = circularSplitSystems n (k-1)
		ss = circularSplitSystems n k
		subs = [ s1 | s1 <- st sskm1, s2 <- st ss, Set.isSubsetOf s1 s2 ]
	in
		Data.Map.fromListWith (Set.union) [ (occurences subs k, sf [k]) | k <- unique subs ]

unique lst = st $ sf $ lst

occurences [] k = 0
occurences (x:xs) k 
	| x == k = 1 + rec
	| otherwise = 0 + rec
	where rec = occurences xs k

fVector :: Int -> [Int]
fVector n = map (Set.size) [ circularSplitSystems n k | k <- [1 .. n*(n-3) `div` 2] ]

main = 
	do 
		(n:k:_) <- getArgs
		print $ Data.Map.map (Set.size) (countLifts (read n) (read k))
