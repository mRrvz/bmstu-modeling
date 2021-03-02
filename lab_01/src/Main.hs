import Text.PrettyPrint.Boxes
import Data.List
import Methods

func :: Double -> Double -> Double
func x y = x ** 2 + y ** 2

alpha :: Double
alpha = 1

h :: Double
h = 1e-4

y0 :: Double
y0 = 0

x0 :: Double
x0 = 0

xmax :: Double
xmax = 2.05

everyf :: Int -> [a] -> [a]
everyf n [] = []
everyf n as  = head as : everyf n (drop n as)

printTable :: [[String]] -> IO ()
printTable rows = printBox $ hsep 2 left (map (vcat left . map text) (transpose rows))

formatResults :: [[Double]] -> [[String]]
formatResults = wrapper [[]]
  where
   wrapper res xs
     | length (head xs) == 0 = res
     | otherwise = map show (map head xs) : (wrapper res $ map tail xs)

main :: IO ()
main = do
    let getEvery = 1

    printTable $
      [[ "X"
      ,  "Пикара (1 порядок)"
      ,  "Пикара (2 порядок)"
      ,  "Пикара (3 порядок)"
      ,  "Пикара (4 порядок)"
      ,  "Эйлера"
      ,  "Рунге-Кутта"
      ]] ++ (formatResults [
          everyf getEvery [x0, x0 + h..xmax]
        , everyf getEvery $ getPicard 1 h x0 xmax
        , everyf getEvery $ getPicard 2 h x0 xmax
        , everyf getEvery $ getPicard 3 h x0 xmax
        , everyf getEvery $ getPicard 4 h x0 xmax
        , everyf getEvery $ getEuler h x0 xmax y0 func
        , everyf getEvery $ getRungeKutta h alpha x0 xmax y0 func
      ])
