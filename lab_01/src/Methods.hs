module Methods
  (
    getPicard
  , getRungeKutta
  , getEuler
  ) where

type Approximation = Int
type Alpha = Double

rungeKutta :: Double -> Double -> Double -> Alpha -> (Double -> Double -> Double) -> Double
rungeKutta x y h alpha f = y + h * ((1 - alpha) * r1 + alpha * r2)
  where r1 = f x y
        r2 = f (x + h / (2 * alpha)) $ y + h / (2 * alpha) * r1

euler :: Double -> Double -> Double -> (Double -> Double -> Double) -> Double
euler x y h f = y + h * f x y

picard :: Approximation -> Double -> Double
picard 1 x = x ** 3 / 3
picard 2 x = x ** 3 / 3 + x ** 7 / 63
picard 3 x = x ** 3 / 3 + x ** 7 / 63 + 2 * x ** 11 / 2079 + x ** 15 / 59535
picard 4 x = x ** 3 / 3 + x ** 7 / 63 + 2 * x ** 11 / 2079 + x ** 15 / 59535 +
  2 * x ** 15 / 93555 + 2 * x ** 19 / 3393495 + 2 * x ** 19 / 2488563 + 2 * x ** 23 / 86266215 +
  x ** 23 / 99411543 + 2 * x ** 27 / 3341878155 + x ** 31 / 109876902975

getPicard :: Approximation -> Double -> Double -> Double -> [Double]
getPicard approx h x0 xm = map (\x -> picard approx x) [x0, x0 + h..xm]

getEuler :: Double -> Double -> Double -> Double -> (Double -> Double -> Double) -> [Double]
getEuler h x0 xm y0 f = snd $ foldl
  (\acc x ->
    (euler x (fst acc) h f, snd acc ++ [euler x (fst acc) h f])) (y0, []) [x0, x0 + h..xm]

getRungeKutta :: Double -> Alpha -> Double -> Double -> Double -> (Double -> Double -> Double) -> [Double]
getRungeKutta h alpha x0 xm y0 f = snd $ foldl
  (\acc x ->
    (rungeKutta x (fst acc) h alpha f, snd acc ++ [rungeKutta x (fst acc) h alpha f])) (y0, []) [x0, x0 + h..xm]
