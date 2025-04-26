module Main where

import Cosmology
import HMF

args :: [Double]
args = []

main :: IO ()
main = do
  putStrLn "Hello, Haskell!"
  MyLib.someFunc
