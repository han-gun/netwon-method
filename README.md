# netwon-method

## Reference 
Netwon Method `(ENG)`
- https://en.wikipedia.org/wiki/Newton%27s_method

뉴턴 방법 `(KOR)`
- https://ko.wikipedia.org/wiki/%EB%89%B4%ED%84%B4_%EB%B0%A9%EB%B2%95

## Newton Method
Newton's method (also known as the Newton–Raphson method), named after Isaac Newton and Joseph Raphson, is a method for finding successively better approximations to the roots (or zeroes) of a real-valued function. It is one example of a [root-finding algorithm](https://en.wikipedia.org/wiki/Root-finding_algorithm).

## The basic iteration Newton method
The method starts with a function `f` defined over the real numbers `x`, the function's derivative `f′`, and an initial guess `x0` for a [root of the function](https://en.wikipedia.org/wiki/Zero_of_a_function) f. If the function satisfies the assumptions made in the derivation of the formula and the initial guess is close, then a better approximation `x1` is

        x1      = x0 - f( x0 ) / f'( x0 )
        
The process is repeated as 

        x_(n+1) = x_(n) - f( x_(n) ) / f'( x_(n) )
        
until a sufficiently accurate value is reached.
